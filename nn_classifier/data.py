"""
Data loading for the homology classifier.

Splits are *by query* so the same query never appears in both train and test —
otherwise the model learns query identity (composition is fixed per query).

FAM + SFAM = positive (homolog).
CROSS      = negative.
FOLD       = held out entirely, surfaced separately for eval.
"""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Sequence
import glob
import json
import os
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import torch
from torch.utils.data import TensorDataset

# Default feature set — drops length-coupled raw features (log_qlen, log_tlen,
# raw bits — bits is kept separately for the baseline only).
_AA_LETTERS  = "ACDEFGHIKLMNPQRSTVWY"
_TDI_LETTERS = "ACDEFGHIKLMNPQRSTVWY"
DEFAULT_FEATURE_COLS: tuple = (
    # Scalar features (9): raw + log of the two scores, plus alignment quality.
    "evalue", "log_evalue", "bits", "log_bits",
    "lddt", "log_alnlen", "gap_frac", "fident", "qtmscore",
) + tuple(f"q_aln_comp_b{i:02d}" for i in range(16)) \
  + tuple(f"t_aln_comp_b{i:02d}" for i in range(16)) \
  + tuple(f"q_aln_aa_{c}"  for c in _AA_LETTERS) \
  + tuple(f"t_aln_aa_{c}"  for c in _AA_LETTERS) \
  + tuple(f"q_aln_3di_{c}" for c in _TDI_LETTERS) \
  + tuple(f"t_aln_3di_{c}" for c in _TDI_LETTERS)


@dataclass
class Splits:
    train: TensorDataset
    val:   TensorDataset
    test:  TensorDataset
    fold_holdout: TensorDataset           # FOLD pairs, never used for training
    feature_cols: list
    feat_mean: np.ndarray
    feat_std:  np.ndarray
    pos_weight: float                     # #neg / #pos in train


def _split_queries(queries: np.ndarray, val_pct: float, test_pct: float,
                   seed: int = 42):
    rng = np.random.default_rng(seed)
    uniq = np.unique(queries)
    rng.shuffle(uniq)
    n_test = int(len(uniq) * test_pct / 100)
    n_val  = int(len(uniq) * val_pct  / 100)
    test = set(uniq[:n_test])
    val  = set(uniq[n_test:n_test + n_val])
    train = set(uniq[n_test + n_val:])
    return train, val, test


def _split_by_fold(queries: np.ndarray, scop_lookup_path: str,
                   val_pct: float, test_pct: float, seed: int = 42):
    """Fold-level (SCOP first 2 levels) train/val/test split, balanced so
    the *query-count* shares match val_pct / test_pct rather than the
    fold-count shares.  Fold variance in size (some folds have 1 domain,
    others 400+) means a naive fold-count split can drift by ~10 pp.

    Greedy packing: shuffle folds, pack into test until query share ≥
    test_pct, then into val until query share ≥ val_pct, rest into train.
    No fold is ever split across sets — leakage is still eliminated.
    """
    scop = pd.read_csv(scop_lookup_path, sep="\t", header=None,
                       names=["domain", "scop"])
    scop["fold"] = scop["scop"].str.split(".").str[:2].str.join(".")
    d2fold = dict(zip(scop["domain"], scop["fold"]))

    qs = np.unique(queries)
    folds = {}
    for q in qs:
        f = d2fold.get(q)
        if f is None:
            continue
        folds.setdefault(f, []).append(q)

    fold_keys = sorted(folds.keys())
    rng = np.random.default_rng(seed)
    rng.shuffle(fold_keys)

    total = sum(len(v) for v in folds.values())
    test_target = total * test_pct / 100
    val_target  = total * val_pct  / 100

    test_folds, val_folds, train_folds = [], [], []
    cum_test = 0
    cum_val = 0
    for f in fold_keys:
        n = len(folds[f])
        if cum_test < test_target:
            test_folds.append(f); cum_test += n
        elif cum_val < val_target:
            val_folds.append(f);  cum_val += n
        else:
            train_folds.append(f)

    test  = set(q for f in test_folds  for q in folds[f])
    val   = set(q for f in val_folds   for q in folds[f])
    train = set(q for f in train_folds for q in folds[f])
    return train, val, test


def load_splits(parquet_dir: str,
                val_pct: float = 10,
                test_pct: float = 10,
                seed: int = 42,
                norm_stats_out: str | None = None,
                feature_cols: Sequence[str] = DEFAULT_FEATURE_COLS,
                cross_top_k_per_query: int = 300,
                cross_random_tail_frac: float = 0.005,
                cross_tail_per_query: int = 0,
                extra_hard_negatives_tsv: str | None = None,
                split_by_fold: bool = False,
                scop_lookup_path: str | None = None,
                exclude_folds: Sequence[str] | None = None,
                explicit_split_json: str | None = None) -> Splits:
    shards = sorted(glob.glob(os.path.join(parquet_dir, "shard_*.parquet")))
    if not shards:
        raise FileNotFoundError(f"no shards in {parquet_dir}")

    feature_cols = list(feature_cols)
    # `target` is needed only when extra-hard-negatives filtering is requested,
    # but it's cheap to always load and keeps the code path uniform.
    needed = ["query", "target", "rel"] + feature_cols
    rng = np.random.default_rng(seed)

    # ----- Pass 1: per-query top-K CROSS threshold ------------------------
    # Random subsampling would discard exactly the rows the classifier needs:
    # the high-scoring CROSS pairs are the FALSE POSITIVES that define the
    # discriminative boundary. We instead keep, per query, the K CROSS pairs
    # with the largest `bits` — the *hardest* negatives for that query —
    # plus a small random tail of the rest for diversity.
    t_p1 = time.time()
    print(f"Pass 1: per-query CROSS thresholds "
          f"(top_k={cross_top_k_per_query}, tail_frac={cross_random_tail_frac:.4f}) …",
          flush=True)
    cross_qb = []
    for i, sh in enumerate(shards):
        d = pq.read_table(sh, columns=["query", "rel", "bits"]).to_pandas()
        cross_qb.append(d.loc[d["rel"] == "CROSS", ["query", "bits"]])
    qb = pd.concat(cross_qb, ignore_index=True)
    del cross_qb
    # K-th highest bits per query; rows above that threshold are "hard".
    # cross_top_k_per_query <= 0 disables bits-based mining (threshold +inf).
    def _kth_thresh(g):
        if cross_top_k_per_query <= 0:
            return float("inf")
        if g.size >= cross_top_k_per_query:
            return float(np.partition(g.values, -cross_top_k_per_query)[-cross_top_k_per_query])
        return float("-inf")
    thresholds = qb.groupby("query")["bits"].agg(_kth_thresh).to_dict()
    print(f"  threshold computed for {len(thresholds):,} queries  "
          f"(median={np.median(list(thresholds.values())):.1f}, "
          f"in {time.time() - t_p1:.1f}s)", flush=True)
    # Per-query random tail: instead of a flat global per-row coin, keep ~N
    # random CROSS per query (Bernoulli with p = N / cross_count[q]). This
    # balances the random negatives across queries to match the per-query
    # averaging of the ROC1/ROC5 eval metrics. 0 = flat global tail (legacy).
    tail_prob_map = None
    if cross_tail_per_query > 0:
        cross_count = qb.groupby("query").size()
        tail_prob_map = (cross_tail_per_query / cross_count).clip(upper=1.0).to_dict()
        print(f"  per-query tail: target {cross_tail_per_query}/query over "
              f"{len(tail_prob_map):,} queries (median p="
              f"{np.median(list(tail_prob_map.values())):.4f})", flush=True)
    del qb

    # ----- Optional: model-confidence hard negatives ---------------------
    # Round-2 mining: pairs the previous model is wrongly confident about
    # (CROSS with high logit). Read as a set of (q, t) keys and include in
    # the kept rows on top of the bits-based hard set.
    extra_hard = None
    if extra_hard_negatives_tsv:
        eh = pd.read_csv(extra_hard_negatives_tsv, sep="\t",
                         usecols=["query", "target"])
        extra_hard = set(zip(eh["query"], eh["target"]))
        print(f"  loaded {len(extra_hard):,} model-confidence hard negatives "
              f"from {extra_hard_negatives_tsv}", flush=True)

    # ----- Pass 2: load full features with hard-CROSS filter --------------
    t_p2 = time.time()
    print(f"Pass 2: reading {len(shards)} shards with filter …", flush=True)
    parts = []
    n_in = 0
    n_hard = 0
    n_tail = 0
    n_extra = 0
    for i, sh in enumerate(shards):
        d = pq.read_table(sh, columns=needed).to_pandas()
        n_in += len(d)
        is_cross = (d["rel"].to_numpy() == "CROSS")
        # Per-row threshold via query map; queries with too few CROSS get -inf
        row_thresh = d["query"].map(thresholds).fillna(float("-inf")).to_numpy()
        bits_arr = d["bits"].to_numpy()
        is_hard = is_cross & (bits_arr >= row_thresh)
        if tail_prob_map is not None:
            p_row = d["query"].map(tail_prob_map).fillna(0.0).to_numpy()
            tail_coin = rng.random(len(d)) < p_row
        else:
            tail_coin = rng.random(len(d)) < cross_random_tail_frac
        is_tail = is_cross & (~is_hard) & tail_coin
        # Model-confidence hard negatives
        is_extra = np.zeros(len(d), dtype=bool)
        if extra_hard is not None:
            q_arr = d["query"].to_numpy()
            t_arr = d["target"].to_numpy()
            is_extra = is_cross & np.array(
                [(q, t) in extra_hard for q, t in zip(q_arr, t_arr)])
        keep = (~is_cross) | is_hard | is_tail | is_extra
        n_hard += int(is_hard.sum())
        n_tail += int(is_tail.sum())
        n_extra += int(is_extra.sum())
        parts.append(d.loc[keep].reset_index(drop=True))
        if (i + 1) % 10 == 0 or (i + 1) == len(shards):
            print(f"  shard {i+1}/{len(shards)}: cum_in={n_in:,}  "
                  f"hard={n_hard:,}  tail={n_tail:,}  extra={n_extra:,}",
                  flush=True)
    df = pd.concat(parts, ignore_index=True)
    del parts
    print(f"loaded {len(df):,} rows of {n_in:,} ({100*len(df)/max(n_in,1):.2f}%) "
          f"in {time.time() - t_p2:.1f}s "
          f"(hard={n_hard:,}  tail={n_tail:,}  extra={n_extra:,})", flush=True)
    print("  relation tally:", dict(df["rel"].value_counts()), flush=True)

    # ----- Optional: drop excluded folds entirely (query OR target) -------
    # SCOP folds with known cross-fold homology (e.g. Rossmann-and-related)
    # inflate false positives; excluding them removes those domains from the
    # train/val/test/holdout universe as both queries and targets.
    if exclude_folds:
        if scop_lookup_path is None:
            raise ValueError("exclude_folds requires scop_lookup_path")
        excl = set(exclude_folds)
        sc = pd.read_csv(scop_lookup_path, sep="\t", header=None,
                         names=["domain", "scop"])
        sc["fold"] = sc["scop"].str.split(".").str[:2].str.join(".")
        d2fold = dict(zip(sc["domain"], sc["fold"]))
        qf = df["query"].map(d2fold); tf = df["target"].map(d2fold)
        drop = qf.isin(excl) | tf.isin(excl)
        print(f"  exclude_folds={sorted(excl)}: dropping {int(drop.sum()):,} of "
              f"{len(df):,} rows ({100*drop.mean():.2f}%)", flush=True)
        df = df[~drop].reset_index(drop=True)
        print("  relation tally after exclude:", dict(df["rel"].value_counts()), flush=True)

    # Hold-out FOLD: same domain set but no FAM/SFAM/CROSS overlap.
    fold_df  = df[df["rel"] == "FOLD"].reset_index(drop=True)
    train_df = df[df["rel"].isin(("FAM", "SFAM", "CROSS"))].reset_index(drop=True)
    train_df["label"] = (train_df["rel"].isin(("FAM", "SFAM"))).astype(np.float32)

    if explicit_split_json:
        with open(explicit_split_json) as fh:
            sp = json.load(fh)
        train_q = set(sp["train"]); val_q = set(sp["val"]); test_q = set(sp["test"])
        print(f"  explicit split ({explicit_split_json}): train_q={len(train_q):,}  "
              f"val_q={len(val_q):,}  test_q={len(test_q):,}", flush=True)
    elif split_by_fold:
        if scop_lookup_path is None:
            raise ValueError("split_by_fold=True requires scop_lookup_path")
        train_q, val_q, test_q = _split_by_fold(
            train_df["query"].to_numpy(), scop_lookup_path,
            val_pct, test_pct, seed=seed)
        print(f"  fold-out split: train_q={len(train_q):,}  "
              f"val_q={len(val_q):,}  test_q={len(test_q):,}",
              flush=True)
    else:
        train_q, val_q, test_q = _split_queries(
            train_df["query"].to_numpy(), val_pct, test_pct, seed=seed)
    train_mask = train_df["query"].isin(train_q)
    val_mask   = train_df["query"].isin(val_q)
    test_mask  = train_df["query"].isin(test_q)

    Xtr = train_df.loc[train_mask, feature_cols].to_numpy(dtype=np.float32)
    Ytr = train_df.loc[train_mask, "label"].to_numpy(dtype=np.float32)
    Xv  = train_df.loc[val_mask,   feature_cols].to_numpy(dtype=np.float32)
    Yv  = train_df.loc[val_mask,   "label"].to_numpy(dtype=np.float32)
    Xte = train_df.loc[test_mask,  feature_cols].to_numpy(dtype=np.float32)
    Yte = train_df.loc[test_mask,  "label"].to_numpy(dtype=np.float32)
    Xfd = fold_df[feature_cols].to_numpy(dtype=np.float32)
    Yfd = np.zeros(len(Xfd), dtype=np.float32)   # nominal "negative" but really intermediate

    # Z-score from the train split only.
    mean = Xtr.mean(axis=0)
    std  = Xtr.std(axis=0)
    std[std < 1e-6] = 1.0
    Xtr = (Xtr - mean) / std
    Xv  = (Xv  - mean) / std
    Xte = (Xte - mean) / std
    Xfd = (Xfd - mean) / std

    pos = float(Ytr.sum())
    neg = float(len(Ytr) - pos)
    pos_weight = neg / max(pos, 1.0)
    print(f"  train: pos={int(pos):,}  neg={int(neg):,}  pos_weight={pos_weight:.1f}",
          flush=True)

    if norm_stats_out:
        os.makedirs(os.path.dirname(norm_stats_out), exist_ok=True)
        with open(norm_stats_out, "w") as fh:
            json.dump({"feature_cols": feature_cols,
                       "mean": mean.tolist(), "std": std.tolist(),
                       "pos_weight": pos_weight}, fh)
        print(f"  saved norm stats → {norm_stats_out}", flush=True)

    to_td = lambda X, Y: TensorDataset(torch.from_numpy(X), torch.from_numpy(Y))
    return Splits(
        train=to_td(Xtr, Ytr), val=to_td(Xv, Yv), test=to_td(Xte, Yte),
        fold_holdout=to_td(Xfd, Yfd),
        feature_cols=feature_cols, feat_mean=mean, feat_std=std,
        pos_weight=pos_weight,
    )

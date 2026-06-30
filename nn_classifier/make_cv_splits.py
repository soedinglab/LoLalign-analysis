#!/usr/bin/env python3
"""Build 5-fold, fold-DISJOINT cross-validation splits.

The CV partitions the SCOP *folds* (first 2 sccs levels) into 5 balanced bins
(LPT: largest-fold-first, assigned to the currently-smallest bin by query count).
Each bin serves as the held-out TEST set once; within each iteration a VAL set
(~val_pct of all queries) is carved from whole folds of the remaining bins, and
everything else is TRAIN. No SCOP fold ever spans train/val/test, so there is no
fold-level leakage -- the same guarantee as phase1_data._split_by_fold, extended
to a complete 5-fold partition.

Query universe = unique non-FOLD queries in the (topo) parquet, identical to
make_split_q.py. Writes, for k in 0..4:
    outputs/cv5/cv{k}_split.json   {"train":[...],"val":[...],"test":[...]}
    outputs/cv5/cv{k}_test.txt     test query list (for the AUPR/ROC awks)
and outputs/cv5/cv_bins.json (bin -> folds, for the record).
"""
import glob, json, os
import numpy as np, pandas as pd, pyarrow.parquet as pq

SCOP = "/cbscratch/lolalignBenchmark/scope/data/scop_lookup.fix.tsv"
PARQ = "/cbscratch/lolalignBenchmark/nn_classifier/outputs/features_scalar_compact_topo_parquet"
OUT  = "/cbscratch/lolalignBenchmark/nn_classifier/outputs/cv5"
K = 5
VAL_PCT = 10
SEED = 43


def main():
    os.makedirs(OUT, exist_ok=True)
    shards = sorted(glob.glob(f"{PARQ}/shard_*.parquet"))
    qs = np.unique(np.concatenate([
        pq.read_table(s, columns=["query", "rel"]).to_pandas()
          .query("rel!='FOLD'")["query"].unique() for s in shards]))
    print(f"non-FOLD query universe: {len(qs):,}")

    sc = pd.read_csv(SCOP, sep="\t", header=None, names=["d", "s"])
    sc["fold"] = sc["s"].str.split(".").str[:2].str.join(".")
    d2f = dict(zip(sc.d, sc.fold))

    folds = {}
    for q in qs:
        f = d2f.get(q)
        if f is None:
            continue
        folds.setdefault(f, []).append(q)
    total = sum(len(v) for v in folds.values())
    print(f"folds: {len(folds):,}  queries with fold: {total:,}")

    # LPT balanced binning: largest fold first -> currently-smallest bin.
    fold_keys = sorted(folds, key=lambda f: (-len(folds[f]), f))
    bins = [[] for _ in range(K)]
    bin_n = [0] * K
    for f in fold_keys:
        j = int(np.argmin(bin_n))
        bins[j].append(f); bin_n[j] += len(folds[f])
    print("bin query counts:", bin_n, f"(target ~{total//K:,})")

    json.dump({f"bin{j}": sorted(bins[j]) for j in range(K)},
              open(f"{OUT}/cv_bins.json", "w"), indent=0)

    val_target = total * VAL_PCT / 100
    for k in range(K):
        test_folds = set(bins[k])
        rest = [f for j in range(K) if j != k for f in bins[j]]
        rng = np.random.default_rng(SEED + k)
        rng.shuffle(rest)
        val_folds, cum = set(), 0
        for f in rest:
            if cum >= val_target:
                break
            val_folds.add(f); cum += len(folds[f])
        train_folds = set(rest) - val_folds

        test  = sorted(q for f in test_folds  for q in folds[f])
        val   = sorted(q for f in val_folds   for q in folds[f])
        train = sorted(q for f in train_folds for q in folds[f])
        assert not (set(test) & set(val)) and not (set(test) & set(train)) \
            and not (set(val) & set(train)), "leak!"
        json.dump({"train": train, "val": val, "test": test},
                  open(f"{OUT}/cv{k}_split.json", "w"))
        with open(f"{OUT}/cv{k}_test.txt", "w") as fh:
            for q in test:
                fh.write(q + "\n")
        print(f"cv{k}: train={len(train):,} val={len(val):,} test={len(test):,} "
              f"(test_folds={len(test_folds)})")


if __name__ == "__main__":
    main()

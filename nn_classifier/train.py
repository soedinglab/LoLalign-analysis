#!/usr/bin/env python3
"""Train HomologyClassifier with BCE + pos_weight."""
from __future__ import annotations
import argparse
import os
import time

import numpy as np
import torch
import torch.nn as nn
import yaml
from torch.utils.data import DataLoader
from sklearn.metrics import roc_auc_score, average_precision_score

from data import load_splits
from model import HomologyClassifier


@torch.no_grad()
def eval_loader(model, loader, device):
    model.eval()
    logits_all, y_all = [], []
    for x, y in loader:
        logits_all.append(model(x.to(device)).cpu())
        y_all.append(y)
    model.train()
    logits = torch.cat(logits_all).numpy()
    y = torch.cat(y_all).numpy()
    if y.sum() == 0 or y.sum() == len(y):
        return {"loss": float("nan"), "roc_auc": float("nan"), "pr_auc": float("nan")}
    bce = nn.functional.binary_cross_entropy_with_logits(
        torch.from_numpy(logits), torch.from_numpy(y)
    ).item()
    return {
        "loss": bce,
        "roc_auc": roc_auc_score(y, logits),
        "pr_auc":  average_precision_score(y, logits),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="phase1_config.yaml")
    args = ap.parse_args()
    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    torch.manual_seed(cfg["train"]["seed"])
    np.random.seed(cfg["train"]["seed"])

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}", flush=True)

    from data import DEFAULT_FEATURE_COLS
    feat_cols = cfg["data"].get("feature_cols") or list(DEFAULT_FEATURE_COLS)
    splits = load_splits(
        parquet_dir=cfg["data"]["parquet_dir"],
        val_pct=cfg["data"]["val_pct"],
        test_pct=cfg["data"]["test_pct"],
        seed=cfg["data"]["seed"],
        norm_stats_out=cfg["data"]["norm_stats_out"],
        feature_cols=feat_cols,
        cross_top_k_per_query=cfg["data"].get("cross_top_k_per_query", 300),
        cross_random_tail_frac=cfg["data"].get("cross_random_tail_frac", 0.005),
        cross_tail_per_query=cfg["data"].get("cross_tail_per_query", 0),
        extra_hard_negatives_tsv=cfg["data"].get("extra_hard_negatives_tsv"),
        split_by_fold=cfg["data"].get("split_by_fold", False),
        scop_lookup_path=cfg["data"].get("scop_lookup_path"),
        exclude_folds=cfg["data"].get("exclude_folds"),
        explicit_split_json=cfg["data"].get("explicit_split_json"),
    )
    in_dim = splits.train.tensors[0].shape[1]
    print(f"in_dim = {in_dim}", flush=True)

    train_loader = DataLoader(splits.train, batch_size=cfg["train"]["batch_size"],
                              shuffle=True, num_workers=cfg["train"]["num_workers"],
                              pin_memory=True, drop_last=True)
    val_loader = DataLoader(splits.val, batch_size=cfg["train"]["batch_size"] * 4,
                            shuffle=False, num_workers=cfg["train"]["num_workers"],
                            pin_memory=True)

    model = HomologyClassifier(in_dim=in_dim,
                               hidden=tuple(cfg["model"]["hidden"]),
                               dropout=cfg["model"]["dropout"]).to(device)
    print(f"Parameters: {sum(p.numel() for p in model.parameters()):,}", flush=True)

    pos_weight = torch.tensor([splits.pos_weight], device=device)
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=pos_weight)

    optimizer = torch.optim.AdamW(model.parameters(),
                                  lr=cfg["train"]["lr"],
                                  weight_decay=cfg["train"]["weight_decay"])
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=cfg["train"]["num_epochs"])

    os.makedirs(cfg["output"]["log_dir"], exist_ok=True)
    log_path = os.path.join(cfg["output"]["log_dir"], "train.log")
    log_f = open(log_path, "w")
    log_f.write("epoch\ttrain_loss\tval_loss\tval_roc_auc\tval_pr_auc\n")

    best_val_auc = -1.0
    t0 = time.time()
    for epoch in range(1, cfg["train"]["num_epochs"] + 1):
        model.train()
        total, n_seen = 0.0, 0
        for x, y in train_loader:
            x = x.to(device, non_blocking=True)
            y = y.to(device, non_blocking=True)
            logits = model(x)
            loss = loss_fn(logits, y)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(),
                                           cfg["train"]["grad_clip"])
            optimizer.step()
            total += float(loss) * y.numel()
            n_seen += y.numel()
        scheduler.step()
        train_loss = total / max(n_seen, 1)

        m = eval_loader(model, val_loader, device)
        msg = (f"[ep {epoch:02d}/{cfg['train']['num_epochs']}] "
               f"train={train_loss:.4f}  val={m['loss']:.4f}  "
               f"ROC={m['roc_auc']:.4f}  PR={m['pr_auc']:.4f}  "
               f"elapsed={time.time()-t0:.0f}s")
        print(msg, flush=True)
        log_f.write(f"{epoch}\t{train_loss:.6f}\t{m['loss']:.6f}\t"
                    f"{m['roc_auc']:.6f}\t{m['pr_auc']:.6f}\n")
        log_f.flush()

        if m["roc_auc"] > best_val_auc:
            best_val_auc = m["roc_auc"]
            ckpt = {
                "model_state": model.state_dict(),
                "in_dim": in_dim,
                "config": cfg,
                "feature_cols": splits.feature_cols,
                "feat_mean": splits.feat_mean.tolist(),
                "feat_std":  splits.feat_std.tolist(),
                "pos_weight": splits.pos_weight,
                "best_val_roc_auc": best_val_auc,
                "epoch": epoch,
            }
            torch.save(ckpt, cfg["output"]["checkpoint"])
            print(f"  ↑ improved, saved → {cfg['output']['checkpoint']}",
                  flush=True)

    log_f.close()
    print(f"Done. Best val ROC AUC: {best_val_auc:.4f}. "
          f"Total {time.time()-t0:.0f}s.", flush=True)


if __name__ == "__main__":
    main()

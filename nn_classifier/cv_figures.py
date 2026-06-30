#!/usr/bin/env python3
"""Benchmark-style CV plots for the v20 homology classifier (5-fold,
fold-disjoint). Held-out test folds are pooled (each query scored by the model
that did NOT train on it). Four panels in the LoLalign-benchmark style:
  a  ROC1 sensitivity-to-first-FP curves (FAM/SFAM/FOLD)
  b  precision-recall curves (mean over folds)
  c  mean sensitivity vs AUPR scatter (per SCOP level)
  d  per-fold AUPR stability
Run on the frontend."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    'figure.figsize': (8, 6),
    'font.size': 12,
    'axes.labelsize': 20,
    'axes.titlesize': 16,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'legend.fontsize': 14,
    'lines.linewidth': 2.5,
    'figure.dpi': 100,
})

FDR = "/cbscratch/lolalignBenchmark/nn_classifier/outputs/cv5"
OUT = "/cbscratch/lolalignBenchmark/nn_classifier"
MODEL = "v20"
NFOLD = 5
# SCOP level -> (legend label, colour, marker)
LEVELS = [("FAM", "Family", "#2ca02c", "o"),
          ("SFAM", "Super Family", "#1f77b4", "^"),
          ("FOLD", "Fold", "#d62728", "d")]
LVL = [l[0] for l in LEVELS]


def test_queries(k):
    return set(l.strip() for l in open(f"{FDR}/cv{k}_test.txt") if l.strip())


def pooled_bench(model):
    """Union of held-out test-query rows across folds (one row per query)."""
    frames = []
    for k in range(NFOLD):
        tq = test_queries(k)
        d = pd.read_csv(f"{FDR}/bench_cv{k}_{model}/bench_model.tsv", sep="\t")
        frames.append(d[d.NAME.isin(tq)])
    return pd.concat(frames, ignore_index=True)


def fold_pr(model, k, lv):
    d = pd.read_csv(f"{FDR}/fdr_cv{k}_{model}.tsv", sep="\t")
    r = d[f"RECALL_{lv}"].to_numpy(float); p = d[f"PREC_{lv}"].to_numpy(float)
    i = np.argsort(r)
    r, p = r[i], p[i]
    r = np.concatenate([[0.0], r]); p = np.concatenate([[p[0]], p])
    return r, p


def aupr(model, k, lv):
    r, p = fold_pr(model, k, lv)
    return float(np.trapezoid(p, r))


GRID = np.linspace(0, 1, 201)


def mean_pr(model, lv):
    curves = [np.interp(GRID, *fold_pr(model, k, lv)) for k in range(NFOLD)]
    C = np.array(curves)
    return C.mean(0), C.std(0)


def panel_letter(ax, letter, x=-0.1):
    ax.text(x, 1.05, letter, transform=ax.transAxes, fontsize=24,
            fontweight='bold', va='top', ha='right')


def style(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# ---- data ----
bench = pooled_bench(MODEL)
# per-fold mean sensitivity (for scatter error bars) + pooled mean
sens_fold = {lv: [] for lv in LVL}
for k in range(NFOLD):
    tq = test_queries(k)
    d = pd.read_csv(f"{FDR}/bench_cv{k}_{MODEL}/bench_model.tsv", sep="\t")
    d = d[d.NAME.isin(tq)]
    for lv in LVL:
        sens_fold[lv].append(float(d[lv].mean()))
aupr_fold = {lv: [aupr(MODEL, k, lv) for k in range(NFOLD)] for lv in LVL}

# Comparison methods for panels a/b: v20 (the NN) vs raw LoL-align bits.
# color = SCOP level, linestyle = method (v20 solid, LoL-align dashed).
METHODS_AB = [("v20", "LoL-align Probability", "-"), ("lol", "LoL-align", "--")]
bench_by_model = {m: pooled_bench(m) for m, *_ in METHODS_AB}


def pooled_pr(model):
    d = pd.read_csv(f"{FDR}/fdr_pooled_{model}.tsv", sep="\t")
    return d.apply(pd.to_numeric, errors="coerce").dropna().reset_index(drop=True)


def ab_legend(ax, loc):
    lh = [Line2D([0], [0], color=c, lw=2.5, marker=mk, markersize=9, label=lab)
          for _, lab, c, mk in LEVELS]
    mh = [Line2D([0], [0], color='0.25', lw=2.5, ls=ls, label=mlab)
          for _, mlab, ls in METHODS_AB]
    ax.legend(handles=lh + mh, frameon=True, fancybox=True, shadow=True,
              loc=loc, fontsize=13)


# ---- a: ROC1 curves (v20 vs LoL-align, per SCOP level) ----
fig, ax = plt.subplots(); style(ax)
for lv, lab, col, mk in LEVELS:
    for mdl, mlab, ls in METHODS_AB:
        v = np.sort(bench_by_model[mdl][lv].to_numpy(float))[::-1]
        x = (np.arange(len(v)) + 1) / len(v)
        ax.plot(x, v, color=col, linestyle=ls, marker=mk, markersize=9,
                markevery=0.1, markeredgewidth=0.5, markeredgecolor='white')
ax.set_xlabel("Fraction of queries")
ax.set_ylabel("Sensitivity up to the first FP")
ax.set_ylim(0, 1.05)
ab_legend(ax, 'lower left')
panel_letter(ax, 'a')
plt.grid(False); plt.tight_layout()
for ext in ("png", "pdf"):
    plt.savefig(f"{OUT}/fig_cv_{MODEL}_roc1.{ext}", dpi=300, bbox_inches='tight')
plt.close()

# ---- b: PR curves (single pooled curve per method, all held-out test pairs) ----
# fdr_pooled_*.tsv = the s1000 fdr awk run once over the global score-sorted
# merge of held-out predictions (slurm/cv_pool_pr_{v20,lol}.sh).
pr_by_model = {m: pooled_pr(m) for m, *_ in METHODS_AB}
pooled_aupr = {m: {} for m, *_ in METHODS_AB}
fig, ax = plt.subplots(); style(ax)
for lv, lab, col, mk in LEVELS:
    for mdl, mlab, ls in METHODS_AB:
        d = pr_by_model[mdl]
        r = d[f"RECALL_{lv}"].to_numpy(float)
        p = d[f"PREC_{lv}"].to_numpy(float)
        i = np.argsort(r); r, p = r[i], p[i]
        r = np.concatenate([[0.0], r]); p = np.concatenate([[p[0]], p])
        pooled_aupr[mdl][lv] = float(np.trapezoid(p, r))
        ax.plot(r, p, color=col, linestyle=ls, marker=mk, markersize=9,
                markevery=0.1, markeredgewidth=0.5, markeredgecolor='white')
ax.set_xlabel("Recall")
ax.set_ylabel("Precision")
ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
panel_letter(ax, 'b')
plt.grid(False); plt.tight_layout()
for ext in ("png", "pdf"):
    plt.savefig(f"{OUT}/fig_cv_{MODEL}_pr.{ext}", dpi=300, bbox_inches='tight')
plt.close()

# ---- c: per-fold AUPR stability (v20 solid vs LoL-align dashed) ----
aupr_fold_by_model = {m: {lv: [aupr(m, k, lv) for k in range(NFOLD)] for lv in LVL}
                      for m, *_ in METHODS_AB}
fig, ax = plt.subplots(); style(ax)
folds = np.arange(1, NFOLD + 1)
for lv, lab, col, mk in LEVELS:
    for mdl, mlab, ls in METHODS_AB:
        ax.plot(folds, aupr_fold_by_model[mdl][lv], color=col, linestyle=ls,
                marker=mk, markersize=11, markeredgewidth=0.5, markeredgecolor='white')
ax.set_xlabel("CV fold")
ax.set_ylabel("Precision-Recall AUC")
ax.set_xticks(folds); ax.set_ylim(0, 1)
ab_legend(ax, 'upper right')
panel_letter(ax, 'c')
plt.grid(False); plt.tight_layout()
for ext in ("png", "pdf"):
    plt.savefig(f"{OUT}/fig_cv_{MODEL}_folds.{ext}", dpi=300, bbox_inches='tight')
plt.close()

# ---- d: mean sensitivity vs AUPR scatter (supplementary) ----
fig, ax = plt.subplots(); style(ax)
for lv, lab, col, mk in LEVELS:
    xm = np.mean(sens_fold[lv]); xs = np.std(sens_fold[lv])
    ym = np.mean(aupr_fold[lv]); ys = np.std(aupr_fold[lv])
    ax.errorbar(xm, ym, xerr=xs, yerr=ys, color=col, marker=mk, markersize=15,
                markeredgewidth=0.5, markeredgecolor='white', capsize=4,
                elinewidth=1.5, label=lab, linestyle='none')
ax.set_xlabel("Avg. sensitivity up to the first FP")
ax.set_ylabel("Precision-Recall AUC")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.legend(frameon=True, fancybox=True, shadow=True, loc='lower right', fontsize=15)
panel_letter(ax, 'd')
plt.grid(False); plt.tight_layout()
for ext in ("png", "pdf"):
    plt.savefig(f"{OUT}/fig_cv_{MODEL}_scatter.{ext}", dpi=300, bbox_inches='tight')
plt.close()

# ---- audit table ----
print(f"# 5-fold CV (pooled held-out test)  n_queries={len(bench)}")
print(f"{'level':6} {'method':10} {'sens(pooled)':14} {'AUPR(pooled)':12}")
for lv, lab, *_ in LEVELS:
    for mdl, mlab, _ in METHODS_AB:
        pooled_sens = float(bench_by_model[mdl][lv].to_numpy(float).mean())
        print(f"{lv:6} {mlab:10} {pooled_sens:.4f}        {pooled_aupr[mdl][lv]:.4f}")
print("v20 AUPR per-fold mean±sd: " +
      ", ".join(f"{lv} {np.mean(aupr_fold[lv]):.3f}±{np.std(aupr_fold[lv]):.3f}" for lv in LVL))
print("wrote fig_cv_v20_roc1 / _pr / _scatter / _folds (.png/.pdf)")

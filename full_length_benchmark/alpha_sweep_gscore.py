#!/usr/bin/env python3
"""Sweep alpha as a score parameter. Each method keeps its native per-query
ranking; the per-hit score becomes

    S(alpha) = aln_lddt^alpha * q_lddt^(1-alpha)

At alpha=0.5  -> S = sqrt(aln*q)   (the current G-score)
At alpha=0    -> S = q_lddt
At alpha=1    -> S = aln_lddt

Benchmark for each method+alpha: walk each query's hits in the method's
native order, accumulate (S - 0.5) / 0.5 until the first S < 0.5
(patience 1). AUC = trapezoid of the per-query-max-normalized curve.

Native order per method:
  Foldseek, LoL cs3/s3, Dali, TM-align -> intrinsic file row order
  LoL v20                              -> sorted by v20_prob desc
  SoftAlign                            -> pre-sorted by S(alpha) desc
                                          (no intrinsic ranking; resort per alpha)
"""
import os, sys, time
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Match the full-length benchmark figure style.
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

t0 = time.time()

# (name, tsv_path, native_sort_kind)
#   native_sort_kind = 'intrinsic' | 'v20' | 'softalign'
JOBS_FULL = [
    ('Foldseek',         '/cbscratch/rad/e_bench_archive/foldseek/fs_results/fs_calc_full', 'intrinsic'),
    ('LoL-align (rs1)',  '/cbscratch/lolalignBenchmark/lol_calc_full_cs1.tsv', 'intrinsic'),
    ('LoL-align (cs3)',  '/cbscratch/rad/e_bench_archive/lol_aln/results/lol_calc_full_md_f', 'intrinsic'),
    ('LoL-align (rs5)',  '/cbscratch/lolalignBenchmark/lol_calc_full_cs5.tsv', 'intrinsic'),
    ('LoL-align (rs10)', '/cbscratch/lolalignBenchmark/lol_calc_full_cs10.tsv', 'intrinsic'),
    ('LoL-align (v20)',  '/cbscratch/lolalignBenchmark/lol_calc_full_md_cs3_v20pred.tsv', 'v20'),
    ('Dali',             '/cbscratch/rad/e_bench_archive/dali/final_results/dali_calc_full', 'intrinsic'),
    ('TM-align',         '/cbscratch/rad/e_bench_archive/tm_aln/final_results/tm_aln_qtm_final', 'intrinsic'),
    ('SoftAlign',        '/cbscratch/lolalignBenchmark/softalign_score_calc_full.tsv', 'softalign'),
]
JOBS_LOW = [
    ('Foldseek',         '/cbscratch/lolalignBenchmark/fs_calc_full_low_plddt.tsv', 'intrinsic'),
    ('LoL-align (rs1)',  '/cbscratch/lolalignBenchmark/lol_calc_full_low_plddt_seed1.tsv', 'intrinsic'),
    ('LoL-align (s3)',   '/cbscratch/lolalignBenchmark/lol_calc_full_low_plddt.tsv', 'intrinsic'),
    ('LoL-align (rs5)',  '/cbscratch/lolalignBenchmark/lol_calc_full_low_plddt_seed5.tsv', 'intrinsic'),
    ('LoL-align (rs10)', '/cbscratch/lolalignBenchmark/lol_calc_full_low_plddt_seed10.tsv', 'intrinsic'),
    ('LoL-align (v20)',  '/cbscratch/lolalignBenchmark/lol_calc_full_low_plddt_md_cs3_v20pred.tsv', 'v20'),
    ('Dali',             '/cbscratch/lolalignBenchmark/dali_lol_score_calc_low_plddt.tsv', 'intrinsic'),
    ('TM-align',         '/cbscratch/lolalignBenchmark/tm_calc_full_low_plddt_qtm.tsv', 'intrinsic'),
    ('SoftAlign',        '/cbscratch/lolalignBenchmark/softalign_score_calc_low_plddt.tsv', 'softalign'),
]
COL = {'Foldseek': '#00AA00',
       'LoL-align (rs1)': '#FFB3B3', 'LoL-align (cs3)': '#FF0000', 'LoL-align (s3)': '#FF0000',
       'LoL-align (rs5)': '#B30000', 'LoL-align (rs10)': '#660000',
       'LoL-align (v20)': '#FF6B8A', 'Dali': '#FFA500', 'TM-align': '#138BD1',
       'SoftAlign': '#6A0DAD'}
MK  = {'Foldseek': 'd',
       'LoL-align (rs1)': '*', 'LoL-align (cs3)': '*', 'LoL-align (s3)': '*',
       'LoL-align (rs5)': '*', 'LoL-align (rs10)': '*',
       'LoL-align (v20)': 'v', 'Dali': 'o', 'TM-align': '^', 'SoftAlign': 'h'}
# legend display names: cs3/s3 -> rs3, rename v20 -> probability
LAB = {'LoL-align (cs3)': 'LoL-align (rs3)', 'LoL-align (s3)': 'LoL-align (rs3)',
       'LoL-align (v20)': 'LoL-align probability'}

ALPHAS = np.linspace(0.0, 1.0, 21)
THR = 0.5
THRESHOLDS = [0.4, 0.5, 0.6]
# Panel letters for the threshold-0.4/0.5/0.6 sweep figures.
LETTERS = {0.4: 'p', 0.5: 'q', 0.6: 'r'}


def load_softmax(path):
    """SoftAlign's native per-query ranking = its softmax score.
    Low-pLDDT scores live in a merged TSV; high-pLDDT in per-query CSVs."""
    import glob
    if 'low_plddt' in path:
        return pd.read_csv('/cbscratch/lolalignBenchmark/softalign_softmax_scores.tsv',
                           sep='\t', header=None, names=['query', 'target', 'softmax_score'])
    rows = []
    hp = '/cbscratch/lolalignBenchmark/e_bench/softalign/softmax_scores_highplddt/output'
    for f in glob.glob(hp + '/scores_sorted_*.csv'):
        q = os.path.basename(f)[len('scores_sorted_'):-len('.csv')]
        d = pd.read_csv(f); d.columns = ['target', 'softmax_score']; d['query'] = q
        rows.append(d[['query', 'target', 'softmax_score']])
    return pd.concat(rows, ignore_index=True)


def load(name, path, kind):
    df = pd.read_csv(path, sep='\t', low_memory=False)
    df = df[df['query'] != df['target']].copy()
    df['__a__'] = df['lddt_all_filtered'].astype(np.float64).clip(lower=1e-12)
    df['__q__'] = (df['lddt_all_filtered'] * df['aln_len_all_filtered'] / df['qlen_all_filtered']).astype(np.float64).clip(lower=1e-12)
    # All rankings are fixed once (no per-alpha resort).
    if kind == 'v20':
        df = df.sort_values(['query', 'v20_prob'], ascending=[True, False], kind='stable')
    elif kind == 'softalign':
        # native ranking = SoftAlign softmax score (replaces the old G-score stand-in)
        sm = load_softmax(path)
        df = df.merge(sm, on=['query', 'target'], how='left')
        nmiss = df['softmax_score'].isna().sum()
        if nmiss:
            print(f"    [softalign] {nmiss}/{len(df)} hits without a softmax score (ranked last)", flush=True)
        df = df.sort_values(['query', 'softmax_score'], ascending=[True, False],
                            kind='stable', na_position='last')
    else:  # intrinsic
        df = df.sort_values('query', kind='stable')
    return df.reset_index(drop=True)


def first_fp_cum(S, q_ix, threshold=THR):
    denom = 1.0 - threshold
    n = len(S)
    out = {}
    i = 0
    while i < n:
        q = q_ix[i]; j = i; total = 0.0; broke = False
        while j < n and q_ix[j] == q:
            s = S[j]
            if not broke and s < threshold:
                broke = True
            if not broke:
                total += (s - threshold) / denom
            j += 1
        out[q] = total
        i = j
    return out


def auc_for_methods(per_q_dict_by_method):
    R = pd.DataFrame(per_q_dict_by_method).fillna(0)
    mx = R.max(axis=1) + 1e-10
    rel = R.div(mx, axis=0)
    out = {}
    for m in rel.columns:
        v = np.sort(rel[m].values)[::-1]
        x = (np.arange(len(v)) + 1) / len(v)
        out[m] = float(np.trapezoid(v, x))
    return out


def sweep_auc(loaded, thr):
    AUC = {}
    for ai, a in enumerate(ALPHAS):
        per_method = {}
        for name, df, kind in loaded:
            ad = df['__a__'].to_numpy(); qd = df['__q__'].to_numpy()
            S = np.power(ad, a) * np.power(qd, 1.0 - a)
            q_ix = df['query'].to_numpy()
            per_method[name] = first_fp_cum(S, q_ix, threshold=thr)
        AUC[a] = auc_for_methods(per_method)
        if ai % 5 == 0 or ai == len(ALPHAS) - 1:
            print(f"    thr={thr:.1f} alpha={a:.2f}  ({time.time()-t0:.0f}s)", flush=True)
    return AUC


def run_half(label, jobs, out_base, thresholds=THRESHOLDS, letter='e'):
    print(f"\n=== {label} ===", flush=True)
    loaded = [(name, load(name, p, kind), kind) for name, p, kind in jobs]
    for name, df, _ in loaded:
        print(f"  {name:18s} rows={len(df):,} ({time.time()-t0:.0f}s)", flush=True)

    for thr in thresholds:
        AUC = sweep_auc(loaded, thr)
        suffix = f"_{int(round(thr * 100)):02d}"
        plot_sweep(AUC, jobs, out_base + suffix + ".png",
                   out_base + suffix + ".tsv", letter=LETTERS.get(thr, letter))


def plot_sweep(AUC, jobs, out_png, out_tsv, letter='e'):
    method_names = [name for name, *_ in jobs]
    out = pd.DataFrame(
        [[AUC[a][m] for m in method_names] for a in ALPHAS],
        columns=method_names, index=[f"{a:.2f}" for a in ALPHAS])
    out.index.name = 'alpha'
    out.to_csv(out_tsv, sep='\t')
    print(f"  wrote {out_tsv}", flush=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for m in method_names:
        is_prob = m.endswith('(v20)')
        is_lol = m.startswith('LoL-align')
        ax.plot(ALPHAS, [AUC[a][m] for a in ALPHAS],
                marker=MK[m], color=COL[m],
                linewidth=3.0 if is_prob else 2.5,
                markersize=15 if is_lol else 10,
                markeredgewidth=0.5, markeredgecolor='white',
                linestyle='--' if is_prob else '-',
                alpha=0.8 if is_prob else 1.0,
                zorder=1 if is_prob else 3, label=LAB.get(m, m))
    ax.set_ylabel('Area under the Curve')
    ax.set_xlabel(' ')
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(['sensitivity', 'G-score', 'precision'], fontsize=18)
    # Vertical dotted guides at the three named x positions.
    for xv in [0, 0.5, 1.0]:
        ax.axvline(xv, color='#999999', linestyle=':', linewidth=2.2, zorder=0)
    for lbl, ha in zip(ax.get_xticklabels(), ['left', 'center', 'right']):
        lbl.set_ha(ha)
    ax.text(-0.1, 1.05, letter, transform=ax.transAxes, fontsize=24, fontweight='bold', va='top', ha='right')
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.savefig(out_png.replace('.png', '.pdf'), dpi=300, bbox_inches='tight', format='pdf')
    plt.close()
    print(f"  wrote {out_png}", flush=True)


which = sys.argv[1] if len(sys.argv) > 1 else 'both'
if which in ('full', 'both'):
    run_half('High pLDDT', JOBS_FULL,
             '/cbscratch/lolalignBenchmark/alpha_sweep_full',
             thresholds=[0.4, 0.5, 0.6])
if which in ('low', 'both'):
    run_half('Low pLDDT', JOBS_LOW,
             '/cbscratch/lolalignBenchmark/alpha_sweep_low',
             thresholds=[0.4, 0.5, 0.6])
print(f"DONE ({time.time()-t0:.0f}s)")

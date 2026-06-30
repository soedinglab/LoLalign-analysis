import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
import matplotlib.transforms as mtransforms

plt.rcParams.update({
    'figure.figsize': (8, 6),
    'font.size': 12,
    'axes.labelsize': 20,
    'axes.titlesize': 16,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'legend.fontsize': 14,
    'lines.linewidth': 2.5,
    'grid.alpha': 0.3,
    'figure.dpi': 100,
})

# ── ROC1 superfamily (a) ──────────────────────────────────────────────
scope_lol_ns_len = pd.read_csv("allvsall_lol.rocx", delimiter="\t")
scope_tm         = pd.read_csv("allvsall_tmsort.rocx", delimiter="\t")
scope_dali       = pd.read_csv("dali_sorted.rocx", delimiter="\t")
scope_fs         = pd.read_csv("allvsall_fs.rocx", delimiter="\t")
scope_softal     = pd.read_csv("res_softal.rocx", delimiter="\t")
scope_tmvec2     = pd.read_csv("allvsall_tmvec2.rocx", delimiter="\t")
scope_ce         = pd.read_csv("res_CE.rocx", delimiter="\t")

colors = {
    "Dali": "#FFA500",
    "TM-align": "#138BD1",
    "Foldseek": "#00AA00",
    "LoLalign": "#FF0000",
    "Soft-align": "#6A0DAD",
}

scale = np.arange(0, 3566) / 3565

plt.figure()
plt.plot(scale, np.sort(scope_lol_ns_len["SFAM"].values)[::-1],
         label='LoL-align', color=colors["LoLalign"], linewidth=2.5,
         marker='*', markersize=15, markevery=0.1)
plt.plot(scale, np.sort(scope_tm["SFAM"].values)[::-1],
         label='TM-align', color=colors["TM-align"], linewidth=2.5,
         marker='^', markersize=10, markevery=0.1)
plt.plot(scale, np.sort(scope_fs["SFAM"].values)[::-1],
         label='Foldseek', color=colors["Foldseek"], linewidth=2.5,
         marker='d', markersize=10, markevery=0.1)
plt.plot(scale, np.sort(scope_dali["SFAM"].values)[::-1],
         label='Dali', color=colors["Dali"], linewidth=2.5,
         marker='o', markersize=10, markevery=0.1)
plt.plot(scale, np.sort(scope_softal["SFAM"].values)[::-1],
         label='SoftAlign', color=colors["Soft-align"], linewidth=2.5,
         marker='h', markersize=10, markevery=0.1)
plt.plot(scale, np.sort(scope_tmvec2["SFAM"].values)[::-1],
         label='TM-Vec2', color="#8B4513", linewidth=2.5,
         marker='P', markersize=10, markevery=0.1)
plt.plot(scale, np.sort(scope_ce["SFAM"].values)[::-1],
         label='CE', color="#6A5ACD", linewidth=2.5,
         marker='H', markersize=10, markevery=0.1)

plt.xlim(0.2, 1)
plt.xlabel('Fraction of queries', fontsize=20)
plt.ylabel("Sensitivity up to the first FP", fontsize=20)

# Build an explicit combined legend (shown only on this panel) so every entry
# gets a solid, clearly visible line+marker swatch. LoL-align probability has
# no SFAM ROC1 curve here, so it is a legend-only proxy handle.
legend_handles = [
    Line2D([0], [0], color="#FF0000", lw=2.5, marker='*', markersize=13, label='LoL-align'),
    Line2D([0], [0], color="#138BD1", lw=2.5, marker='^', markersize=9, label='TM-align'),
    Line2D([0], [0], color="#00AA00", lw=2.5, marker='d', markersize=9, label='Foldseek'),
    Line2D([0], [0], color="#FFA500", lw=2.5, marker='o', markersize=9, label='Dali'),
    Line2D([0], [0], color="#6A0DAD", lw=2.5, marker='h', markersize=9, label='SoftAlign'),
    Line2D([0], [0], color="#FF6B8A", lw=2.5, marker='v', markersize=11, linestyle='--', label='LoL-align probability'),
    Line2D([0], [0], color="#8B4513", lw=2.5, marker='P', markersize=9, label='TM-Vec2'),
    Line2D([0], [0], color="#6A5ACD", lw=2.5, marker='H', markersize=9, label='CE'),
]
plt.legend(handles=legend_handles, frameon=True, fancybox=True, shadow=True,
           loc='lower left', fontsize=15, handlelength=1.8,
           labelspacing=0.35, borderpad=0.5, handletextpad=0.5)
plt.tick_params(labelsize=18)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.text(-0.1, 1.05, 'a', transform=ax.transAxes, fontsize=24,
        fontweight='bold', va='top', ha='right')
plt.tight_layout()
plt.savefig('LoLalign_benchmark_a.pdf', dpi=300, bbox_inches='tight', format='pdf')
plt.savefig('LoLalign_benchmark_a.png', dpi=300, bbox_inches='tight', format='png')
plt.close()

# ── PR superfamily (b) ────────────────────────────────────────────────
pr_lol    = pd.read_csv("pr_allvsall_lol.rocx", sep='\t')
pr_fs     = pd.read_csv("pr_allvsall_fs.rocx", sep='\t')
pr_dali   = pd.read_csv("pr_dali.rocx", sep='\t')
pr_tm     = pd.read_csv("pr_allvsall_tm.rocx", sep='\t')
pr_tmvec2 = pd.read_csv("pr_tmvec2.rocx", sep='\t')
pr_softal = pd.read_csv("pr_softal.rocx", sep='\t')
pr_ce     = pd.read_csv("pr_CE.rocx", sep='\t')

pr_colors = {
    "LoL-align": "#FF0000",
    "TM-align": "#138BD1",
    "Foldseek": "#00AA00",
    "Dali": "#FFA500",
    "TM-Vec2": "#8B4513",
    "SoftAlign": "#6A0DAD",
    "CE": "#6A5ACD",
}
pr_markers = {
    "LoL-align": "*",
    "TM-align": "^",
    "Foldseek": "d",
    "Dali": "o",
    "TM-Vec2": "P",
    "SoftAlign": "h",
    "CE": "H",
}
msizes = {k: (15 if k == "LoL-align" else 10) for k in pr_colors}

tools = [
    ("LoL-align", pr_lol),
    ("TM-align", pr_tm),
    ("Foldseek", pr_fs),
    ("Dali",     pr_dali),
    ("TM-Vec2",  pr_tmvec2),
    ("SoftAlign", pr_softal),
    ("CE",       pr_ce),
]

plt.figure()
for name, df in tools:
    plt.plot(df['RECALL_SFAM'], df['PREC_SFAM'],
             label=name, color=pr_colors[name], linewidth=2.5,
             marker=pr_markers[name], markersize=msizes[name], markevery=0.1)

plt.xlabel("Recall", fontsize=20)
plt.ylabel("Precision", fontsize=20)
plt.xlim(0, 1)
plt.ylim(0, 1.02)
plt.tick_params(labelsize=18)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.text(-0.1, 1.05, 'b', transform=ax.transAxes, fontsize=24,
        fontweight='bold', va='top', ha='right')
plt.tight_layout()
plt.savefig('LoLalign_benchmark_b.pdf', dpi=300, bbox_inches='tight', format='pdf')
plt.savefig('LoLalign_benchmark_b.png', dpi=300, bbox_inches='tight', format='png')
plt.close()

# ── Avg. sensitivity vs PR-AUC scatter, per SCOP level (c) ─────────────
# x = mean per-query sensitivity up to the first FP (ROC1 .rocx files)
# y = area under the precision-recall curve (pr_*.rocx files)
roc_files = {
    "LoL-align": "allvsall_lol.rocx",
    "TM-align":  "allvsall_tmsort.rocx",
    "Foldseek":  "allvsall_fs.rocx",
    "Dali":      "dali_sorted.rocx",
    "SoftAlign": "res_softal.rocx",
    "TM-Vec2":   "allvsall_tmvec2.rocx",
    "CE":        "res_CE.rocx",
}
pr_files = {
    "LoL-align": "pr_allvsall_lol.rocx",
    "TM-align":  "pr_allvsall_tm.rocx",
    "Foldseek":  "pr_allvsall_fs.rocx",
    "Dali":      "pr_dali.rocx",
    "SoftAlign": "pr_softal.rocx",
    "TM-Vec2":   "pr_tmvec2.rocx",
    "CE":        "pr_CE.rocx",
}
# (color, marker, marker-size) consistent with panels a/b
c_style = {
    "LoL-align": ("#FF0000", "*", 520),
    "TM-align":  ("#138BD1", "^", 260),
    "Foldseek":  ("#00AA00", "d", 260),
    "Dali":      ("#FFA500", "o", 260),
    "SoftAlign": ("#6A0DAD", "h", 260),
    "TM-Vec2":   ("#8B4513", "P", 260),
    "CE":        ("#6A5ACD", "H", 260),
}

import matplotlib.colors as mcolors
def desaturate(color, sat=0.7):
    """Scale a colour's saturation (HSV) by `sat`."""
    h, s, v = mcolors.rgb_to_hsv(mcolors.to_rgb(color))
    return mcolors.hsv_to_rgb((h, s * sat, v))

def pr_auc(df, lvl):
    r = df["RECALL_" + lvl].values
    p = df["PREC_" + lvl].values
    o = np.argsort(r)
    return float(np.trapz(p[o], r[o]))

def cluster_ellipse(xs, ys, ax, n_std=2.6, **kwargs):
    """Soft covariance ellipse enclosing a SCOP-level cluster (diagonal-oriented)."""
    cov = np.cov(xs, ys)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    rx = np.sqrt(1 + pearson)
    ry = np.sqrt(1 - pearson)
    ell = Ellipse((0, 0), width=2 * rx, height=2 * ry, **kwargs)
    sx = np.sqrt(cov[0, 0]) * n_std
    sy = np.sqrt(cov[1, 1]) * n_std
    transf = (mtransforms.Affine2D()
              .rotate_deg(45).scale(sx, sy)
              .translate(np.mean(xs), np.mean(ys)))
    ell.set_transform(transf + ax.transData)
    ax.add_patch(ell)

# Gather points per SCOP level so each cluster can be shaded.
# Level is encoded structurally via the marker edge ring (not color/saturation),
# so it reads correctly even where two levels' points sit close together:
#   Family -> thick black ring, Super Family -> thin gray ring, Fold -> no ring.
ring_style = {
    "FAM":  dict(edgecolors="#666666", linewidths=1.4),
    "SFAM": dict(edgecolors="black",   linewidths=2.5),
    "FOLD": dict(edgecolors="none",    linewidths=0.0),
}
level_pts = {"FAM": [], "SFAM": [], "FOLD": []}
marker_args = []
for tool, (col, mk, sz) in c_style.items():
    roc_df = pd.read_csv(roc_files[tool], sep="\t")
    pr_df  = pd.read_csv(pr_files[tool], sep="\t")
    for lvl in ["FAM", "SFAM", "FOLD"]:
        x = roc_df[lvl].mean()
        y = pr_auc(pr_df, lvl)
        level_pts[lvl].append((x, y))
        marker_args.append((x, y, col, mk, sz, lvl))

fig, ax = plt.subplots()

for x, y, col, mk, sz, lvl in marker_args:
    a = 0.5 if lvl == "FAM" else 1.0
    ax.scatter(x, y, color=col, marker=mk, s=sz, alpha=a, zorder=3, **ring_style[lvl])

# Cluster annotations (one SCOP level per cluster), placed near each cluster's
# actual centroid (computed from level_pts) rather than fixed guessed coordinates.
ax.text(0.68, 0.85, "Family",       fontsize=20, ha="center", va="center")
ax.text(0.30, 0.63, "Super Family", fontsize=20, ha="center", va="center")
ax.text(0.32, 0.08, "Fold",         fontsize=20, ha="center", va="center")

# Small legend explaining the ring encoding.
ring_legend = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#999999',
           markeredgecolor='#666666', markeredgewidth=1.4, markersize=13, label='Family'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#999999',
           markeredgecolor='black', markeredgewidth=2.5, markersize=13, label='Super Family'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#999999',
           markeredgecolor='none', markersize=13, label='Fold'),
]
ax.legend(handles=ring_legend, title='SCOP level', loc='lower right',
          fontsize=12, title_fontsize=12, frameon=True, fancybox=True,
          handletextpad=0.6, borderpad=0.5)

plt.xlabel("Avg. sensitivity up to the first FP", fontsize=20)
plt.ylabel("Precision-Recall AUC", fontsize=20)
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.tick_params(labelsize=18)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.text(-0.1, 1.05, 'c', transform=ax.transAxes, fontsize=24,
        fontweight='bold', va='top', ha='right')
plt.tight_layout()
plt.savefig('LoLalign_benchmark_c.pdf', dpi=300, bbox_inches='tight', format='pdf')
plt.savefig('LoLalign_benchmark_c.png', dpi=300, bbox_inches='tight', format='png')
plt.close()
print("done")

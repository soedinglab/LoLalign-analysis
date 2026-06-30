import os.path as osp
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

foldseek_3di  = pd.read_csv('/cbscratch/lolalignBenchmark/fs_calc_full_low_plddt.tsv', sep='\t')
lol_rs1       = pd.read_csv('/cbscratch/lolalignBenchmark/lol_calc_full_low_plddt_seed1.tsv', sep='\t')
foldseek_lol  = pd.read_csv('/cbscratch/lolalignBenchmark/lol_calc_full_low_plddt.tsv', sep='\t')
lol_rs5       = pd.read_csv('/cbscratch/lolalignBenchmark/lol_calc_full_low_plddt_seed5.tsv', sep='\t')
lol_rs10      = pd.read_csv('/cbscratch/lolalignBenchmark/lol_calc_full_low_plddt_seed10.tsv', sep='\t')
foldseek_lol_v20 = pd.read_csv('/cbscratch/lolalignBenchmark/lol_calc_full_low_plddt_md_cs3_v20pred.tsv', sep='\t', low_memory=False)
dali_aln      = pd.read_csv('/cbscratch/lolalignBenchmark/dali_lol_score_calc_low_plddt.tsv', sep='\t')
tm_aln_qtm    = pd.read_csv('/cbscratch/lolalignBenchmark/tm_calc_full_low_plddt_qtm.tsv', sep='\t')
softalign_aln = pd.read_csv('/cbscratch/lolalignBenchmark/softalign_score_calc_low_plddt.tsv', sep='\t')

print(f'Foldseek: {len(foldseek_3di)} alignments')
print(f'LoL-align (rs1): {len(lol_rs1)} alignments')
print(f'LoL-align (rs3): {len(foldseek_lol)} alignments')
print(f'LoL-align (rs5): {len(lol_rs5)} alignments')
print(f'LoL-align (rs10): {len(lol_rs10)} alignments')
print(f'LoL-align (v20): {len(foldseek_lol_v20)} alignments')
print(f'DALI: {len(dali_aln)} alignments')
print(f'TM-align: {len(tm_aln_qtm)} alignments')
print(f'SoftAlign: {len(softalign_aln)} alignments')

# Compute derived scores
alignments = [foldseek_3di, lol_rs1, foldseek_lol, lol_rs5, lol_rs10, dali_aln, tm_aln_qtm, softalign_aln]

for df in alignments:
    df['query_lddt_filter'] = (df['lddt_all_filtered'] * df['aln_len_all_filtered'] / df['qlen_all_filtered'])
    df['G_score_lddt_filter'] = np.sqrt(df['lddt_all_filtered'] * df['query_lddt_filter'])

# SoftAlign native ranking = softmax score (merged TSV for low-pLDDT).
_softmax = pd.read_csv('/cbscratch/lolalignBenchmark/softalign_softmax_scores.tsv',
                       sep='\t', header=None, names=['query', 'target', 'softmax_score'])
softalign_aln = softalign_aln.merge(_softmax, on=['query', 'target'], how='left')
nmiss = softalign_aln['softmax_score'].isna().sum()
print(f'SoftAlign: merged softmax scores, {nmiss} of {len(softalign_aln)} rows without a score (ranked last)')
softalign_aln = softalign_aln.sort_values(['query', 'softmax_score'], ascending=[True, False],
                                           kind='stable', na_position='last').reset_index(drop=True)
alignments[7] = softalign_aln

# LoL-align (v20): rs3 alignments re-ranked by predicted homology probability.
foldseek_lol_v20['query_lddt_filter'] = (foldseek_lol_v20['lddt_all_filtered'] * foldseek_lol_v20['aln_len_all_filtered'] / foldseek_lol_v20['qlen_all_filtered'])
foldseek_lol_v20['G_score_lddt_filter'] = np.sqrt(foldseek_lol_v20['lddt_all_filtered'] * foldseek_lol_v20['query_lddt_filter'])
foldseek_lol_v20 = foldseek_lol_v20.sort_values(['query', 'v20_prob'], ascending=[True, False], kind='stable').reset_index(drop=True)
alignments.append(foldseek_lol_v20)


def calculate_query_scores(df_optimized, mode='G_score_lddt_filter', patience=1, threshold_value=0.5):
    score_column = mode
    threshold = threshold_value
    denominator = 1 - threshold_value

    df_work = df_optimized[['query', score_column]].copy()

    def compute_group_score(group):
        scores = group.values
        total = 0
        below_threshold_count = 0
        for score in scores:
            if score < threshold:
                below_threshold_count += 1
                if below_threshold_count >= patience:
                    break
            else:
                total += (score - threshold) / denominator
        return round(float(total), 6)

    query_scores = df_work.groupby('query', sort=False)[score_column].apply(compute_group_score)
    return query_scores.to_dict()


def calculate_query_sensitivity(df_optimized, mode='G_score_lddt_filter'):
    score_column = mode
    threshold = 0.5

    df_work = df_optimized[['query', score_column]].copy()

    def compute_group_score(group):
        scores = group.values
        valid_mask = scores >= threshold
        if not valid_mask.any():
            return 0
        first_invalid = np.argmax(~valid_mask) if not valid_mask.all() else len(scores)
        if first_invalid > 0:
            return len(scores[:first_invalid])
        return 0

    query_scores = df_work.groupby('query', sort=False)[score_column].apply(compute_group_score)
    return query_scores.to_dict()


def calculate_query_aln_quality(df_optimized, mode='G_score_lddt_filter', threshold_value=0.5):
    score_column = mode
    threshold = threshold_value
    denominator = 1 - threshold_value

    df_work = df_optimized[['query', score_column]].copy()

    def compute_group_score(group):
        sorted_scores = group.sort_values(ascending=False).values
        total = 0
        for score in sorted_scores:
            if score >= threshold:
                total += (score - threshold) / denominator
            else:
                break
        return round(float(total), 6)

    query_scores = df_work.groupby('query', sort=False)[score_column].apply(compute_group_score)
    return query_scores.to_dict()


alignment_names = ['foldseek_3di', 'lol_rs1', 'foldseek_lol', 'lol_rs5', 'lol_rs10',
                   'dali_aln', 'tm_aln_qtm', 'softalign_aln', 'foldseek_lol_v20']
display_names = {
    'foldseek_3di':    'Foldseek',
    'lol_rs1':         'LoL-align (rs1)',
    'foldseek_lol':    'LoL-align (rs3)',
    'lol_rs5':         'LoL-align (rs5)',
    'lol_rs10':        'LoL-align (rs10)',
    'dali_aln':        'Dali',
    'tm_aln_qtm':      'TM-align',
    'softalign_aln':   'SoftAlign',
    'foldseek_lol_v20': 'LoL-align probability',
}

def compute_all_scores(alignments, alignment_names, mode, thresholds=[0.4, 0.5, 0.6]):
    results = {}
    for thr in thresholds:
        scores_dict = {}
        for df, name in zip(alignments, alignment_names):
            df_opt = df[['query', 'target', 'lddt_all_filtered', 'query_lddt_filter', 'G_score_lddt_filter']].copy()
            df_opt.drop(df_opt[df_opt['query'] == df_opt['target']].index, inplace=True)
            col_name = display_names[name]
            scores_dict[col_name] = calculate_query_scores(df_opt, mode=mode, patience=1, threshold_value=thr)
        result_df = pd.DataFrame(scores_dict).reset_index().rename(columns={'index': 'query'})
        result_df.fillna(0, inplace=True)
        results[thr] = result_df
    return results

# G-score LDDT
g_score_results = compute_all_scores(alignments, alignment_names, 'G_score_lddt_filter')
g_score_lddt_04 = g_score_results[0.4]
g_score_lddt_05 = g_score_results[0.5]
g_score_lddt_06 = g_score_results[0.6]

# Query LDDT
query_lddt_results = compute_all_scores(alignments, alignment_names, 'query_lddt_filter')
query_lddt_04 = query_lddt_results[0.4]
query_lddt_05 = query_lddt_results[0.5]
query_lddt_06 = query_lddt_results[0.6]

# Alignment LDDT
aln_lddt_results = compute_all_scores(alignments, alignment_names, 'lddt_all_filtered')
aln_lddt_04 = aln_lddt_results[0.4]
aln_lddt_05 = aln_lddt_results[0.5]
aln_lddt_06 = aln_lddt_results[0.6]

colors = {
    'Dali':                  '#FFA500',
    'TM-align':              '#138BD1',
    'Foldseek':              '#00AA00',
    'LoL-align (rs1)':       '#FFB3B3',
    'LoL-align (rs3)':       '#FF0000',
    'LoL-align (rs5)':       '#B30000',
    'LoL-align (rs10)':      '#660000',
    'SoftAlign':             '#6A0DAD',
    'LoL-align probability': '#FF6B8A',
}

markers = {
    'Dali': 'o', 'TM-align': '^', 'Foldseek': 'd',
    'LoL-align (rs1)': '*', 'LoL-align (rs3)': '*',
    'LoL-align (rs5)': '*', 'LoL-align (rs10)': '*',
    'SoftAlign': 'h', 'LoL-align probability': 'v',
}

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
    'figure.dpi': 100
})

SAVE_DIR = '/cbscratch/lolalignBenchmark/plots_low_plddt'
os.makedirs(SAVE_DIR, exist_ok=True)


def _plot_order(cols):
    """Draw the probability curve last so it sits on top of the overlapping LoL cluster."""
    return [c for c in cols if c != 'LoL-align probability'] + [c for c in cols if c == 'LoL-align probability']


def plot_relative_scores(df, title, save_name, ylabel='Relative G-score up to the first FP', label_letter='a'):
    """ROC-style relative score plot."""
    df_plot = df.copy()
    num_cols = df_plot.select_dtypes(include='number').columns
    df_plot['max_score'] = df_plot[num_cols].max(axis=1) + 1e-10
    sens_df = df_plot[num_cols].div(df_plot['max_score'], axis=0)
    
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    for col in _plot_order(sens_df.columns):
        sorted_vals = sens_df[col].sort_values(ascending=False).reset_index(drop=True)
        n_queries = len(sorted_vals)
        fraction_queries = (sorted_vals.index + 1) / n_queries
        is_prob = (col == 'LoL-align probability')
        ms = 15 if col.startswith('LoL-align') else 10
        ax.plot(fraction_queries, sorted_vals.values, label=col,
                color=colors.get(col, '#000000'),
                marker=markers.get(col),
                markevery=max(1, len(sorted_vals) // 12),
                markersize=ms, markeredgewidth=0.5, markeredgecolor='white',
                linewidth=3.0 if is_prob else 2.5,
                linestyle='--' if is_prob else '-',
                alpha=0.8 if is_prob else 1.0,
                zorder=20 if is_prob else 3)
    
    ax.set_xlabel('Fraction of queries')
    ax.set_ylabel(ylabel, y=0.41)
    ax.set_ylim(0, 1.05)
    ax.text(-0.1, 1.05, label_letter, transform=ax.transAxes, fontsize=24, fontweight='bold', va='top', ha='right')
    plt.legend(frameon=True, fancybox=True, shadow=True, loc='lower left', fontsize=13)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, f'{save_name}.pdf'), dpi=300, bbox_inches='tight', format='pdf')
    plt.savefig(os.path.join(SAVE_DIR, f'{save_name}.png'), dpi=300, bbox_inches='tight', format='png')
    plt.show()


def plot_cumulative_scores(df, title, save_name, ylabel='Cumulative G-score up to the first FP', label_letter='b'):
    """Cumulative score plot (log scale)."""
    df_plot = df.drop('query', axis=1).copy()
    df_plot = df_plot + 1
    
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    for col in _plot_order(df_plot.columns):
        sorted_vals = df_plot[col].sort_values(ascending=False).reset_index(drop=True)
        n_queries = len(sorted_vals)
        fraction_queries = (sorted_vals.index + 1) / n_queries
        is_prob = (col == 'LoL-align probability')
        ms = 15 if col.startswith('LoL-align') else 10
        ax.plot(fraction_queries, sorted_vals.values, label=col,
                color=colors.get(col, '#000000'),
                marker=markers.get(col),
                markevery=max(1, len(sorted_vals) // 13),
                markersize=ms, markeredgewidth=0.1, markeredgecolor='white',
                linewidth=3.0 if is_prob else 2.5,
                linestyle='--' if is_prob else '-',
                alpha=0.8 if is_prob else 1.0,
                zorder=20 if is_prob else 3)
    
    ax.set_xlabel('Fraction of queries')
    ax.set_ylabel(ylabel, y=0.41)
    ax.set_yscale('log')
    ax.set_ylim(0, max(df_plot.max()) * 2)
    ax.text(-0.11, 1.05, label_letter, transform=ax.transAxes, fontsize=24, fontweight='bold', va='top', ha='right')
    plt.legend(frameon=True, fancybox=True, shadow=True, loc='upper right', fontsize=13)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, f'{save_name}.pdf'), dpi=300, bbox_inches='tight', format='pdf')
    plt.savefig(os.path.join(SAVE_DIR, f'{save_name}.png'), dpi=300, bbox_inches='tight', format='png')
    plt.show()


plot_relative_scores(g_score_lddt_04, 'G-score LDDT (threshold 0.4)', 'ROC_G_score_04', label_letter='a')
plot_relative_scores(g_score_lddt_05, 'G-score LDDT (threshold 0.5)', 'ROC_G_score_05', label_letter='b')
plot_relative_scores(g_score_lddt_06, 'G-score LDDT (threshold 0.6)', 'ROC_G_score_06', label_letter='c')

plot_cumulative_scores(g_score_lddt_04, 'Cumulative G-score (threshold 0.4)', 'cumulative_G_score_04', label_letter='d')
plot_cumulative_scores(g_score_lddt_05, 'Cumulative G-score (threshold 0.5)', 'cumulative_G_score_05', label_letter='e')
plot_cumulative_scores(g_score_lddt_06, 'Cumulative G-score (threshold 0.6)', 'cumulative_G_score_06', label_letter='f')

plot_relative_scores(query_lddt_04, 'Query LDDT (threshold 0.4)', 'ROC_query_lddt_04', ylabel='Relative query LDDT score', label_letter='g')
plot_relative_scores(query_lddt_05, 'Query LDDT (threshold 0.5)', 'ROC_query_lddt_05', ylabel='Relative query LDDT score', label_letter='h')
plot_relative_scores(query_lddt_06, 'Query LDDT (threshold 0.6)', 'ROC_query_lddt_06', ylabel='Relative query LDDT score', label_letter='i')

plot_relative_scores(aln_lddt_04, 'Alignment LDDT (threshold 0.4)', 'ROC_aln_lddt_04', ylabel='Relative alignment LDDT score', label_letter='j')
plot_relative_scores(aln_lddt_05, 'Alignment LDDT (threshold 0.5)', 'ROC_aln_lddt_05', ylabel='Relative alignment LDDT score', label_letter='k')
plot_relative_scores(aln_lddt_06, 'Alignment LDDT (threshold 0.6)', 'ROC_aln_lddt_06', ylabel='Relative alignment LDDT score', label_letter='l')


def compute_aln_quality_scores(alignments, alignment_names, mode, thresholds=[0.4, 0.5, 0.6]):
    results = {}
    for thr in thresholds:
        scores_dict = {}
        for df, name in zip(alignments, alignment_names):
            df_opt = df[['query', 'target', 'lddt_all_filtered', 'query_lddt_filter', 'G_score_lddt_filter']].copy()
            df_opt.drop(df_opt[df_opt['query'] == df_opt['target']].index, inplace=True)
            col_name = display_names[name]
            scores_dict[col_name] = calculate_query_aln_quality(df_opt, mode=mode, threshold_value=thr)
        result_df = pd.DataFrame(scores_dict).reset_index().rename(columns={'index': 'query'})
        result_df.fillna(0, inplace=True)
        results[thr] = result_df
    return results

aln_quality_results = compute_aln_quality_scores(alignments, alignment_names, 'G_score_lddt_filter')

plot_relative_scores(aln_quality_results[0.4], 'Alignment quality G-score (threshold 0.4)', 'aln_quality_G_score_04', ylabel='Relative alignment quality', label_letter='m')
plot_relative_scores(aln_quality_results[0.5], 'Alignment quality G-score (threshold 0.5)', 'aln_quality_G_score_05', ylabel='Relative alignment quality', label_letter='n')
plot_relative_scores(aln_quality_results[0.6], 'Alignment quality G-score (threshold 0.6)', 'aln_quality_G_score_06', ylabel='Relative alignment quality', label_letter='o')

def compute_aln_quality_scores(alignments, alignment_names, mode, thresholds=[0.4, 0.5, 0.6]):
    results = {}
    for thr in thresholds:
        scores_dict = {}
        for df, name in zip(alignments, alignment_names):
            df_opt = df[['query', 'target', 'lddt_all_filtered', 'query_lddt_filter', 'G_score_lddt_filter']].copy()
            df_opt.drop(df_opt[df_opt['query'] == df_opt['target']].index, inplace=True)
            col_name = display_names[name]
            scores_dict[col_name] = calculate_query_aln_quality(df_opt, mode=mode, threshold_value=thr)
        result_df = pd.DataFrame(scores_dict).reset_index().rename(columns={'index': 'query'})
        result_df.fillna(0, inplace=True)
        results[thr] = result_df
    return results

aln_quality_results = compute_aln_quality_scores(alignments, alignment_names, 'G_score_lddt_filter')

plot_relative_scores(aln_quality_results[0.4], 'Alignment quality G-score (threshold 0.4)', 'aln_quality_G_score_04', ylabel='Relative alignment quality', label_letter='m')
plot_relative_scores(aln_quality_results[0.5], 'Alignment quality G-score (threshold 0.5)', 'aln_quality_G_score_05', ylabel='Relative alignment quality', label_letter='n')
plot_relative_scores(aln_quality_results[0.6], 'Alignment quality G-score (threshold 0.6)', 'aln_quality_G_score_06', ylabel='Relative alignment quality', label_letter='o')

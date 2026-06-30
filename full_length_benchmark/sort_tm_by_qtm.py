import pandas as pd

df = pd.read_csv('/cbscratch/lolalignBenchmark/=tm_calc_full_low_plddt.tsv', sep='\t')
print(f"Loaded {len(df)} rows")

# Sort by query (preserve order), then by qtmscore descending within each query
df_sorted = df.sort_values(['query', 'qtmscore'], ascending=[True, False], kind='stable')

# But we need to preserve the original query appearance order
# Actually, re-sort: keep query groups together in original order, sort within group by qtmscore
query_order = df['query'].drop_duplicates().tolist()
df_sorted = pd.concat([df[df['query'] == q].sort_values('qtmscore', ascending=False) for q in query_order])

df_sorted.to_csv('/cbscratch/lolalignBenchmark/tm_calc_full_low_plddt_qtm.tsv', sep='\t', index=False)
print(f"Saved sorted file with {len(df_sorted)} rows")

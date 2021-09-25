# coding: utf-8
import pandas as pd
df = pd.read_csv('sgrna_counts.tsv', sep='\t')
df.head()
get_ipython().system('head sgrna_counts_dropped5717dups.tsv')
df.duplicated
df.duplicated().head()
df.duplicated().value_counts()
df.head()
df['sgRNA'].duplicated().value_counts()
df[df['sgRNA'].duplicated()].head()
f = dict.fromkeys(df, 'mean')
f.update(
    dict.fromkeys(df.columns[df.dtypes.eq(object)], 'first'))
f = dict.fromkeys(df, 'first')
f.update(
    dict.fromkeys(['gene'], lambda g: g.iloc[0] if len(g) == 1 else f'multiKO_{"_".join(g)}'))
f
f.groupby(sgRNA).agg(f).head()
df.groupby(sgRNA).agg(f).head()
df.groupby('sgRNA').agg(f).head()
df.groupby('sgRNA').agg(f).shape
df.shape
df.groupby('sgRNA').agg(f).loc['AAAAATGCACACACTTCACT']
df.groupby('sgRNA').agg(f).to_csv('./sgrna_counts_merged.tsv', sep='\t', index=False)

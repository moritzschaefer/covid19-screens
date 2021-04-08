# coding: utf-8
df = pd.read_csv('cs_sgrnas.csv')
import pandas as pd
df = pd.read_csv('cs_sgrnas.csv')
df.head()
df.groupby('grna.sequence.without.PAM').agg({'gene_name': lambda g: f'dKO_{"_".join(g)}', 'SeqID': lambda v: v[0]}).head()
df.groupby('grna.sequence.without.PAM').agg({'gene_name': lambda g: f'dKO_{"_".join(g)}', 'SeqID': len}).head()
df.groupby('grna.sequence.without.PAM').agg({'gene_name': lambda g: f'dKO_{"_".join(g)}', 'SeqID': len})['SeqID'].value_counts()
df.groupby('grna.sequence.without.PAM').agg({'gene_name': lambda g: g.iloc[0] if len(g) == 1 else f'multiKO_{"_".join(g)}', 'SeqID': g.iloc[0]}).head()
df.groupby('grna.sequence.without.PAM').agg({'gene_name': lambda g: g.iloc[0] if len(g) == 1 else f'multiKO_{"_".join(g)}', 'SeqID': lambda g: g.iloc[0]}).head()
df.groupby('grna.sequence.without.PAM').agg({'gene_name': lambda g: g.iloc[0] if len(g) == 1 else f'multiKO_{"_".join(g)}', 'SeqID': lambda g: g.iloc[0]}).reset_index().head()
df.head()
df.groupby('grna.sequence.without.PAM').agg({'gene_name': lambda g: g.iloc[0] if len(g) == 1 else f'multiKO_{"_".join(g)}', 'SeqID': lambda g: g.iloc[0]}).reset_index()[['SeqID', 'grna.sequence.without.PAM', 'gene_name']].to_csv('./cs_sgrnas_grouped_kos.csv', index=False)
get_ipython().run_line_magic('save', 'group_duplicates.py 1-13')

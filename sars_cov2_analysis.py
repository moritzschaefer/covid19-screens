# ** Code block **
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LogNorm
# ** Code block **
# ** Code block **
from matplotlib_venn import venn2, venn3

from moritzsphd.util import remote_file

df_meta = pd.read_csv(f'{ATT_DIR}/genesets.csv')
df_meta['selection'] = df_meta['Title'].apply(lambda t: 'negative' if 'negative' in t.lower() else 'positive')
df_meta['author'] = df_meta[['Title', 'Detailed Title']].astype(str).agg(','.join, axis=1).str.extract('(\w+) et al')
df_meta.loc[[6, 8, 9], 'Organism'] = 'Human'
df_meta.loc[[6, 8, 9], 'Cell/Tissue'] = 'A549-ACE2 cells'
df_meta['#(reported hits)'] = df_meta['Gene Set'].apply(lambda v: len(v.split(',')))
df_meta.loc[0, 'selection'] = 'positive'  # there were reported differently from the others (ace2 appeared in the negatives..)
df_meta.loc[1, 'selection'] = 'negative'

df_meta['label'] = df_meta.reset_index()[['index', 'author', 'Cell/Tissue']].astype(str).agg('_'.join, axis=1)
def _extract_temperature(v):
    try:
        return int(re.search('(\d\d) deg C', v).groups()[0])
    except AttributeError:
        return 37
df_meta['temperature'] = df_meta['Title'].apply(_extract_temperature)
df_meta[['author', 'Cell/Tissue', 'selection', 'temperature', '#(reported hits)']]

# ** Code block **
def split_genes(gs):
    return set([v.strip(' ') for v in gs.split(',')])
# overlap heatmap
def score(gs1, gs2):
    '''
    Return the mean of gs1 in gs2 and gs2 in gs1
    '''
    gs1 = split_genes(gs1)
    gs2 = split_genes(gs2)
    # return len(gs1&gs2) / len(gs1|gs2)

    return int(100 * ((len(gs1 & gs2) / len(gs1)) + (len(gs1 & gs2) / len(gs2))) / 2)

def overlap_count(gs1, gs2):
    gs1 = split_genes(gs1)
    gs2 = split_genes(gs2)
    return len(gs1&gs2)

pos_index = df_meta.index[df_meta['selection'] == 'positive']
neg_index = df_meta.index[df_meta['selection'] == 'negative']

def plot_heatmap(index, title, score_fn=score, annot=False, ax=None):
    if not ax:
        fig, ax = plt.subplots(figsize=(7.5, 6))
    plot_df_meta = pd.DataFrame(
        {df_meta.loc[i, 'label']: [score_fn(df_meta.loc[i, 'Gene Set'], df_meta.loc[i2, 'Gene Set'])
             for i2 in index] for i in index
        },
        index=df_meta.loc[index, 'label'])
    sns.heatmap(plot_df_meta, annot=annot, fmt='d', ax=ax) # , norm=LogNorm(vmin=1, vmax=plot_df_meta.max().max()))
    ax.set_title(title)
    ax.set_ylabel('')

fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
plot_heatmap(pos_index, 'Overlap % of positively selected genes in sars-cov2 screens', annot=True, ax=axes[0])
plot_heatmap(pos_index, 'Overlap count of positively selected genes in sars-cov2 screens', ax=axes[1], annot=True, score_fn=overlap_count)

# get a consensus per cell line (genes that come up in at least 1/3 of data sets)

# venn diagram with our analyses to highlight overlapping genes (for convidence) and specific genes (for novelty).

# ** Code block **
fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
plot_heatmap(neg_index, 'Overlap % of negatively selected genes in sars-cov2 screens', ax=axes[0], annot=True)
plot_heatmap(neg_index, 'Overlap count of negatively selected genes in sars-cov2 screens', ax=axes[1], annot=True, score_fn=overlap_count)

genecount = defaultdict(int)
def _count(gs):
    for g in gs.split(','):
        genecount[g.strip(' ')] += 1
        
_ = df_meta.loc[pos_index, 'Gene Set'].apply(_count)
','.join(pd.Series(genecount).sort_values(ascending=False).iloc[:20].index)

# ** Code block **
genecount = defaultdict(int)
_ = df_meta.loc[neg_index, 'Gene Set'].apply(_count)
s = pd.Series(genecount).sort_values(ascending=False)
s[s >= 2]

# ** Code block **
datasets = [('daniloski', 'daniloski_moi0.01_rra'),
            ('daniloski', 'daniloski_moi0.3_rra'),
            ('hoffmann', 'hoffmann_33deg_rra'),
            ('hoffmann', 'hoffmann_37deg_rra'),
            ('schneider', 'schneider_33deg_rra'),
            ('schneider', 'schneider_37deg_rra'),
            ('wang', 'wang_cov_combined'),
            ]

def plot_datasets(datasets, selection, field, reload_file=False):
    '''
    field: One of fdr, p-value
    '''

    fig, axes = plt.subplots(len(datasets), 2, figsize=(12, 32))
    for dataset, (ax0, ax1) in zip(datasets, axes):
        temp = 33 if '33' in dataset[1] else 37
        try:
            # get reported hits
            dataset_meta = df_meta.set_index(['author', 'temperature', 'selection']).loc[(dataset[0].title(), temp, selection)]
            dataset_genes = [v.strip(' ') for v in dataset_meta['Gene Set'].split(',')]

            # get our hits
            mageck_df = pd.read_csv(
                remote_file(f'cclab:~/sars/{dataset[0]}/workflow/results/test/{dataset[1]}.gene_summary.txt', reload=reload_file),
                sep='\t'
            )
            test_col = f'{selection[:3]}|{field}'

            sns.distplot(mageck_df.loc[mageck_df['id'].isin(dataset_genes), test_col], ax=ax0)
            # draw venn diagram with overlap: 0.05, 0.1 and reported hits. also add all our IDs to see if they all exist
            venn3((set(dataset_genes), set(mageck_df.loc[mageck_df[test_col] < 0.1, 'id']), set(mageck_df.loc[mageck_df[test_col] < 0.05, 'id'])), ('Reported top-hits', f'Our analysis ({field} < 0.1)', f'Our analysis ({field} < 0.05)'), ax=ax1)

            ax0.set_xlabel('Distribution of FDR from our analysis for the reported top picks')
            ax0.set_ylabel(dataset[1])
        except KeyError:
            pass

# ** Code block **
plot_datasets(datasets, 'positive', 'fdr', False)

# ** Code block **
plot_datasets(datasets, 'positive', 'p-value')


def load_dataset(dataset):
    mageck_df = pd.read_csv(
        remote_file(f'cclab:~/sars/{dataset[0]}/workflow/results/test/{dataset[1]}.gene_summary.txt'),
        sep='\t'
    )
    return mageck_df 


def select_by_threshold(df, field, threshold):
    return set(df[df[field] < threshold]['id'])

def select_top_n(df, field, n, ascending=True):
    return set(df.sort_values(field, ascending=ascending).iloc[:n]['id'])
    
hoffmann33 = datasets[2]
hoffmann37 = datasets[3]
schneider33 = datasets[4]
schneider37 = datasets[5]

# ** Code block **
fig, axes = plt.subplots(2, 2, figsize=(8, 7))
THRESHOLD = 0.1
# TODO threshold selection
for y, x, seta, setb, labels, field in [
    (0, 0, hoffmann33, hoffmann37, ('hoffmann33', 'hoffmann37'), 'pos|fdr'),
    (0, 1, schneider33, schneider37, ('schneider33', 'schneider37'), 'pos|fdr'),
    (1, 0, hoffmann33, schneider33, ('hoffmann33', 'schneider37'), 'pos|fdr'),
    (1, 1, hoffmann37, schneider37, ('hoffmann33', 'schneider37'), 'pos|fdr'),
]:
    venn2((select_by_threshold(load_dataset(seta), field, THRESHOLD),
          select_by_threshold(load_dataset(setb), field, THRESHOLD)),
          labels, ax=axes[y, x])
    axes[y, x].set_title(f'{field} < {THRESHOLD}')

# ** Code block **
# get our hits
def gene_set_matrix(datasets, selector_function, field, threshold):
    gene_sets = {}
    for dataset in datasets:
        gene_sets[dataset[1]] = selector_function(load_dataset(dataset), field, threshold)
    index = pd.Index(set.union(*gene_sets.values()))
    sub_df = pd.DataFrame(index=index, data={key: [gene in gene_set for gene in index] for key, gene_set in gene_sets.items()})
    return sub_df

dfs = []
threshold = 0.01
field = 'pos|fdr'
selector_function = select_by_threshold
for threshold, field, selector_function in [
    (0.001, 'pos|fdr', select_by_threshold), (0.01, 'pos|fdr', select_by_threshold), (0.1, 'pos|fdr', select_by_threshold),
    (20, 'pos|fdr', select_top_n), (50, 'pos|fdr', select_top_n), (100, 'pos|fdr', select_top_n), (200, 'pos|fdr', select_top_n),
    
    (0.001, 'pos|p-value', select_by_threshold), (0.01, 'pos|p-value', select_by_threshold), (0.1, 'pos|p-value', select_by_threshold),
    (20, 'pos|p-value', select_top_n), (50, 'pos|p-value', select_top_n), (100, 'pos|p-value', select_top_n), (200, 'pos|p-value', select_top_n),
    ]:
    sub_df = gene_set_matrix(datasets, selector_function, field, threshold)
    sub_df = sub_df.sum(axis=1).value_counts().sort_index().to_frame().reset_index()
    sub_df.columns = ['dataset_count', 'overlap_count']
    # sub_df.index.name = 'dataset_count'
    if threshold < 1.0:
        sub_df['selection'] = f'Selected by {field} < {threshold}'
    else:
        sub_df['selection'] = f'Selected top {threshold} {field} hits'
    dfs.append(sub_df)
plot_df = pd.concat(dfs)

# ** Code block **
sns.barplot(data=plot_df[plot_df['selection'].str.contains('pos\|p-value')], x='selection', y='overlap_count', hue='dataset_count')
plt.xticks(rotation=30, ha='right')
plt.yscale('log')
plt.title('Number of targets found in 1-6 of the analyzed datasets')

# ** Code block **
gene_set1 = gene_set_matrix(datasets, select_by_threshold, 'pos|p-value', 0.01).sum(axis=1)
gene_set1.name = '#(datasets with p-value < 0.01)'
gene_set1[gene_set1 >= 4].sort_values(ascending=False)
# gene_set2 = gene_set_matrix(datasets, select_by_threshold, 'pos|p-value', 0.1)


# gene_set1 = gene_set1[gene_set1 >= 4]
# gene_set3 = gene_set2.index[gene_set2.sum(axis=1) >= 5]
# gene_set2 = gene_set2.index[gene_set2.sum(axis=1) >= 4]

# ** Code block **
gene_set2 = gene_set_matrix(datasets, select_by_threshold, 'pos|p-value', 0.1).sum(axis=1)
gene_set2.name = '#(datasets with p-value < 0.1)'
gene_set2[gene_set2 >= 4].sort_values(ascending=False)
# gene_set2 = gene_set_matrix(datasets, select_by_threshold, 'pos|p-value', 0.1)


# gene_set1 = gene_set1[gene_set1 >= 4]
# gene_set3 = gene_set2.index[gene_set2.sum(axis=1) >= 5]
# gene_set2 = gene_set2.index[gene_set2.sum(axis=1) >= 4]

# ** Code block **
our_datasets = [
    ('brunello1', 'HEK_CRISPRko1_rra_cov_1st'),
    ('brunello1', 'HEK_CRISPRko1_rra_cov_2nd_0.01'),
    ('brunello1', 'HEK_CRISPRko1_rra_cov_2nd_0.1'),
    ('brunello1', 'HEK_CRIPSRko1_rra_cov_2nd_ext_0.01'),
    ('brunello1', 'HEK_CRIPSRko1_rra_cov_2nd_ext_0.1'),
    ('brunello2', 'HEK_CRISPRko2_rra_2nd'),
    ('brunello2', 'HEK_CRISPRko2_rra_1st')
]

def our_selection(selection_function, threshold, field, experiments=['brunello1', 'brunello2']):
    return set.union(*[
        selection_function(load_dataset(dataset), field, threshold)
        for dataset in our_datasets
        if dataset[0] in experiments
    ])

# ** Code block **
our_genes = our_selection(select_by_threshold, 0.05, 'pos|fdr', ['brunello1'])  # since we take the union of hits, we have to use the FDR (instead of the p-value). Because there is no increase in statistical significance by combining sets

venn3((our_genes, set(gene_set1.index[gene_set1 >= 3]), set(gene_set2.index[gene_set2 >= 3])), ('Our positive selection (FDR<0.05)', 'p-value < 0.01 in at least 3 publ. datasets', 'p-value < 0.1 in at least 3 publ. datasets'))

# ** Code block **
','.join(set(gene_set2.index[gene_set2 >= 3]).intersection(our_genes))

# ** Code block **
','.join(our_genes - set(gene_set_matrix(datasets, select_by_threshold, 'pos|p-value', 0.1).sum(axis=1).index))

# ** Code block **
','.join(our_genes)

# ** Code block **
our_genes = our_selection(select_by_threshold, 0.05, 'pos|fdr', ['brunello2'])  # since we take the union of hits, we have to use the FDR (instead of the p-value). Because there is no increase in statistical significance by combining sets

venn3((our_genes, set(gene_set1.index[gene_set1 >= 3]), set(gene_set2.index[gene_set2 >= 3])), ('Our positive selection (FDR<0.05)', 'p-value < 0.01 in at least 3 publ. datasets', 'p-value < 0.1 in at least 3 publ. datasets'))

# ** Code block **
','.join(set(gene_set2.index[gene_set2 >= 3]).intersection(our_genes))

# ** Code block **
','.join(our_genes - set(gene_set_matrix(datasets, select_by_threshold, 'pos|p-value', 0.1).sum(axis=1).index))

# ** Code block **
','.join(our_genes)

# ** Code block **
our_genes = our_selection(select_by_threshold, 0.05, 'pos|fdr')  # since we take the union of hits, we have to use the FDR (instead of the p-value). Because there is no increase in statistical significance by combining sets

venn3((our_genes, set(gene_set1.index[gene_set1 >= 3]), set(gene_set2.index[gene_set2 >= 3])), ('Our positive selection (FDR<0.05)', 'p-value < 0.01 in at least 3 publ. datasets', 'p-value < 0.1 in at least 3 publ. datasets'))

# ** Code block **
','.join(set(gene_set2.index[gene_set2 >= 3]).intersection(our_genes))

# ** Code block **
','.join(our_genes - set(gene_set_matrix(datasets, select_by_threshold, 'pos|p-value', 0.1).sum(axis=1).index))

# ** Code block **
','.join(our_genes)

# ** Code block **
datasets = [('daniloski', 'daniloski_moi0.01_rra'),
            ('daniloski', 'daniloski_moi0.3_rra'),
            ('hoffmann', 'hoffmann_33deg_rra'),
            ('hoffmann', 'hoffmann_37deg_rra'),
            ('schneider', 'schneider_33deg_rra'),
            ('schneider', 'schneider_37deg_rra'),
            ('wang', 'wang_cov_combined'),
            ]

# ** Code block **
plot_datasets(datasets, 'negative', 'fdr', False)

# ** Code block **
plot_datasets(datasets, 'negative', 'p-value')


def load_dataset(dataset):
    mageck_df = pd.read_csv(
        remote_file(f'cclab:~/sars/{dataset[0]}/workflow/results/test/{dataset[1]}.gene_summary.txt'),
        sep='\t'
    )
    return mageck_df 


def select_by_threshold(df, field, threshold):
    return set(df[df[field] < threshold]['id'])

def select_top_n(df, field, n, ascending=True):
    return set(df.sort_values(field, ascending=ascending).iloc[:n]['id'])
    
hoffmann33 = datasets[2]
hoffmann37 = datasets[3]
schneider33 = datasets[4]
schneider37 = datasets[5]

# ** Code block **
fig, axes = plt.subplots(2, 2, figsize=(8, 7))
THRESHOLD = 0.1
# TODO threshold selection
for y, x, seta, setb, labels, field in [
    (0, 0, hoffmann33, hoffmann37, ('hoffmann33', 'hoffmann37'), 'neg|fdr'),
    (0, 1, schneider33, schneider37, ('schneider33', 'schneider37'), 'neg|fdr'),
    (1, 0, hoffmann33, schneider33, ('hoffmann33', 'schneider37'), 'neg|fdr'),
    (1, 1, hoffmann37, schneider37, ('hoffmann33', 'schneider37'), 'neg|fdr'),
]:
    venn2((select_by_threshold(load_dataset(seta), field, THRESHOLD),
          select_by_threshold(load_dataset(setb), field, THRESHOLD)),
          labels, ax=axes[y, x])
    axes[y, x].set_title(f'{field} < {THRESHOLD}')

# ** Code block **
# get our hits
def gene_set_matrix(datasets, selector_function, field, threshold):
    gene_sets = {}
    for dataset in datasets:
        gene_sets[dataset[1]] = selector_function(load_dataset(dataset), field, threshold)
    index = pd.Index(set.union(*gene_sets.values()))
    sub_df = pd.DataFrame(index=index, data={key: [gene in gene_set for gene in index] for key, gene_set in gene_sets.items()})
    return sub_df

dfs = []
threshold = 0.01
field = 'neg|fdr'
selector_function = select_by_threshold
for threshold, field, selector_function in [
    (0.001, 'neg|fdr', select_by_threshold), (0.01, 'neg|fdr', select_by_threshold), (0.1, 'neg|fdr', select_by_threshold),
    (20, 'neg|fdr', select_top_n), (50, 'neg|fdr', select_top_n), (100, 'neg|fdr', select_top_n), (200, 'neg|fdr', select_top_n),
    
    (0.001, 'neg|p-value', select_by_threshold), (0.01, 'neg|p-value', select_by_threshold), (0.1, 'neg|p-value', select_by_threshold),
    (20, 'neg|p-value', select_top_n), (50, 'neg|p-value', select_top_n), (100, 'neg|p-value', select_top_n), (200, 'neg|p-value', select_top_n),
    ]:
    sub_df = gene_set_matrix(datasets, selector_function, field, threshold)
    sub_df = sub_df.sum(axis=1).value_counts().sort_index().to_frame().reset_index()
    sub_df.columns = ['dataset_count', 'overlap_count']
    # sub_df.index.name = 'dataset_count'
    if threshold < 1.0:
        sub_df['selection'] = f'Selected by {field} < {threshold}'
    else:
        sub_df['selection'] = f'Selected top {threshold} {field} hits'
    dfs.append(sub_df)
plot_df = pd.concat(dfs)

# ** Code block **
sns.barplot(data=plot_df[plot_df['selection'].str.contains('neg\|p-value')], x='selection', y='overlap_count', hue='dataset_count')
plt.xticks(rotation=30, ha='right')
plt.yscale('log')
plt.title('Number of targets found in 1-6 of the analyzed datasets')

# ** Code block **
gene_set1 = gene_set_matrix(datasets, select_by_threshold, 'neg|p-value', 0.01).sum(axis=1)
gene_set1.name = '#(datasets with p-value < 0.01)'
gene_set1[gene_set1 >= 3].sort_values(ascending=False)
# gene_set2 = gene_set_matrix(datasets, select_by_threshold, 'neg|p-value', 0.1)


# gene_set1 = gene_set1[gene_set1 >= 4]
# gene_set3 = gene_set2.index[gene_set2.sum(axis=1) >= 5]
# gene_set2 = gene_set2.index[gene_set2.sum(axis=1) >= 4]

# ** Code block **
gene_set2 = gene_set_matrix(datasets, select_by_threshold, 'neg|p-value', 0.1).sum(axis=1)
gene_set2.name = '#(datasets with p-value < 0.1)'
gene_set2[gene_set2 >= 4].sort_values(ascending=False)
# gene_set2 = gene_set_matrix(datasets, select_by_threshold, 'neg|p-value', 0.1)


# gene_set1 = gene_set1[gene_set1 >= 4]
# gene_set3 = gene_set2.index[gene_set2.sum(axis=1) >= 5]
# gene_set2 = gene_set2.index[gene_set2.sum(axis=1) >= 4]

# ** Code block **
our_genes = our_selection(select_by_threshold, 0.01, 'neg|fdr')  # since we take the union of hits, we have to use the FDR (instead of the p-value). Because there is no increase in statistical significance by combining sets

venn3((our_genes, set(gene_set1.index[gene_set1 >= 3]), set(gene_set2.index[gene_set2 >= 3])), ('Our negative selection (FDR<0.01)', 'p-value < 0.01 in at least 3 publ. datasets', 'p-value < 0.1 in at least 3 publ. datasets'))

# ** Code block **
','.join(set(gene_set2.index[gene_set2 >= 3].intersection(our_genes)))

# ** Code block **
','.join(our_genes - set(gene_set_matrix(datasets, select_by_threshold, 'neg|p-value', 0.1).sum(axis=1).index))

# ** Code block **
','.join(our_genes)

# ** Code block **
datasets = [('hoffmann', 'hoffmann_33deg_oc43_rra'),
            ('schneider', 'schneider_33deg_oc43_rra'),
            ('wang', 'wang_oc43'),
            ]

# ** Code block **
plot_datasets(datasets, 'positive', 'fdr', False)

# ** Code block **
plot_datasets(datasets, 'positive', 'p-value')

# ** Code block **
# get our hits
def gene_set_matrix(datasets, selector_function, field, threshold):
    gene_sets = {}
    for dataset in datasets:
        gene_sets[dataset[1]] = selector_function(load_dataset(dataset), field, threshold)
    index = pd.Index(set.union(*gene_sets.values()))
    sub_df = pd.DataFrame(index=index, data={key: [gene in gene_set for gene in index] for key, gene_set in gene_sets.items()})
    return sub_df

dfs = []
threshold = 0.01
field = 'pos|fdr'
selector_function = select_by_threshold
for threshold, field, selector_function in [
    (0.001, 'pos|fdr', select_by_threshold), (0.01, 'pos|fdr', select_by_threshold), (0.1, 'pos|fdr', select_by_threshold),
    (20, 'pos|fdr', select_top_n), (50, 'pos|fdr', select_top_n), (100, 'pos|fdr', select_top_n), (200, 'pos|fdr', select_top_n),
    
    (0.001, 'pos|p-value', select_by_threshold), (0.01, 'pos|p-value', select_by_threshold), (0.1, 'pos|p-value', select_by_threshold),
    (20, 'pos|p-value', select_top_n), (50, 'pos|p-value', select_top_n), (100, 'pos|p-value', select_top_n), (200, 'pos|p-value', select_top_n),
    ]:
    sub_df = gene_set_matrix(datasets, selector_function, field, threshold)
    sub_df = sub_df.sum(axis=1).value_counts().sort_index().to_frame().reset_index()
    sub_df.columns = ['dataset_count', 'overlap_count']
    # sub_df.index.name = 'dataset_count'
    if threshold < 1.0:
        sub_df['selection'] = f'Selected by {field} < {threshold}'
    else:
        sub_df['selection'] = f'Selected top {threshold} {field} hits'
    dfs.append(sub_df)
plot_df = pd.concat(dfs)

# ** Code block **
sns.barplot(data=plot_df[plot_df['selection'].str.contains('pos\|p-value')], x='selection', y='overlap_count', hue='dataset_count')
plt.xticks(rotation=30, ha='right')
plt.yscale('log')
plt.title('Number of targets found in 1-6 of the analyzed datasets')

# ** Code block **
gene_set1 = gene_set_matrix(datasets, select_by_threshold, 'pos|p-value', 0.01).sum(axis=1)
gene_set1.name = '#(datasets with p-value < 0.01)'
gene_set1[gene_set1 >= 2].sort_values(ascending=False)
# gene_set2 = gene_set_matrix(datasets, select_by_threshold, 'pos|p-value', 0.1)


# gene_set1 = gene_set1[gene_set1 >= 4]
# gene_set3 = gene_set2.index[gene_set2.sum(axis=1) >= 5]
# gene_set2 = gene_set2.index[gene_set2.sum(axis=1) >= 4]

# ** Code block **
gene_set2 = gene_set_matrix(datasets, select_by_threshold, 'pos|p-value', 0.1).sum(axis=1)
gene_set2.name = '#(datasets with p-value < 0.1)'
gene_set2[gene_set2 >= 3].sort_values(ascending=False)
# gene_set2 = gene_set_matrix(datasets, select_by_threshold, 'pos|p-value', 0.1)


# gene_set1 = gene_set1[gene_set1 >= 4]
# gene_set3 = gene_set2.index[gene_set2.sum(axis=1) >= 5]
# gene_set2 = gene_set2.index[gene_set2.sum(axis=1) >= 4]

# ** Code block **
our_datasets = [
    ('brunello1', 'HEK_CRIPSRko1_rra_oc43_1st'),
    ('brunello1', 'HEK_CRIPSRko1_rra_oc43_2nd_0.1'),
    ('brunello1', 'HEK_CRIPSRko1_rra_oc43_2nd_0.01'),
]

def our_selection(selection_function, threshold, field, experiments=['brunello1', 'brunello2']):
    return set.union(*[
        selection_function(load_dataset(dataset), field, threshold)
        for dataset in our_datasets
        if dataset[0] in experiments
    ])

# ** Code block **
our_genes = our_selection(select_by_threshold, 0.05, 'pos|fdr', ['brunello1'])  # since we take the union of hits, we have to use the FDR (instead of the p-value). Because there is no increase in statistical significance by combining sets

venn3((our_genes, set(gene_set1.index[gene_set1 >= 2]), set(gene_set2.index[gene_set2 >= 3])), ('Our positive selection (FDR<0.05)', 'p-value < 0.01 in at least 3 publ. datasets', 'p-value < 0.1 in at least 3 publ. datasets'))

# ** Code block **
','.join(set(gene_set1.index[gene_set1 >= 2]).intersection(our_genes))

# ** Code block **
','.join(our_genes - set(gene_set_matrix(datasets, select_by_threshold, 'pos|p-value', 0.1).sum(axis=1).index))

# ** Code block **
','.join(our_genes)

# ** Code block **
datasets = [('hoffmann', 'hoffmann_33deg_oc43_rra'),
            ('schneider', 'schneider_33deg_oc43_rra'),
            ('wang', 'wang_oc43'),
            ]

# ** Code block **
plot_datasets(datasets, 'negative', 'fdr', False)

# ** Code block **
plot_datasets(datasets, 'negative', 'p-value')

# ** Code block **
# get our hits
def gene_set_matrix(datasets, selector_function, field, threshold):
    gene_sets = {}
    for dataset in datasets:
        gene_sets[dataset[1]] = selector_function(load_dataset(dataset), field, threshold)
    index = pd.Index(set.union(*gene_sets.values()))
    sub_df = pd.DataFrame(index=index, data={key: [gene in gene_set for gene in index] for key, gene_set in gene_sets.items()})
    return sub_df

dfs = []
threshold = 0.01
field = 'neg|fdr'
selector_function = select_by_threshold
for threshold, field, selector_function in [
    (0.001, 'neg|fdr', select_by_threshold), (0.01, 'neg|fdr', select_by_threshold), (0.1, 'neg|fdr', select_by_threshold),
    (20, 'neg|fdr', select_top_n), (50, 'neg|fdr', select_top_n), (100, 'neg|fdr', select_top_n), (200, 'neg|fdr', select_top_n),
    
    (0.001, 'neg|p-value', select_by_threshold), (0.01, 'neg|p-value', select_by_threshold), (0.1, 'neg|p-value', select_by_threshold),
    (20, 'neg|p-value', select_top_n), (50, 'neg|p-value', select_top_n), (100, 'neg|p-value', select_top_n), (200, 'neg|p-value', select_top_n),
    ]:
    sub_df = gene_set_matrix(datasets, selector_function, field, threshold)
    sub_df = sub_df.sum(axis=1).value_counts().sort_index().to_frame().reset_index()
    sub_df.columns = ['dataset_count', 'overlap_count']
    # sub_df.index.name = 'dataset_count'
    if threshold < 1.0:
        sub_df['selection'] = f'Selected by {field} < {threshold}'
    else:
        sub_df['selection'] = f'Selected top {threshold} {field} hits'
    dfs.append(sub_df)
plot_df = pd.concat(dfs)

# ** Code block **
sns.barplot(data=plot_df[plot_df['selection'].str.contains('neg\|p-value')], x='selection', y='overlap_count', hue='dataset_count')
plt.xticks(rotation=30, ha='right')
plt.yscale('log')
plt.title('Number of targets found in 1-6 of the analyzed datasets')

# ** Code block **
gene_set1 = gene_set_matrix(datasets, select_by_threshold, 'neg|p-value', 0.01).sum(axis=1)
gene_set1.name = '#(datasets with p-value < 0.01)'
gene_set1[gene_set1 >= 2].sort_values(ascending=False)
# gene_set2 = gene_set_matrix(datasets, select_by_threshold, 'neg|p-value', 0.1)


# gene_set1 = gene_set1[gene_set1 >= 4]
# gene_set3 = gene_set2.index[gene_set2.sum(axis=1) >= 5]
# gene_set2 = gene_set2.index[gene_set2.sum(axis=1) >= 4]

# ** Code block **
gene_set2 = gene_set_matrix(datasets, select_by_threshold, 'neg|p-value', 0.1).sum(axis=1)
gene_set2.name = '#(datasets with p-value < 0.1)'
gene_set2[gene_set2 >= 2].sort_values(ascending=False)
# gene_set2 = gene_set_matrix(datasets, select_by_threshold, 'neg|p-value', 0.1)


# gene_set1 = gene_set1[gene_set1 >= 4]
# gene_set3 = gene_set2.index[gene_set2.sum(axis=1) >= 5]
# gene_set2 = gene_set2.index[gene_set2.sum(axis=1) >= 4]

# ** Code block **
our_datasets = [
    ('brunello1', 'HEK_CRIPSRko1_rra_oc43_1st'),
    ('brunello1', 'HEK_CRIPSRko1_rra_oc43_2nd_0.1'),
    ('brunello1', 'HEK_CRIPSRko1_rra_oc43_2nd_0.01'),
]

def our_selection(selection_function, threshold, field, experiments=['brunello1', 'brunello2']):
    return set.union(*[
        selection_function(load_dataset(dataset), field, threshold)
        for dataset in our_datasets
        if dataset[0] in experiments
    ])

# ** Code block **
our_genes = our_selection(select_by_threshold, 0.05, 'neg|fdr', ['brunello1'])  # since we take the union of hits, we have to use the FDR (instead of the p-value). Because there is no increase in statistical significance by combining sets

venn3((our_genes, set(gene_set1.index[gene_set1 >= 2]), set(gene_set2.index[gene_set2 >= 2])), ('Our negative selection (FDR<0.05)', 'p-value < 0.01 in at least 3 publ. datasets', 'p-value < 0.1 in at least 3 publ. datasets'))

# ** Code block **
','.join(set(gene_set2.index[gene_set2 >= 2]).intersection(our_genes))

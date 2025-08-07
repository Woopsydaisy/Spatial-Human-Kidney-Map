import csv
import anndata as ad
import gzip
import os
import scipy.io
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from google.colab import drive
import leidenalg as la
from pathlib import Path

##first piecharts of center and surrounding
adata=sc.read_h5ad('/content/drive/MyDrive/Bernhard/revision_nature/revision_object.h5ad')
colors_neighbors={
    'Immune': '#EF785F',
    'Fibroblast': '#d279af',
    'iPT': '#d1c80f',
    'iTAL': '#88bd2c',
    'TAL':'#8d584e',
    'PT':'#28a237',
    'glom':'#67bc93',
    'ec': '#1a79b7',
    'other': '#2B2B2B',
}

adata_subset=adata[adata.obs['iPT_iLOH_ME_20um']!='Unknown'].copy()

import matplotlib.pyplot as plt
import pandas as pd

# Assuming new_adata_subset is your AnnData object
samples = adata_subset.obs['iPT_iLOH_ME']
cluster_labels = adata_subset.obs['annotation_updated']

# Create a DataFrame from the relevant columns
data = pd.DataFrame({'iPT_iLOH_ME': samples, 'annotation_updated': cluster_labels})
data=data[data['annotation_updated']!='Unknown']
# Group by 'sample' and 'cluster_labels_annotated' and count the occurrences
grouped_data = data.groupby(['iPT_iLOH_ME', 'annotation_updated']).size().unstack(fill_value=0)
grouped_data['glom']=grouped_data['EC_glom']+grouped_data['MC1']+grouped_data['PEC']+grouped_data['Podo']
grouped_data=grouped_data.drop(['EC_glom', 'MC1', 'PEC', 'Podo'], axis=1)
grouped_data['ec']=grouped_data['EC_DVR']+grouped_data['EC_Peritub']+grouped_data['VSMC']
grouped_data=grouped_data.drop(['EC_DVR', 'EC_Peritub', 'VSMC'], axis=1)
grouped_data['other']=grouped_data['CNT']+grouped_data['DCT']+grouped_data['DTL_ATL']+grouped_data['EC_Lymph']+grouped_data['IC A']+grouped_data['IC B']+grouped_data['PC']
grouped_data=grouped_data.drop(['CNT', 'DCT', 'DTL_ATL', 'IC A', 'IC B', 'PC', 'EC_Lymph'], axis=1)
grouped_data=grouped_data.drop('Unknown', axis=0)
grouped_data

# Extract the relevant columns
samples = adata_subset.obs['iPT_iLOH_ME_20um']
cluster_labels = adata_subset.obs['annotation_updated']

# Create a DataFrame from the relevant columns
data = pd.DataFrame({'iPT_iLOH_ME_20um': samples, 'annotation_updated': cluster_labels})
# Group by 'sample' and 'cluster_labels_annotated' and count the occurrences
grouped_data_all = data.groupby(['iPT_iLOH_ME_20um', 'annotation_updated']).size().unstack(fill_value=0)
grouped_data_all['glom']=grouped_data_all['EC_glom']+grouped_data_all['MC1']+grouped_data_all['PEC']+grouped_data_all['Podo']
grouped_data_all=grouped_data_all.drop(['EC_glom', 'MC1', 'PEC', 'Podo'], axis=1)
grouped_data_all['ec']=grouped_data_all['EC_DVR']+grouped_data_all['EC_Peritub']+grouped_data_all['VSMC']
grouped_data_all=grouped_data_all.drop(['EC_DVR', 'EC_Peritub', 'VSMC'], axis=1)
grouped_data_all['other']=grouped_data_all['CNT']+grouped_data_all['DCT']+grouped_data_all['DTL_ATL']+grouped_data_all['EC_Lymph']+grouped_data_all['IC A']+grouped_data_all['IC B']+grouped_data_all['PC']
grouped_data_all=grouped_data_all.drop(['CNT', 'DCT', 'DTL_ATL', 'IC A', 'IC B', 'PC', 'EC_Lymph'], axis=1)
grouped_data_all

grouped_data_all = grouped_data_all-grouped_data
grouped_data_all

colors = [colors_neighbors.get(col, 'grey') for col in grouped_data_all.columns]
# Plot pie charts for each subcluster
for subcluster in grouped_data_all.index:
    subcluster_data = grouped_data_all.loc[subcluster]
    print(subcluster)
    # Plot pie chart
    plt.figure(figsize=(6, 6))
    plt.pie(subcluster_data, startangle=0,colors=colors)
    plt.title(f'')
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig(f'Pie Chart {subcluster} outer circle.png', dpi=450, bbox_inches='tight')
    plt.show()



import matplotlib.pyplot as plt
import pandas as pd
inner_circle=adata_subset[adata_subset.obs['iPT_iLOH_ME']!='Unknown'].copy()
# Assuming new_adata_subset is your AnnData object
samples = inner_circle.obs['iPT_iLOH_ME']
cluster_labels = inner_circle.obs['annotation_updated']

# Create a DataFrame from the relevant columns
data = pd.DataFrame({'iPT_iLOH_ME': samples, 'annotation_updated': cluster_labels})
data=data[data['annotation_updated']!='Unknown']
# Group by 'sample' and 'cluster_labels_annotated' and count the occurrences
grouped_data = data.groupby(['iPT_iLOH_ME', 'annotation_updated']).size().unstack(fill_value=0)
grouped_data

colors = [colors_neighbors.get(col, 'grey') for col in grouped_data.columns]
# Plot pie charts for each subcluster
for subcluster in grouped_data.index:
    subcluster_data = grouped_data.loc[subcluster]
    print(subcluster)
    # Plot pie chart
    plt.figure(figsize=(6, 6))
    plt.pie(subcluster_data, startangle=0,colors=colors)
    plt.title(f'')
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig(f'Pie Chart {subcluster} inner circle.png', dpi=450, bbox_inches='tight')
    plt.show()


##second comparisons between celltypes across microenvironments

adata_subset.obs['iPT_iLOH_ME_20um'].value_counts()

# Extract the relevant columns
samples = adata_subset.obs['iPT_iLOH_ME_20um']
cluster_labels = adata_subset.obs['annotation_updated']

# Create a DataFrame from the relevant columns
data = pd.DataFrame({'iPT_iLOH_ME_20um': samples, 'annotation_updated': cluster_labels})
# Group by 'sample' and 'cluster_labels_annotated' and count the occurrences
grouped_data_all = data.groupby(['iPT_iLOH_ME_20um', 'annotation_updated']).size().unstack(fill_value=0)
grouped_data_all['glom']=grouped_data_all['EC_glom']+grouped_data_all['MC1']+grouped_data_all['PEC']+grouped_data_all['Podo']
grouped_data_all=grouped_data_all.drop(['EC_glom', 'MC1', 'PEC', 'Podo'], axis=1)
grouped_data_all['ec']=grouped_data_all['EC_DVR']+grouped_data_all['EC_Peritub']+grouped_data_all['VSMC']
grouped_data_all=grouped_data_all.drop(['EC_DVR', 'EC_Peritub', 'VSMC'], axis=1)
grouped_data_all['other']=grouped_data_all['CNT']+grouped_data_all['DCT']+grouped_data_all['DTL_ATL']+grouped_data_all['EC_Lymph']+grouped_data_all['IC A']+grouped_data_all['IC B']+grouped_data_all['PC']
grouped_data_all=grouped_data_all.drop(['CNT', 'DCT', 'DTL_ATL', 'IC A', 'IC B', 'PC', 'EC_Lymph'], axis=1)
grouped_data_all

samples = adata_subset.obs['iPT_iLOH_ME_20um']
cluster_labels = adata_subset.obs['unique_sample_identifier']
data = pd.DataFrame({'unique_sample_identifier': cluster_labels,'iPT_iLOH_ME_20um': samples  })
grouped_data_all_counts = data.groupby([ 'unique_sample_identifier','iPT_iLOH_ME_20um']).size().unstack(fill_value=0)
grouped_data_all_counts

fibroblasts=adata[adata.obs['annotation_updated']=='Fibroblast']
fibroblasts=fibroblasts[fibroblasts.obs['iPT_iLOH_ME_20um']!='Unknown'].copy()
fibroblasts

samples = fibroblasts.obs['unique_sample_identifier']
cluster_labels = fibroblasts.obs['iPT_iLOH_ME_20um']
data = pd.DataFrame({'iPT_iLOH_ME_20um': cluster_labels, 'unique_sample_identifier': samples})
grouped_data = data.groupby([ 'unique_sample_identifier', 'iPT_iLOH_ME_20um']).size().unstack(fill_value=0)
grouped_data

grouped_data=grouped_data/grouped_data_all_counts
grouped_data

grouped_data=grouped_data.fillna(0)

import seaborn as sns
import matplotlib.pyplot as plt

# Reset index to bring 'unique_sample_identifier' into a column
grouped_data_reset = grouped_data.reset_index()

# Melt to long format
melted = grouped_data_reset.melt(id_vars='unique_sample_identifier', 
                                 var_name='ME', 
                                 value_name='Fraction')

# Make violin plot
plt.figure(figsize=(7, 4))
sns.violinplot(data=melted, x='ME', y='Fraction', inner='box', cut=0)
sns.stripplot(data=melted, x='ME', y='Fraction', color='k', size=2, alpha=0.5)

plt.xticks(rotation=0, ha='center')
plt.title('Fraction of Fibroblasts per ME per Sample')
plt.tight_layout()
plt.savefig('Fraction of Fibroblasts per ME per Sample.png', dpi=900, bbox_inches='tight')
plt.show()


from scipy.stats import mannwhitneyu

results = []
reference_group = 'Profibrotic ME'

# Loop over other MEs
for other_me in melted['ME'].unique():
    if other_me == reference_group:
        continue

    vals_profibrotic = melted[melted['ME'] == reference_group]['Fraction']
    vals_other = melted[melted['ME'] == other_me]['Fraction']

    stat, pval = mannwhitneyu(vals_profibrotic, vals_other, alternative='two-sided')

    results.append({
        'reference_ME': reference_group,
        'compared_ME': other_me,
        'U_statistic': stat,
        'p_value': pval
    })

# Create DataFrame
stats_df = pd.DataFrame(results)

# Adjust p-values (Benjamini-Hochberg)
from statsmodels.stats.multitest import multipletests
pvals = stats_df['p_value'].values
_, pvals_adj, _, _ = multipletests(pvals, method='fdr_bh')
stats_df['adjusted_p_value'] = pvals_adj

# View results
stats_df = stats_df.sort_values(by='adjusted_p_value')
stats_df




pops=['EC_glom', 'MC1', 'PEC','Podo']
mask=adata.obs['annotation_updated'].isin(pops)
glom=adata[mask].copy()
glom=glom[glom.obs['iPT_iLOH_ME_20um']!='Unknown'].copy()
glom

samples = glom.obs['unique_sample_identifier']
cluster_labels = glom.obs['iPT_iLOH_ME_20um']
data = pd.DataFrame({'iPT_iLOH_ME_20um': cluster_labels, 'unique_sample_identifier': samples})
grouped_data = data.groupby([ 'unique_sample_identifier', 'iPT_iLOH_ME_20um']).size().unstack(fill_value=0)
grouped_data

grouped_data_all_counts

grouped_data=grouped_data/grouped_data_all_counts
grouped_data=grouped_data.fillna(0)
grouped_data

import seaborn as sns
import matplotlib.pyplot as plt

# Reset index to bring 'unique_sample_identifier' into a column
grouped_data_reset = grouped_data.reset_index()

# Melt to long format
melted = grouped_data_reset.melt(id_vars='unique_sample_identifier',
                                 var_name='ME',
                                 value_name='Fraction')

# Make violin plot
plt.figure(figsize=(7, 4))
sns.violinplot(data=melted, x='ME', y='Fraction', inner='box', cut=0)
sns.stripplot(data=melted, x='ME', y='Fraction', color='k', size=2, alpha=0.5)

plt.xticks(rotation=0, ha='center')
plt.title('Fraction of Glomerular Cells per ME per Sample')
plt.tight_layout()
plt.savefig('Fraction of Glomerular Cells per ME per Sample.png', dpi=900, bbox_inches='tight')
plt.show()


from scipy.stats import mannwhitneyu

results = []
reference_group = 'Periglomerular ME'

# Loop over other MEs
for other_me in melted['ME'].unique():
    if other_me == reference_group:
        continue

    vals_profibrotic = melted[melted['ME'] == reference_group]['Fraction']
    vals_other = melted[melted['ME'] == other_me]['Fraction']

    stat, pval = mannwhitneyu(vals_profibrotic, vals_other, alternative='two-sided')

    results.append({
        'reference_ME': reference_group,
        'compared_ME': other_me,
        'U_statistic': stat,
        'p_value': pval
    })

# Create DataFrame
stats_df = pd.DataFrame(results)

# Adjust p-values (Benjamini-Hochberg)
from statsmodels.stats.multitest import multipletests
pvals = stats_df['p_value'].values
_, pvals_adj, _, _ = multipletests(pvals, method='fdr_bh')
stats_df['adjusted_p_value'] = pvals_adj

# View results
stats_df = stats_df.sort_values(by='adjusted_p_value')
stats_df




pops=['PT']
mask=adata.obs['annotation_updated'].isin(pops)
pt=adata[mask].copy()
pt=pt[pt.obs['iPT_iLOH_ME_20um']!='Unknown'].copy()
pt

samples = pt.obs['unique_sample_identifier']
cluster_labels = pt.obs['iPT_iLOH_ME_20um']
data = pd.DataFrame({'iPT_iLOH_ME_20um': cluster_labels, 'unique_sample_identifier': samples})
grouped_data = data.groupby([ 'unique_sample_identifier', 'iPT_iLOH_ME_20um']).size().unstack(fill_value=0)
grouped_data=grouped_data/grouped_data_all_counts
grouped_data=grouped_data.fillna(0)
grouped_data

import seaborn as sns
import matplotlib.pyplot as plt

# Reset index to bring 'unique_sample_identifier' into a column
grouped_data_reset = grouped_data.reset_index()

# Melt to long format
melted = grouped_data_reset.melt(id_vars='unique_sample_identifier',
                                 var_name='ME',
                                 value_name='Fraction')

# Make violin plot
plt.figure(figsize=(7, 4))
sns.violinplot(data=melted, x='ME', y='Fraction', inner='box', cut=0)
sns.stripplot(data=melted, x='ME', y='Fraction', color='k', size=2, alpha=0.5)

plt.xticks(rotation=0, ha='center')
plt.title('Fraction of PT Cells per ME per Sample')
plt.tight_layout()
plt.savefig('Fraction of PT Cells per ME per Sample.png', dpi=900, bbox_inches='tight')
plt.show()


from scipy.stats import mannwhitneyu

results = []
reference_group = 'Proximal ME'

# Loop over other MEs
for other_me in melted['ME'].unique():
    if other_me == reference_group:
        continue

    vals_profibrotic = melted[melted['ME'] == reference_group]['Fraction']
    vals_other = melted[melted['ME'] == other_me]['Fraction']

    stat, pval = mannwhitneyu(vals_profibrotic, vals_other, alternative='two-sided')

    results.append({
        'reference_ME': reference_group,
        'compared_ME': other_me,
        'U_statistic': stat,
        'p_value': pval
    })

# Create DataFrame
stats_df = pd.DataFrame(results)

# Adjust p-values (Benjamini-Hochberg)
from statsmodels.stats.multitest import multipletests
pvals = stats_df['p_value'].values
_, pvals_adj, _, _ = multipletests(pvals, method='fdr_bh')
stats_df['adjusted_p_value'] = pvals_adj

# View results
stats_df = stats_df.sort_values(by='adjusted_p_value')
stats_df






pops=['TAL']
mask=adata.obs['annotation_updated'].isin(pops)
tal=adata[mask].copy()
tal=tal[tal.obs['iPT_iLOH_ME_20um']!='Unknown'].copy()
tal

samples = tal.obs['unique_sample_identifier']
cluster_labels = tal.obs['iPT_iLOH_ME_20um']
data = pd.DataFrame({'iPT_iLOH_ME_20um': cluster_labels, 'unique_sample_identifier': samples})
grouped_data = data.groupby([ 'unique_sample_identifier', 'iPT_iLOH_ME_20um']).size().unstack(fill_value=0)
grouped_data=grouped_data/grouped_data_all_counts
grouped_data=grouped_data.fillna(0)
grouped_data

import seaborn as sns
import matplotlib.pyplot as plt

# Reset index to bring 'unique_sample_identifier' into a column
grouped_data_reset = grouped_data.reset_index()

# Melt to long format
melted = grouped_data_reset.melt(id_vars='unique_sample_identifier',
                                 var_name='ME',
                                 value_name='Fraction')

# Make violin plot
plt.figure(figsize=(7, 4))
sns.violinplot(data=melted, x='ME', y='Fraction', inner='box', cut=0)
sns.stripplot(data=melted, x='ME', y='Fraction', color='k', size=2, alpha=0.5)

plt.xticks(rotation=0, ha='center')
plt.title('Fraction of TAL Cells per ME per Sample')
plt.tight_layout()
plt.savefig('Fraction of TAL Cells per ME per Sample.png', dpi=900, bbox_inches='tight')
plt.show()


from scipy.stats import mannwhitneyu

results = []
reference_group = 'Distal ME'

# Loop over other MEs
for other_me in melted['ME'].unique():
    if other_me == reference_group:
        continue

    vals_profibrotic = melted[melted['ME'] == reference_group]['Fraction']
    vals_other = melted[melted['ME'] == other_me]['Fraction']

    stat, pval = mannwhitneyu(vals_profibrotic, vals_other, alternative='two-sided')

    results.append({
        'reference_ME': reference_group,
        'compared_ME': other_me,
        'U_statistic': stat,
        'p_value': pval
    })

# Create DataFrame
stats_df = pd.DataFrame(results)

# Adjust p-values (Benjamini-Hochberg)
from statsmodels.stats.multitest import multipletests
pvals = stats_df['p_value'].values
_, pvals_adj, _, _ = multipletests(pvals, method='fdr_bh')
stats_df['adjusted_p_value'] = pvals_adj

# View results
stats_df = stats_df.sort_values(by='adjusted_p_value')
stats_df




pops=['EC_Peritub', 'VSMC', 'EC_DVR']
mask=adata.obs['annotation_updated'].isin(pops)
ec=adata[mask].copy()
ec=ec[ec.obs['iPT_iLOH_ME_20um']!='Unknown'].copy()
ec

ec.obs['annotation_updated'].value_counts()

samples = ec.obs['unique_sample_identifier']
cluster_labels = ec.obs['iPT_iLOH_ME_20um']
data = pd.DataFrame({'iPT_iLOH_ME_20um': cluster_labels, 'unique_sample_identifier': samples})
grouped_data = data.groupby([ 'unique_sample_identifier', 'iPT_iLOH_ME_20um']).size().unstack(fill_value=0)
grouped_data=grouped_data/grouped_data_all_counts
grouped_data=grouped_data.fillna(0)
grouped_data

import seaborn as sns
import matplotlib.pyplot as plt

# Reset index to bring 'unique_sample_identifier' into a column
grouped_data_reset = grouped_data.reset_index()

# Melt to long format
melted = grouped_data_reset.melt(id_vars='unique_sample_identifier',
                                 var_name='ME',
                                 value_name='Fraction')

# Make violin plot
plt.figure(figsize=(7, 4))
sns.violinplot(data=melted, x='ME', y='Fraction', inner='box', cut=0)
sns.stripplot(data=melted, x='ME', y='Fraction', color='k', size=2, alpha=0.5)

plt.xticks(rotation=0, ha='center')
plt.title('Fraction of EC Cells per ME per Sample')
plt.tight_layout()
plt.savefig('Fraction of EC Cells per ME per Sample.png', dpi=900, bbox_inches='tight')
plt.show()


from scipy.stats import mannwhitneyu

results = []
reference_group = 'Perivascular ME'

# Loop over other MEs
for other_me in melted['ME'].unique():
    if other_me == reference_group:
        continue

    vals_profibrotic = melted[melted['ME'] == reference_group]['Fraction']
    vals_other = melted[melted['ME'] == other_me]['Fraction']

    stat, pval = mannwhitneyu(vals_profibrotic, vals_other, alternative='two-sided')

    results.append({
        'reference_ME': reference_group,
        'compared_ME': other_me,
        'U_statistic': stat,
        'p_value': pval
    })

# Create DataFrame
stats_df = pd.DataFrame(results)

# Adjust p-values (Benjamini-Hochberg)
from statsmodels.stats.multitest import multipletests
pvals = stats_df['p_value'].values
_, pvals_adj, _, _ = multipletests(pvals, method='fdr_bh')
stats_df['adjusted_p_value'] = pvals_adj

# View results
stats_df = stats_df.sort_values(by='adjusted_p_value')
stats_df



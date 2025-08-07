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

###first making pie charts for center cell and 20um environment:

adata=sc.read_h5ad('/content/drive/MyDrive/Bernhard/revision_nature/revision_object.h5ad')
adata_subset=adata[adata.obs['immune_ME']!='Unknown'].copy()
colors_immune={
    'Immune': '#2B2B2B',
    'CD4+': '#60BEB2',
    'Baso_Mast': '#2B2B2B',
    'B': '#ABA6CE',
    'CD8+':'#9A1A2D',
    'Neutrophil': '#2B2B2B',
    'Macro': '#FFD600',
    'NK': '#2B2B2B',
    'cDC': '#2B2B2B',
    'mDC': '#2B2B2B',
    'pDC':'#2B2B2B',
    'Plasma': '#F07E22',
    'other': '#2B2B2B',
}


# Assuming new_adata_subset is your AnnData object
# Extract the relevant columns
samples = adata_subset.obs['immune_ME']
cluster_labels = adata_subset.obs['immune_cell_annotation_combined']

# Create a DataFrame from the relevant columns
data = pd.DataFrame({'immune_ME': samples, 'immune_cell_annotation_combined': cluster_labels})
# Group by 'sample' and 'cluster_labels_annotated' and count the occurrences
grouped_data = data.groupby(['immune_ME', 'immune_cell_annotation_combined']).size().unstack(fill_value=0)
grouped_data['other']=grouped_data['Baso_Mast']+grouped_data['cDC']+grouped_data['mDC']+grouped_data['pDC']+grouped_data['NK']+grouped_data['Neutrophil']
grouped_data=grouped_data.drop(['Baso_Mast', 'cDC', 'mDC', 'pDC', 'NK', 'Neutrophil'], axis=1)

colors = [colors_immune.get(col, 'grey') for col in grouped_data.columns]
# Plot pie charts for each subcluster
for subcluster in grouped_data.index:
    subcluster_data = grouped_data.loc[subcluster]
    print(subcluster)
    # Plot pie chart
    plt.figure(figsize=(6, 6))
    plt.pie(subcluster_data, startangle=0,colors=colors)
    plt.title(f'')
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig(f'Pie Chart {subcluster} inner circle.png', dpi=900, bbox_inches='tight',transparent=True)
    plt.show()


grouped_data

adata

adata_subset=adata[adata.obs['Immune_ME_20um']!='Unknown'].copy()

colors_neighbors={
    'Immune': '#EF785F',
    'Fibroblast': '#D279AF',
    'PT_TAL': '#8D574C',
    'iPT_iTAL': '#88BD2C',
    'glom':'#9068AA',
    'ec': '#1B79B6',
    'other': '#2B2B2B',
}


# Extract the relevant columns
samples = adata_subset.obs['Immune_ME_20um']
cluster_labels = adata_subset.obs['annotation_updated']

# Create a DataFrame from the relevant columns
data = pd.DataFrame({'Immune_ME_20um': samples, 'annotation_updated': cluster_labels})
# Group by 'sample' and 'cluster_labels_annotated' and count the occurrences
grouped_data = data.groupby(['Immune_ME_20um', 'annotation_updated']).size().unstack(fill_value=0)
grouped_data
grouped_data['PT_TAL']=grouped_data['PT']+grouped_data['TAL']
grouped_data=grouped_data.drop(['PT', 'TAL'], axis=1)
grouped_data['iPT_iTAL']=grouped_data['iPT']+grouped_data['iTAL']
grouped_data=grouped_data.drop(['iPT', 'iTAL'], axis=1)
grouped_data['glom']=grouped_data['EC_glom']+grouped_data['MC1']+grouped_data['PEC']+grouped_data['Podo']
grouped_data=grouped_data.drop(['EC_glom', 'MC1', 'PEC', 'Podo'], axis=1)
grouped_data['ec']=grouped_data['EC_DVR']+grouped_data['EC_Peritub']+grouped_data['VSMC']
grouped_data=grouped_data.drop(['EC_DVR', 'EC_Peritub', 'VSMC'], axis=1)
grouped_data['other']=grouped_data['CNT']+grouped_data['DCT']+grouped_data['DTL_ATL']+grouped_data['EC_Lymph']+grouped_data['IC A']+grouped_data['IC B']+grouped_data['PC']
grouped_data=grouped_data.drop(['CNT', 'DCT', 'DTL_ATL', 'IC A', 'IC B', 'PC', 'EC_Lymph'], axis=1)

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
    plt.savefig(f'Pie Chart {subcluster} entire neighborhood.png', dpi=900, bbox_inches='tight',transparent=True)
    plt.show()


###second making diagrams for cell fractions and statistical comparison:
###first tubular immune ME
excluded_conditions = [
    'DM', 'DM/HTN', 'CKD', 'C3GN', 'AA amyloid', 'IgA',
    'FSGS', 'MN', 'HTN', 'DKD+FSGS', 'TMA'
]

# Build annotation DataFrame
df_annotations = adata.obs[['Condition', 'unique_sample_identifier', 'annotation_updated', 'Immune_ME_20um', 'immune_cell_annotation_combined']].copy()
df_annotations.columns = ['type', 'sample', 'annotation', 'niche', 'cluster_labels']

# Filter for valid disease types
df_annotations = df_annotations[~df_annotations['type'].isin(excluded_conditions)]
df_annotations['type'] = df_annotations['type'].cat.remove_unused_categories()
df_annotations['sample'] = df_annotations['sample'].cat.remove_unused_categories()

# Filter for niche of interest
niches_to_keep = ['Inj. Tubular Immune ME']
df_filtered = df_annotations[df_annotations['niche'].isin(niches_to_keep)].copy()
df_filtered['niche'] = df_filtered['niche'].cat.remove_unused_categories()
df_filtered['sample'] = df_filtered['sample'].cat.remove_unused_categories()
df_filtered = df_filtered.drop(columns='niche')

# Total cell counts per sample
total_counts = (
    df_filtered
    .groupby('sample', as_index=False)
    .size()
    .rename(columns={'size': 'total_cell_count'})
)

# Subset to immune cells only
df_immune = df_filtered[df_filtered['annotation'] == 'Immune']

# Count immune cluster composition per sample
cluster_counts = df_immune.groupby(['sample', 'cluster_labels']).size().unstack(fill_value=0)

# Merge total cell count and cluster counts
merged = pd.merge(cluster_counts, total_counts, on='sample')

df_type = df_annotations[['sample', 'type']].drop_duplicates().set_index('sample')
merged = merged.join(df_type, on='sample')

# Drop 'Unknown' cluster label if present
merged = merged.drop(columns='Unknown', errors='ignore')

celltypes = ['B', 'CD4+', 'CD8+', 'Macro', 'NK', 'Neutrophil', 'Plasma']

# Compute fraction per cell type
for cell in celltypes:
    merged[f'{cell}_fraction'] = merged.get(cell, 0) / merged['total_cell_count']

# Melt for visualization
df_melted = merged.melt(
    id_vars='type',
    value_vars=[f'{c}_fraction' for c in celltypes],
    var_name='celltype',
    value_name='fraction'
)

# Compute control medians
control_medians = (
    df_melted[df_melted['type'] == 'Control']
    .groupby('celltype')['fraction']
    .median()
)

# Center all values on Control median
df_melted['fraction_centered'] = df_melted.apply(
    lambda row: row['fraction'] - control_medians[row['celltype']], axis=1
)

from scipy.stats import ttest_ind

stats = []
for cell in df_melted['celltype'].unique():
    subset = df_melted[df_melted['celltype'] == cell]
    dkd_vals = subset[subset['type'] == 'DKD']['fraction']
    ctrl_vals = subset[subset['type'] == 'Control']['fraction']

    stat, pval = ttest_ind(dkd_vals, ctrl_vals, equal_var=False)
    stats.append({
        'celltype': cell,
        'DKD_mean': dkd_vals.mean(),
        'Control_mean': ctrl_vals.mean(),
        't_statistic': stat,
        'p_value': pval
    })

stats_df = pd.DataFrame(stats).sort_values(by='p_value')


plt.figure(figsize=(7, 5))
ax = sns.boxplot(
    data=df_melted,
    x='celltype',
    y='fraction_centered',
    hue='type',
    showcaps=False,
    whiskerprops={'linewidth': 0},
    fliersize=0,
    linewidth=1,
    linecolor='black',
    palette='Set2',
)

# Add horizontal reference line
#plt.axhline(0, linestyle='--', color='black')

# Set axis and spine styles
plt.xticks(rotation=90)
plt.ylabel('Centered Celltype Fraction\n(Relative to Control)')
plt.ylim(-0.019, 0.031)
plt.xlabel('')
plt.title('Relative Immune Cell Enrichment in Tubular Microenvironment')

# Remove top, right, and bottom spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

# Remove ticks from both axes
ax.tick_params(axis='x', bottom=False)
ax.tick_params(axis='y', left=False)

# Move legend outside
plt.legend(title='Type', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

plt.tight_layout()
plt.savefig(f'boxplot_fraction_pval_tubular_immune_me.png', dpi=900, transparent=True)
plt.show()


#now glom
sample_metadata = pd.read_csv('/content/diagnosis.csv', sep=';').set_index('Sample ID')
sample_metadata = sample_metadata.drop(columns=['Disease', 'Race', 'Sex', 'DM', 'HTN'])

# Exclude undesired conditions
excluded_conditions = ['DM', 'DM/HTN', 'CKD', 'C3GN', 'AA amyloid', 'IgA', 'FSGS', 'MN', 'HTN', 'DKD+FSGS', 'TMA']
sample_metadata_filtered = sample_metadata[~sample_metadata['Condition'].isin(excluded_conditions)]
sample_metadata_filtered = sample_metadata_filtered.drop(columns=['Condition', 'Age'])

df_annotations = adata.obs[['Condition', 'unique_sample_identifier', 'annotation_updated', 'niches_annotation_based', 'immune_cell_annotation_combined']].copy()
df_annotations.columns = ['type', 'sample', 'annotation', 'niche', 'cluster_labels']

# Remove excluded conditions
df_annotations = df_annotations[~df_annotations['type'].isin(excluded_conditions)]
df_annotations['type'] = df_annotations['type'].cat.remove_unused_categories()
df_annotations['sample'] = df_annotations['sample'].cat.remove_unused_categories()

# Focus on Glomerular niche
df_glom = df_annotations[df_annotations['niche'] == 'Glomerular niche'].copy()
df_glom['niche'] = df_glom['niche'].cat.remove_unused_categories()
df_glom['sample'] = df_glom['sample'].cat.remove_unused_categories()
df_glom = df_glom.drop(columns='niche')

# Total glomerular cells per sample
total_counts = df_glom.groupby('sample').size().reset_index(name='total_cell_count')

# Filter to immune cells
df_glom_immune = df_glom[df_glom['annotation'] == 'Immune']

# Immune cluster counts
cluster_counts = df_glom_immune.groupby(['sample', 'cluster_labels']).size().unstack(fill_value=0)

# Merge with total counts
df_merged = pd.merge(cluster_counts, total_counts, on='sample')
df_merged = df_merged[df_merged['total_cell_count'] > 50]

df_types = df_annotations[['sample', 'type']].drop_duplicates().set_index('sample')
df_merged = df_merged.join(df_types, on='sample')

# Drop unknown cluster label if present
df_merged = df_merged.drop(columns='Unknown', errors='ignore')

celltypes = ['B', 'CD4+', 'CD8+', 'Macro', 'NK', 'Neutrophil', 'Plasma']
for ctype in celltypes:
    df_merged[f'{ctype}_fraction'] = df_merged.get(ctype, 0) / df_merged['total_cell_count']

df_melted = df_merged.melt(
    id_vars='type',
    value_vars=[f'{ctype}_fraction' for ctype in celltypes],
    var_name='celltype',
    value_name='fraction'
)

control_medians = df_melted[df_melted['type'] == 'Control'].groupby('celltype')['fraction'].median()
df_melted['fraction_centered'] = df_melted.apply(
    lambda row: row['fraction'] - control_medians.get(row['celltype'], 0), axis=1
)

from scipy.stats import ttest_ind

df_stats = []
for celltype in df_melted['celltype'].unique():
    group = df_melted[df_melted['celltype'] == celltype]
    dkd = group[group['type'] == 'DKD']['fraction']
    control = group[group['type'] == 'Control']['fraction']

    stat, pval = ttest_ind(dkd, control, equal_var=False)
    df_stats.append({
        'celltype': celltype,
        'DKD_mean': dkd.mean(),
        'Control_mean': control.mean(),
        't_statistic': stat,
        'p_value': pval
    })

stats_df = pd.DataFrame(df_stats).sort_values('p_value')

plt.figure(figsize=(7, 5))
ax = sns.boxplot(
    data=df_melted,
    x='celltype',
    y='fraction_centered',
    hue='type',
    showcaps=False,
    whiskerprops={'linewidth': 0},
    fliersize=0,
    linewidth=1,
    linecolor='black',
    palette='Set2',
)

# Add horizontal reference line
#plt.axhline(0, linestyle='--', color='black')

# Set axis and spine styles
plt.xticks(rotation=90)
plt.ylabel('Centered Celltype Fraction\n(Relative to Control)')
plt.ylim(-0.003, 0.02)
plt.xlabel('')
plt.title('Relative Immune Cell Enrichment in Glomerular Microenvironment')

# Remove top, right, and bottom spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

# Remove ticks from both axes
ax.tick_params(axis='x', bottom=False)
ax.tick_params(axis='y', left=False)

# Move legend outside
plt.legend(title='Type', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

plt.tight_layout()
plt.savefig(f'boxplot_fraction_pval_glomerular_immune_me.png', dpi=900, transparent=True)
plt.show()

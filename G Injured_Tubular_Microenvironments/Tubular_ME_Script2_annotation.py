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
import leidenalg as la
from pathlib import Path
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap
import seaborn as sns
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import MinMaxScaler
from scipy.cluster.hierarchy import dendrogram, linkage

sc.settings.figdir = "/home/bcd/revision_nature/iPT_iLOH_ME/plots_script2"  # Update this path as needed
sc.settings.set_figure_params(dpi=900)

new_adata=sc.read_h5ad('/home/bcd/revision_nature/iPT_iLOH_ME/xenium_cosmx_neighbors_as_genes_iTAL_iPT_subset.h5ad')
print(new_adata.obs['unique_sample_identifier'].unique())
for identifier in new_adata.obs['unique_sample_identifier'].unique():
    print(identifier)
new_adata.obs["cluster_labels"].value_counts()
cell_identities = {0: 'Profibrotic ME', 1: 'Distal ME', 2: 'Proximal ME', 3: 'Perivascular ME',
                   4: 'Periglomerular ME', 5: 'Distal ME'}
new_adata.obs["cluster_labels_annotated_iPT_iLOH"] = new_adata.obs['cluster_labels'].map(cell_identities).astype('category')
sc.pl.umap(new_adata, color='cluster_labels', save='_cluster_labels.png')
sc.pl.umap(new_adata, color='cluster_labels_annotated_iPT_iLOH', save='_cluster_labels_annotated_iPT_iLOH.png')

for cluster in new_adata.obs['cluster_labels_annotated_iPT_iLOH'].unique():
    subset=new_adata[new_adata.obs['cluster_labels_annotated_iPT_iLOH']==cluster].copy()
    print(cluster)
    print(subset.obs['tech'].value_counts())
    sc.pl.umap(new_adata, color = "cluster_labels_annotated_iPT_iLOH", groups=cluster, save=f'_cluster_labels_annotated_iPT_iLOH_{cluster}.png')
    sc.pl.umap(subset, color = "annotation_updated", save=f'_annotation_updated_{cluster}.png')

# Extract the relevant columns
annotations = new_adata.obs["annotation_updated"]
cluster_labels = new_adata.obs["cluster_labels_annotated_iPT_iLOH"]

# Create a contingency table (cross-tabulation)
contingency_table = pd.crosstab(cluster_labels, annotations)

# Normalize the data across the columns (annotations)
normalized_table = contingency_table.div(contingency_table.sum(axis=0), axis=1)

# Scale the data from 0 to 1 across the columns
scaler = MinMaxScaler()
scaled_table = pd.DataFrame(scaler.fit_transform(normalized_table),
                            index=normalized_table.index,
                            columns=normalized_table.columns)

desired_column_order = ['PT','EC_Peritub','iPT','TAL','DTL_ATL','iTAL','CNT', 'IC A', 'IC B', 'PC','DCT', 'EC_DVR', 'VSMC',
                  'Fibroblast',  'Immune','EC_Lymph',   'EC_glom','PEC',
                    'Podo',  'MC1']  # Replace with actual annotation labels

desired_index_order = ['Profibrotic ME','Distal ME','Proximal ME', 'Perivascular ME',
                   'Periglomerular ME']  # Replace with actual cluster labels
scaled_table = scaled_table.reindex(index=desired_index_order, columns=desired_column_order)

# Plot the heatmap
plt.figure(figsize=(7, 10))
sns.heatmap(scaled_table.T, annot=False, fmt=".2f", cmap="Blues")
plt.title('iPT iLOH Microenvironments')
plt.tick_params(axis='x', which='both', length=0)
plt.tick_params(axis='y', which='both', length=0)
plt.xlabel('Niche Labels')
plt.ylabel('Cell Types')
plt.savefig('/home/bcd/revision_nature/iPT_iLOH_ME/plots_script2/injured_ME.png', dpi=900, bbox_inches='tight')
plt.close()

samples = new_adata.obs['unique_sample_identifier']
cluster_labels = new_adata.obs['cluster_labels_annotated_iPT_iLOH']

# Create a DataFrame from the relevant columns
data = pd.DataFrame({'unique_sample_identifier': samples, 'cluster_labels_annotated_iPT_iLOH': cluster_labels})

# Group by 'sample' and 'cluster_labels_annotated' and count the occurrences
grouped_data = data.groupby(['unique_sample_identifier', 'cluster_labels_annotated_iPT_iLOH']).size().unstack(fill_value=0)

# Get sample totals in the same order as grouped_data's index
sample_totals = data['unique_sample_identifier'].value_counts()
sample_totals = sample_totals.loc[grouped_data.index]  # Match the index order

# Normalize the counts
normalized_count_df = grouped_data.div(sample_totals, axis='index')
ax = normalized_count_df.plot(kind='bar', stacked=True, figsize=(32, 8), colormap='tab20')
ax.grid(False)
# Set plot labels and title
plt.xlabel('')
plt.ylabel('Relative Amount')
plt.title('Sample Niche Composition')
plt.legend(title='Niche Labels', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
plt.tick_params(axis='x', which='both', length=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('/home/bcd/revision_nature/iPT_iLOH_ME/plots_script2/iPT_iLOH_ME_Samples.png', dpi=900, bbox_inches='tight')
plt.close()


new_adata.write('/home/bcd/revision_nature/iPT_iLOH_ME/xenium_cosmx_neighbors_as_genes_iTAL_iPT_subset.h5ad')
adata = sc.read_h5ad('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad')

adata.obs['iPT_iLOH_ME']=new_adata.obs['cluster_labels_annotated_iPT_iLOH'].copy()
adata.obs['iPT_iLOH_ME']=adata.obs['iPT_iLOH_ME'].astype('category')
adata.obs['iPT_iLOH_ME'] = adata.obs['iPT_iLOH_ME'].cat.add_categories('Unknown')
adata.obs['iPT_iLOH_ME']=adata.obs['iPT_iLOH_ME'].fillna('Unknown')
sc.pl.umap(adata, color = "iPT_iLOH_ME", save='_iPT_iLOH_ME.png', palette='tab10')
adata.write('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad')
print(adata)

samples = adata.obs['unique_sample_identifier']
cluster_labels = adata.obs['iPT_iLOH_ME']

# Create a DataFrame from the relevant columns
data = pd.DataFrame({'unique_sample_identifier': samples, 'iPT_iLOH_ME': cluster_labels})

# Group by 'sample' and 'niches_annotation_based' and count the occurrences
grouped_data = data.groupby(['unique_sample_identifier', 'iPT_iLOH_ME']).size().unstack(fill_value=0)

# Get sample totals in the same order as grouped_data's index
sample_totals = data['unique_sample_identifier'].value_counts()
sample_totals = sample_totals.loc[grouped_data.index]  # Match the index order

# Normalize the counts
normalized_count_df = grouped_data.div(sample_totals, axis='index')
ax = normalized_count_df.plot(kind='bar', stacked=True, figsize=(32, 8), colormap='tab20')
ax.grid(False)
# Set plot labels and title
plt.xlabel('')
plt.ylabel('Relative Amount')
plt.title('Sample Niche Composition')
plt.legend(title='Niche Labels', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
plt.tick_params(axis='x', which='both', length=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('/home/bcd/revision_nature/iPT_iLOH_ME/plots_script2/iPT_iLOH_with_unknown.png', dpi=900, bbox_inches='tight')
plt.close()

# CosMx Analysis
adata_cosmx = adata[adata.obs['tech'] == 'CosMx'].copy()
adata_cosmx_subset = adata_cosmx[adata_cosmx.obs['iPT_iLOH_ME'] != 'Unknown'].copy()
adata_cosmx_subset.X = adata_cosmx_subset.layers['counts'].copy()
sc.pp.normalize_total(adata_cosmx_subset, target_sum=1e4)
sc.pp.log1p(adata_cosmx_subset)
sc.tl.rank_genes_groups(adata_cosmx_subset, groupby='iPT_iLOH_ME', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups(adata_cosmx_subset, n_genes=25, sharey=False, save='_iPT_iLOH_ME_cosmx_all.png')

# CosMx iPT
adata_cosmx_ipt = adata_cosmx_subset[adata_cosmx_subset.obs['annotation_postscanvi_level2'] == 'iPT'].copy()
adata_cosmx_ipt.X = adata_cosmx_ipt.layers['counts'].copy()
sc.pp.normalize_total(adata_cosmx_ipt, target_sum=1e4)
sc.pp.log1p(adata_cosmx_ipt)
sc.tl.rank_genes_groups(adata_cosmx_ipt, groupby='iPT_iLOH_ME', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups(adata_cosmx_ipt, n_genes=25, sharey=False, save='_iPT_iLOH_ME_cosmx_iPT.png')

# CosMx iTAL
adata_cosmx_ital = adata_cosmx_subset[adata_cosmx_subset.obs['annotation_postscanvi_level2'] == 'iTAL'].copy()
adata_cosmx_ital.X = adata_cosmx_ital.layers['counts'].copy()
sc.pp.normalize_total(adata_cosmx_ital, target_sum=1e4)
sc.pp.log1p(adata_cosmx_ital)
sc.tl.rank_genes_groups(adata_cosmx_ital, groupby='iPT_iLOH_ME', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups(adata_cosmx_ital, n_genes=25, sharey=False, save='_iPT_iLOH_ME_cosmx_iTAL.png')

# Xenium Analysis
adata_xenium = adata[adata.obs['tech'] == 'Xenium'].copy()
adata_xenium_subset = adata_xenium[adata_xenium.obs['iPT_iLOH_ME'] != 'Unknown'].copy()
adata_xenium_subset.X = adata_xenium_subset.layers['counts'].copy()
sc.pp.normalize_total(adata_xenium_subset, target_sum=1e4)
sc.pp.log1p(adata_xenium_subset)
sc.tl.rank_genes_groups(adata_xenium_subset, groupby='iPT_iLOH_ME', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups(adata_xenium_subset, n_genes=25, sharey=False, save='_iPT_iLOH_ME_xenium_all.png')

# Xenium iPT
adata_xenium_ipt = adata_xenium_subset[adata_xenium_subset.obs['annotation_postscanvi_level2'] == 'iPT'].copy()
adata_xenium_ipt.X = adata_xenium_ipt.layers['counts'].copy()
sc.pp.normalize_total(adata_xenium_ipt, target_sum=1e4)
sc.pp.log1p(adata_xenium_ipt)
sc.tl.rank_genes_groups(adata_xenium_ipt, groupby='iPT_iLOH_ME', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups(adata_xenium_ipt, n_genes=25, sharey=False, save='_iPT_iLOH_ME_xenium_iPT.png')

# Xenium iTAL
adata_xenium_ital = adata_xenium_subset[adata_xenium_subset.obs['annotation_postscanvi_level2'] == 'iTAL'].copy()
adata_xenium_ital.X = adata_xenium_ital.layers['counts'].copy()
sc.pp.normalize_total(adata_xenium_ital, target_sum=1e4)
sc.pp.log1p(adata_xenium_ital)
sc.tl.rank_genes_groups(adata_xenium_ital, groupby='iPT_iLOH_ME', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups(adata_xenium_ital, n_genes=25, sharey=False, save='_iPT_iLOH_ME_xenium_iTAL.png')

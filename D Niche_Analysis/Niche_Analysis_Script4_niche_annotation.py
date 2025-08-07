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

sc.settings.figdir = "/home/bcd/revision_nature/neighbors_neighborhood/plots_script5"  # Update this path as needed
sc.settings.set_figure_params(dpi=900)

new_adata=sc.read_h5ad('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_neighbors_as_genes.h5ad')
print(new_adata.obs['unique_sample_identifier'].unique())
for identifier in new_adata.obs['unique_sample_identifier'].unique():
    print(identifier)
new_adata.obs["cluster_labels"].value_counts()
cell_identities = {0: 'Glomerular niche', 1: 'PC niche', 2: 'DCT niche', 3: 'Fibroblast niche',
                   4: 'PT niche', 5: 'LOH niche', 6: 'Glomerular niche', 7: 'Immune niche',
                   8: 'Glomerular niche', 9: 'iPT niche', 10:'Vascular niche', 11:'Vascular niche',
                   12:'Fibroblast niche', 13:'LOH niche', 14:'iLOH niche', 15:'CNT niche',16:'PT niche'}
new_adata.obs["cluster_labels_annotated"] = new_adata.obs['cluster_labels'].map(cell_identities).astype('category')
sc.pl.umap(new_adata, color='cluster_labels', save='_cluster_lables.png')
sc.pl.umap(new_adata, color='cluster_labels_annotated', save='_cluster_labels_annotated.png')

for cluster in new_adata.obs['cluster_labels_annotated'].unique():
    subset=new_adata[new_adata.obs['cluster_labels_annotated']==cluster].copy()
    print(cluster)
    print(subset.obs['tech'].value_counts())
    sc.pl.umap(new_adata, color = "cluster_labels_annotated", groups=cluster, save=f'_cluster_labels_annotated_{cluster}.png')
    sc.pl.umap(subset, color = "annotation_updated", save=f'_annotation_updated_{cluster}.png')

# Extract the relevant columns
annotations = new_adata.obs["annotation_updated"]
cluster_labels = new_adata.obs["cluster_labels_annotated"]

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

desired_index_order = ['PT niche', 'iPT niche','LOH niche', 'iLOH niche',
                       'CNT niche','PC niche','DCT niche', 'Vascular niche',
                   'Fibroblast niche','Immune niche','Glomerular niche']  # Replace with actual cluster labels
scaled_table = scaled_table.reindex(index=desired_index_order, columns=desired_column_order)
# Plot the heatmap
plt.figure(figsize=(7, 10))
sns.heatmap(scaled_table.T, annot=False, fmt=".2f", cmap="Blues")
plt.title('Kidney Niches')
plt.tick_params(axis='x', which='both', length=0)
plt.tick_params(axis='y', which='both', length=0)
plt.xlabel('Niche Labels')
plt.ylabel('Cell Types')
plt.savefig('/home/bcd/revision_nature/neighbors_neighborhood/plots_script5/Kindey_Niches_first.png', dpi=900, bbox_inches='tight')
plt.close()

samples = new_adata.obs['unique_sample_identifier']
cluster_labels = new_adata.obs['cluster_labels_annotated']

# Create a DataFrame from the relevant columns
data = pd.DataFrame({'unique_sample_identifier': samples, 'cluster_labels_annotated': cluster_labels})

# Group by 'sample' and 'cluster_labels_annotated' and count the occurrences
grouped_data = data.groupby(['unique_sample_identifier', 'cluster_labels_annotated']).size().unstack(fill_value=0)

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
plt.savefig('/home/bcd/revision_nature/neighbors_neighborhood/plots_script5/Kindey_Niches_Samples.png', dpi=900, bbox_inches='tight')
plt.close()

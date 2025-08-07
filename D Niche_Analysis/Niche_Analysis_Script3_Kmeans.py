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

sc.settings.figdir = "/home/bcd/revision_nature/neighbors_neighborhood/plots_script4"  # Update this path as needed
sc.settings.set_figure_params(dpi=900)
adata = sc.read_h5ad('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad')
print(adata)
print(adata.obs['unique_sample_identifier'].value_counts())
print(adata.obs['annotation_updated'].value_counts())

df_concat = pd.read_csv('/home/bcd/revision_nature/neighbors_neighborhood/neighborhood_dataframe.csv', index_col=0)
df_concat.index.name = None 
print(df_concat)

import anndata
# Convert merged_df to a sparse matrix
sparse_matrix = scipy.sparse.csr_matrix(df_concat.values)
# Create a new AnnData object with the sparse matrix
new_adata = anndata.AnnData(X=sparse_matrix)
# Assign merged_df columns to .var
new_adata.var = pd.DataFrame(index=df_concat.columns)
new_adata.var['feature_names'] = df_concat.columns
# Assign merged_df index to .obs
new_adata.obs = pd.DataFrame(index=df_concat.index)
# Verify the new AnnData object
print(new_adata)

# Subset adata to only include cells present in new_adata
adata_subset = adata[adata.obs_names.isin(new_adata.obs_names)].copy()
# Reorder adata_subset to match new_adata
adata_subset = adata_subset[new_adata.obs_names]

assert all(adata_subset.obs_names == new_adata.obs_names), "Cell ordering mismatch!"
# Now safely transfer metadata
obs_to_map = [
    'n_neighbors_20um', 'n_neighbors_40um', 'n_neighbors_60um', 'n_neighbors_80um', 'n_neighbors_200um',
    'immune_cell_annotation_xenium', 'b_cell_annotation_xenium', 'unique_sample_identifier', 'annotation_updated', 'tech',
]
for obs in obs_to_map:
    new_adata.obs[obs] = adata_subset.obs[obs].copy()

# Copy UMAP coordinates
new_adata.obsm['X_umap'] = adata_subset.obsm['X_umap'].copy()
# Save AnnData
new_adata.write('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_neighbors_as_genes.h5ad')

# Normalize the data
scaler = StandardScaler()
features_normalized = scaler.fit_transform(df_concat)

# Run PCA again (if not already)
pca = PCA()
pca.fit(features_normalized)

# Plot explained variance
plt.figure(figsize=(8, 5))
plt.plot(np.cumsum(pca.explained_variance_ratio_), marker='o')
plt.xlabel('Number of Principal Components')
plt.ylabel('Cumulative Explained Variance')
plt.title('Explained Variance by PCA Components')
plt.grid(True)
plt.tight_layout()
plt.axhline(y=0.9, color='r', linestyle='--', label='90% variance')
plt.legend()
plt.savefig('/home/bcd/revision_nature/neighbors_neighborhood/plots_script4/pca_variance_curve.png', dpi=900)
plt.close()

pca = PCA(n_components=25)  
principal_components = pca.fit_transform(features_normalized)

inertias = []
k = range(1, 20)
# dissimilarity would not be defined for a single cluster, thus, minimum number of clusters should be 2
for K in k:
    # Create a MiniBatchKMeans object, arbitrarily choosing 10 clusters
    kmeans = MiniBatchKMeans(n_clusters=K, random_state=42)
    # Fit the model to the data
    kmeans.fit(principal_components)
    inertias.append(kmeans.inertia_)

plt.plot(k, inertias, 'bx-')
plt.xlabel('Values of K')
plt.ylabel('Inertia')
plt.title('The Elbow Method using Inertia')
plt.tight_layout()
plt.savefig('/home/bcd/revision_nature/neighbors_neighborhood/plots_script4/elbow_kmeans.png', dpi=900)
plt.close()

# Perform MiniBatchKMeans clustering
kmeans = MiniBatchKMeans(n_clusters=17, random_state=42)
labels = kmeans.fit_predict(principal_components)
df_concat['cluster_labels'] = labels

new_adata.obs['cluster_labels']=df_concat['cluster_labels'].copy()
new_adata.obs['cluster_labels']=new_adata.obs['cluster_labels'].astype('category')
sc.tl.rank_genes_groups(new_adata, groupby='cluster_labels', method='wilcoxon', pts = True)
sc.pl.rank_genes_groups(new_adata, n_genes=25, sharey=False, save='_cluster_labels.png')
new_adata.write('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_neighbors_as_genes.h5ad')

aggregated_df = df_concat.groupby('cluster_labels').mean()

# Scale the aggregated data
scaler = StandardScaler()
aggregated_normalized = scaler.fit_transform(aggregated_df)
print(aggregated_df)
print(aggregated_normalized)
# Compute the linkage matrix using the Ward method
linkage_matrix = linkage(aggregated_normalized, method='ward')

# Plot the dendrogram
plt.figure(figsize=(15, 10))
dendrogram(linkage_matrix, labels=aggregated_df.index.astype(str).tolist(), leaf_rotation=90, leaf_font_size=10)
plt.title('Hierarchical Clustering Dendrogram (Ward)')
plt.xlabel('Cluster label')
plt.ylabel('Distance')
plt.tight_layout()
plt.savefig('/home/bcd/revision_nature/neighbors_neighborhood/plots_script4/dendrogram.png', dpi=900)
plt.close()


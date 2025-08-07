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

sc.settings.figdir = "/home/bcd/revision_nature/iPT_iLOH_ME/plots_script1"  # Update this path as needed
sc.settings.set_figure_params(dpi=900)

df_concat = pd.read_csv('/home/bcd/revision_nature/neighbors_neighborhood/neighborhood_dataframe.csv', index_col=0)
df_concat.index.name = None 
print(df_concat)

new_adata=sc.read_h5ad('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_neighbors_as_genes.h5ad')
print(new_adata.obs['cluster_labels_annotated'].value_counts())
print(new_adata)
obs=['iPT', 'iTAL']
mask=new_adata.obs['annotation_updated'].isin(obs)
new_adata=new_adata[mask].copy()
print(new_adata)

# Subset df_concat to include only indexes that are present in new_adata.obs_names
df_concat = df_concat.loc[df_concat.index.intersection(new_adata.obs_names)].copy()

print(df_concat)
df_concat['glomerular_20']=df_concat['EC_glom_20']+df_concat['Podo_20']+df_concat['MC1_20']+df_concat['PEC_20']
df_concat['glomerular_40']=df_concat['EC_glom_40']+df_concat['Podo_40']+df_concat['MC1_40']+df_concat['PEC_40']
df_concat['glomerular_60']=df_concat['EC_glom_60']+df_concat['Podo_60']+df_concat['MC1_60']+df_concat['PEC_60']
df_concat['glomerular_80']=df_concat['EC_glom_80']+df_concat['Podo_80']+df_concat['MC1_80']+df_concat['PEC_80']
df_concat=df_concat.drop(['EC_glom_20', 'EC_glom_40','EC_glom_60','EC_glom_80','Podo_20','Podo_40','Podo_60','Podo_80',
'MC1_20','MC1_40','MC1_60','MC1_80','PEC_20','PEC_40','PEC_60','PEC_80'], axis=1)
df_concat['injured_20']=df_concat['iPT_20']+df_concat['iTAL_20']
df_concat['injured_40']=df_concat['iPT_40']+df_concat['iTAL_40']
df_concat['injured_60']=df_concat['iPT_60']+df_concat['iTAL_60']
df_concat['injured_80']=df_concat['iPT_80']+df_concat['iTAL_80']
df_concat=df_concat.drop(['iPT_20', 'iTAL_20','iPT_40', 'iTAL_40','iPT_60', 'iTAL_60','iPT_80', 'iTAL_80',], axis=1)
df_concat['CD_20']=df_concat['IC A_20']+df_concat['IC B_20']+df_concat['PC_20']
df_concat['CD_40']=df_concat['IC A_40']+df_concat['IC B_40']+df_concat['PC_40']
df_concat['CD_60']=df_concat['IC A_60']+df_concat['IC B_60']+df_concat['PC_60']
df_concat['CD_80']=df_concat['IC A_80']+df_concat['IC B_80']+df_concat['PC_80']
df_concat=df_concat.drop(['PC_20', 'IC A_20','IC B_20','PC_40', 'IC A_40','IC B_40','PC_60', 'IC A_60','IC B_60','PC_80', 'IC A_80','IC B_80', ], axis=1)
df_concat['LOH_20']=df_concat['TAL_20']+df_concat['DTL_ATL_20']
df_concat['LOH_40']=df_concat['TAL_40']+df_concat['DTL_ATL_40']
df_concat['LOH_60']=df_concat['TAL_60']+df_concat['DTL_ATL_60']
df_concat['LOH_80']=df_concat['TAL_80']+df_concat['DTL_ATL_80']
df_concat=df_concat.drop(['TAL_20', 'DTL_ATL_20','TAL_40', 'DTL_ATL_40','TAL_60', 'DTL_ATL_60','TAL_80', 'DTL_ATL_80'], axis=1)
df_concat['EC_20']=df_concat['EC_Peritub_20']+df_concat['EC_Lymph_20']
df_concat['EC_40']=df_concat['EC_Peritub_40']+df_concat['EC_Lymph_40']
df_concat['EC_60']=df_concat['EC_Peritub_60']+df_concat['EC_Lymph_60']
df_concat['EC_80']=df_concat['EC_Peritub_80']+df_concat['EC_Lymph_80']
df_concat=df_concat.drop(['EC_Peritub_20', 'EC_Lymph_20','EC_Peritub_40', 'EC_Lymph_40','EC_Peritub_60', 'EC_Lymph_60','EC_Peritub_80', 'EC_Lymph_80',], axis=1)
print(df_concat)
for column in df_concat:
    print(column)
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
plt.savefig('/home/bcd/revision_nature/iPT_iLOH_ME/plots_script1/pca_variance_curve.png', dpi=900)
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
plt.savefig('/home/bcd/revision_nature/iPT_iLOH_ME/plots_script1/elbow_kmeans.png', dpi=900)
plt.close()

# Perform MiniBatchKMeans clustering
kmeans = MiniBatchKMeans(n_clusters=6, random_state=42)
labels = kmeans.fit_predict(principal_components)
df_concat['cluster_labels'] = labels

new_adata.obs['cluster_labels']=df_concat['cluster_labels'].copy()
new_adata.obs['cluster_labels']=new_adata.obs['cluster_labels'].astype('category')
sc.tl.rank_genes_groups(new_adata, groupby='cluster_labels', method='wilcoxon', pts = True)
sc.pl.rank_genes_groups(new_adata, n_genes=25, sharey=False, save='_cluster_labels.png')
new_adata.write('/home/bcd/revision_nature/iPT_iLOH_ME/xenium_cosmx_neighbors_as_genes_iTAL_iPT_subset.h5ad')

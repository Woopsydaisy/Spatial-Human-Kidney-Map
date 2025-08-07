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
import anndata

##create adata with neighborhood information as gene features for easier downstream analysis
df_concat = pd.read_csv('/home/bcd/revision_nature/immune_ME/neighborhood_dataframe_immune_annotations.csv', index_col=0)
df_concat.index.name = None 

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
    'immune_cell_annotation_combined','immune_annotation_with_other_celltypes'
]
for obs in obs_to_map:
    new_adata.obs[obs] = adata_subset.obs[obs].copy()

# Copy UMAP coordinates
new_adata.obsm['X_umap'] = adata_subset.obsm['X_umap'].copy()
# Save AnnData

print(new_adata)
new_adata=new_adata[new_adata.obs['immune_cell_annotation_combined']!='Unknown'].copy()
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

df_concat['healthy_tubule_20']=df_concat['CNT_20']+df_concat['PT_20']+df_concat['DCT_20']+df_concat['CD_20']+df_concat['LOH_20']
df_concat['healthy_tubule_40']=df_concat['CNT_40']+df_concat['PT_40']+df_concat['DCT_40']+df_concat['CD_40']+df_concat['LOH_40']
df_concat['healthy_tubule_60']=df_concat['CNT_60']+df_concat['PT_60']+df_concat['DCT_60']+df_concat['CD_60']+df_concat['LOH_60']
df_concat['healthy_tubule_80']=df_concat['CNT_80']+df_concat['PT_80']+df_concat['DCT_80']+df_concat['CD_80']+df_concat['LOH_80']
df_concat=df_concat.drop(['CNT_20', 'PT_20','DCT_20', 'CD_20','LOH_20',
'CNT_40', 'PT_40','DCT_40', 'CD_40','LOH_40',
'CNT_60', 'PT_60','DCT_60', 'CD_60','LOH_60',
'CNT_80', 'PT_80','DCT_80', 'CD_80','LOH_80',], axis=1)

print(df_concat)
for column in df_concat:
    print(column)
# Normalize the data
scaler = StandardScaler()
features_normalized = scaler.fit_transform(df_concat)

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
plt.savefig('/home/bcd/revision_nature/immune_ME/plots_script3/elbow_kmeans.png', dpi=900)
plt.close()

# Perform MiniBatchKMeans clustering
kmeans = MiniBatchKMeans(n_clusters=10, random_state=42)
labels = kmeans.fit_predict(principal_components)
df_concat['cluster_labels'] = labels

new_adata.obs['cluster_labels']=df_concat['cluster_labels'].copy()
new_adata.obs['cluster_labels']=new_adata.obs['cluster_labels'].astype('category')
cell_identities = {0: 'Inj. Tubular Immune ME', 1: 'Residential Immune ME', 2: 'B predom. Immune ME', 3: 'Vascular Immune ME',
                   4: 'T predom. Immune ME', 5: 'T predom. Immune ME', 6:'Glomerular Immune ME',7:'Fibro Immune ME',8:'Glomerular Immune ME',9:'T predom. Immune ME'}
new_adata.obs["cluster_labels_annotated_immune"] = new_adata.obs['cluster_labels'].map(cell_identities).astype('category')
new_adata.write('/home/bcd/revision_nature/immune_ME/xenium_cosmx_neighbors_as_genes_immune_subset.h5ad')

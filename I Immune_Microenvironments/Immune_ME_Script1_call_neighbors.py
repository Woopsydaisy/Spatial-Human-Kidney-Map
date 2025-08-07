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

sc.settings.figdir = "/home/bcd/revision_nature/immune_ME/plots_script1"  # Update this path as needed
sc.settings.set_figure_params(dpi=900)

adata = sc.read_h5ad('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad')

# Make a copy to avoid modifying the original annotation directly (optional but good practice)
adata.obs['immune_annotation_with_other_celltypes'] = adata.obs['immune_cell_annotation_combined'].astype(str)
# Replace 'Unknown' with corresponding values from 'annotation_updated'
mask = adata.obs['immune_annotation_with_other_celltypes'] == 'Unknown'
adata.obs.loc[mask, 'immune_annotation_with_other_celltypes'] = adata.obs.loc[mask, 'annotation_updated'].astype(str)
adata.obs['immune_annotation_with_other_celltypes'] = adata.obs['immune_annotation_with_other_celltypes'].astype('category')

all_dfs = []  # List to store DataFrames for each AnnData object

for i in adata.obs["unique_sample_identifier"].unique():
    adata_sample = adata[adata.obs["unique_sample_identifier"] == i]
    # Create a temporary DataFrame for the current AnnData object
    temp_df = pd.DataFrame(index=adata_sample.obs.index)

    for pop in adata_sample.obs['immune_annotation_with_other_celltypes'].unique():
        neighbor_type_indices = np.where(adata_sample.obs['immune_annotation_with_other_celltypes'] == pop)[0]
        number_neighbors = pd.DataFrame(adata_sample.obsp["20_micron_connectivities"][:, neighbor_type_indices].sum(axis=1), index=adata_sample.obs.index)
        temp_df[pop] = number_neighbors[0]  # Assure correct formatting
    # Append the temporary DataFrame to the list
    all_dfs.append(temp_df)

# Concatenate the DataFrames from each AnnData object
df_for_mapping_20 = pd.concat(all_dfs, axis=0)

# Optionally, sort the index if needed
df_for_mapping_20.sort_index(inplace=True)
df_for_mapping_20=df_for_mapping_20.fillna(0)

##now 40um
all_dfs = []

for i in adata.obs["unique_sample_identifier"].unique():
    adata_sample = adata[adata.obs["unique_sample_identifier"] == i]
    # Create a temporary DataFrame for the current AnnData object
    temp_df = pd.DataFrame(index=adata_sample.obs.index)

    for pop in adata_sample.obs['immune_annotation_with_other_celltypes'].unique():
        neighbor_type_indices = np.where(adata_sample.obs['immune_annotation_with_other_celltypes'] == pop)[0]
        number_neighbors = pd.DataFrame(adata_sample.obsp["40_micron_connectivities"][:, neighbor_type_indices].sum(axis=1), index=adata_sample.obs.index)
        temp_df[pop] = number_neighbors[0]  # Assure correct formatting
    # Append the temporary DataFrame to the list
    all_dfs.append(temp_df)

# Concatenate the DataFrames from each AnnData object
df_for_mapping_40 = pd.concat(all_dfs, axis=0)

# Optionally, sort the index if needed
df_for_mapping_40.sort_index(inplace=True)
df_for_mapping_40=df_for_mapping_40.fillna(0)


##now 60um
all_dfs = []
for i in adata.obs["unique_sample_identifier"].unique():
    adata_sample = adata[adata.obs["unique_sample_identifier"] == i]
    # Create a temporary DataFrame for the current AnnData object
    temp_df = pd.DataFrame(index=adata_sample.obs.index)

    for pop in adata_sample.obs['immune_annotation_with_other_celltypes'].unique():
        neighbor_type_indices = np.where(adata_sample.obs['immune_annotation_with_other_celltypes'] == pop)[0]
        number_neighbors = pd.DataFrame(adata_sample.obsp["60_micron_connectivities"][:, neighbor_type_indices].sum(axis=1), index=adata_sample.obs.index)
        temp_df[pop] = number_neighbors[0]  # Assure correct formatting
    # Append the temporary DataFrame to the list
    all_dfs.append(temp_df)

# Concatenate the DataFrames from each AnnData object
df_for_mapping_60 = pd.concat(all_dfs, axis=0)

# Optionally, sort the index if needed
df_for_mapping_60.sort_index(inplace=True)
df_for_mapping_60=df_for_mapping_60.fillna(0)


##now 80um
all_dfs = []
for i in adata.obs["unique_sample_identifier"].unique():
    adata_sample = adata[adata.obs["unique_sample_identifier"] == i]
    # Create a temporary DataFrame for the current AnnData object
    temp_df = pd.DataFrame(index=adata_sample.obs.index)

    for pop in adata_sample.obs['immune_annotation_with_other_celltypes'].unique():
        neighbor_type_indices = np.where(adata_sample.obs['immune_annotation_with_other_celltypes'] == pop)[0]
        number_neighbors = pd.DataFrame(adata_sample.obsp["80_micron_connectivities"][:, neighbor_type_indices].sum(axis=1), index=adata_sample.obs.index)
        temp_df[pop] = number_neighbors[0]  # Assure correct formatting
    # Append the temporary DataFrame to the list
    all_dfs.append(temp_df)

# Concatenate the DataFrames from each AnnData object
df_for_mapping_80 = pd.concat(all_dfs, axis=0)

# Optionally, sort the index if needed
df_for_mapping_80.sort_index(inplace=True)
df_for_mapping_80=df_for_mapping_80.fillna(0)


df_for_mapping_80=df_for_mapping_80-df_for_mapping_60
df_for_mapping_60=df_for_mapping_60-df_for_mapping_40
df_for_mapping_40=df_for_mapping_40-df_for_mapping_20

df_for_mapping_80['sum']=df_for_mapping_80.sum(axis=1)
df_for_mapping_60['sum']=df_for_mapping_60.sum(axis=1)
df_for_mapping_40['sum']=df_for_mapping_40.sum(axis=1)
df_for_mapping_20['sum']=df_for_mapping_20.sum(axis=1)
print(df_for_mapping_80)
# Step 1: Filter out rows where 'sum' == 0
df_for_mapping_20 = df_for_mapping_20[df_for_mapping_20['sum'] > 0]
df_for_mapping_40 = df_for_mapping_40[df_for_mapping_40['sum'] > 0]
df_for_mapping_60 = df_for_mapping_60[df_for_mapping_60['sum'] > 0]
df_for_mapping_80 = df_for_mapping_80[df_for_mapping_80['sum'] > 0]
#filter 5% lowest

threshold = df_for_mapping_80['sum'].quantile(0.05)
# Keep only rows with 'sum' greater than this threshold
df_for_mapping_80 = df_for_mapping_80[df_for_mapping_80['sum'] > threshold]

# Step 2: Find common indices
common_index = (
    df_for_mapping_20.index
    .intersection(df_for_mapping_40.index)
    .intersection(df_for_mapping_60.index)
    .intersection(df_for_mapping_80.index)
)
print('after filter:')
print(df_for_mapping_80)
df_for_mapping_80 = df_for_mapping_80.div(df_for_mapping_80['sum'], axis=0)
df_for_mapping_60 = df_for_mapping_60.div(df_for_mapping_60['sum'], axis=0)
df_for_mapping_40 = df_for_mapping_40.div(df_for_mapping_40['sum'], axis=0)
df_for_mapping_20 = df_for_mapping_20.div(df_for_mapping_20['sum'], axis=0)
print(df_for_mapping_80)
df_for_mapping_80=df_for_mapping_80.drop('sum', axis=1)
df_for_mapping_60=df_for_mapping_60.drop('sum', axis=1)
df_for_mapping_40=df_for_mapping_40.drop('sum', axis=1)
df_for_mapping_20=df_for_mapping_20.drop('sum', axis=1)
print(df_for_mapping_80)

df_for_mapping_20.columns = [f"{col}_20" for col in df_for_mapping_20.columns]
df_for_mapping_40.columns = [f"{col}_40" for col in df_for_mapping_40.columns]
df_for_mapping_60.columns = [f"{col}_60" for col in df_for_mapping_60.columns]
df_for_mapping_80.columns = [f"{col}_80" for col in df_for_mapping_80.columns]

# Subset all mapping DataFrames to only shared cells
df_for_mapping_20 = df_for_mapping_20.loc[common_index]
df_for_mapping_40 = df_for_mapping_40.loc[common_index]
df_for_mapping_60 = df_for_mapping_60.loc[common_index]
df_for_mapping_80 = df_for_mapping_80.loc[common_index]

df_concat = pd.concat(
    [df_for_mapping_20, df_for_mapping_40, df_for_mapping_60, df_for_mapping_80],
    axis=1
)

print(df_concat.shape)
print(df_concat.head())

rows_with_nan = df_concat[df_concat.isna().any(axis=1)]
print(rows_with_nan)
df_concat=df_concat.fillna(0)
df_concat.to_csv('/home/bcd/revision_nature/immune_ME/neighborhood_dataframe_immune_annotations.csv', index=True)
df_concat=pd.read_csv('/home/bcd/revision_nature/immune_ME/neighborhood_dataframe_immune_annotations.csv')
print(df_concat)

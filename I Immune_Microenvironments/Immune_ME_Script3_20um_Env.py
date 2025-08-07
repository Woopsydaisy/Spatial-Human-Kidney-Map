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
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from scipy.sparse import csr_matrix
import gc


adata=sc.read_h5ad('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad')
print(adata)
all_dfs = []

for i in adata.obs["unique_sample_identifier"].unique():
    adata_sample = adata[adata.obs["unique_sample_identifier"] == i]
    connectivity_matrix = csr_matrix(adata_sample.obsp["20_micron_connectivities"])
    # Create a temporary DataFrame for the current AnnData object
    temp_df = pd.DataFrame(index=adata_sample.obs.index)

    # Extract neighbor cell IDs for each cell
    neighbor_cell_ids = []
    for cell_index in range(connectivity_matrix.shape[0]):
        # Extract the indices of non-zero elements (neighbor indices) in the row
        neighbor_indices = connectivity_matrix[cell_index, :].nonzero()[1]

        neighbor_ids = adata_sample.obs.index[neighbor_indices].tolist()
        neighbor_cell_ids.append(neighbor_ids)

    temp_df["neighbor_ids"] = neighbor_cell_ids

    # Append the temporary DataFrame to the list
    all_dfs.append(temp_df)

# Concatenate the DataFrames from each AnnData object
df_neighbor_ids = pd.concat(all_dfs, axis=0)

# Optionally, sort the index if needed
df_neighbor_ids.sort_index(inplace=True)

# Display the DataFrame
print(df_neighbor_ids)

df_neighbor_ids['cluster_labels']=adata.obs['immune_ME']

subset_df=df_neighbor_ids[df_neighbor_ids['cluster_labels']!='Unknown']
subset_df

# DataFrame for the index and cluster_label
index_cluster_df = subset_df.reset_index().rename(columns={'index': 'cell_id'})[['cell_id', 'cluster_labels']]

# DataFrame for the neighbor_ids and cluster_label
neighbors_cluster_df = subset_df.reset_index(drop=True)[['neighbor_ids', 'cluster_labels']]


# Explode the 'neighbor_ids' column to have each neighbor_id in its own row
exploded_neighbors_cluster_df = neighbors_cluster_df.explode('neighbor_ids')

# Drop duplicates in the 'neighbor_ids' column
exploded_neighbors_cluster_df = exploded_neighbors_cluster_df.drop_duplicates(subset='neighbor_ids')


# Rename 'neighbor_ids' column to 'cell_id' in exploded_neighbors_cluster_df for consistency
exploded_neighbors_cluster_df = exploded_neighbors_cluster_df.rename(columns={'neighbor_ids': 'cell_id'})

# Concatenate the DataFrames vertically
combined_df = pd.concat([index_cluster_df, exploded_neighbors_cluster_df], ignore_index=True)

combined_df = combined_df.drop_duplicates(subset=['cell_id'])
combined_df

cluster_label_dict = combined_df.set_index('cell_id')['cluster_labels'].to_dict()
adata.obs['Immune_ME_20um'] = adata.obs_names.map(cluster_label_dict)

adata.obs['Immune_ME_20um']=adata.obs['Immune_ME_20um'].astype('category')
adata.obs['Immune_ME_20um'] = adata.obs['Immune_ME_20um'].cat.add_categories(['Unknown'])
adata.obs["Immune_ME_20um"]=adata.obs["Immune_ME_20um"].fillna('Unknown')

print(adata.obs['Immune_ME_20um'].value_counts())
print(adata.obs['Immune_ME_20um'].value_counts())
adata.write('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad')

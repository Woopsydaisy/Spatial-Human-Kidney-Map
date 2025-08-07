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
import umap
import seaborn as sns
from scipy.sparse import csr_matrix

sc.settings.figdir = "/home/bcd/revision_nature/iPT_iLOH_ME/plots_script4"  # Update this path as needed
sc.settings.set_figure_params(dpi=900)

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

df_neighbor_ids['cluster_labels']=adata.obs['iPT_iLOH_ME']

subset_df=df_neighbor_ids[df_neighbor_ids['cluster_labels']!='Unknown']
subset_df

# DataFrame for the index and cluster_label
index_cluster_df = subset_df.reset_index().rename(columns={'index': 'cell_id'})[['cell_id', 'cluster_labels']]

# DataFrame for the neighbor_ids and cluster_label
neighbors_cluster_df = subset_df.reset_index(drop=True)[['neighbor_ids', 'cluster_labels']]

print(index_cluster_df.head())  # Preview the DataFrame with cell_id and cluster_label
print(neighbors_cluster_df.head())  # Preview the DataFrame with neighbor_ids and cluster_label


# Explode the 'neighbor_ids' column to have each neighbor_id in its own row
exploded_neighbors_cluster_df = neighbors_cluster_df.explode('neighbor_ids')

# Drop duplicates in the 'neighbor_ids' column
exploded_neighbors_cluster_df = exploded_neighbors_cluster_df.drop_duplicates(subset='neighbor_ids')

print(exploded_neighbors_cluster_df)  # Preview the exploded and deduplicated DataFrame


# Rename 'neighbor_ids' column to 'cell_id' in exploded_neighbors_cluster_df for consistency
exploded_neighbors_cluster_df = exploded_neighbors_cluster_df.rename(columns={'neighbor_ids': 'cell_id'})

# Concatenate the DataFrames vertically
combined_df = pd.concat([index_cluster_df, exploded_neighbors_cluster_df], ignore_index=True)

print(combined_df)  # Preview the combined DataFrame


combined_df = combined_df.drop_duplicates(subset=['cell_id'])
combined_df
print(combined_df)
cluster_label_dict = combined_df.set_index('cell_id')['cluster_labels'].to_dict()
adata.obs['iPT_iLOH_ME_20um'] = adata.obs_names.map(cluster_label_dict)

adata.obs['iPT_iLOH_ME_20um']=adata.obs['iPT_iLOH_ME_20um'].astype('category')
adata.obs['iPT_iLOH_ME_20um'] = adata.obs['iPT_iLOH_ME_20um'].cat.add_categories(['Unknown'])
adata.obs["iPT_iLOH_ME_20um"]=adata.obs["iPT_iLOH_ME_20um"].fillna('Unknown')

print(adata.obs['iPT_iLOH_ME'].value_counts())
print(adata.obs['iPT_iLOH_ME_20um'].value_counts())
adata.write('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad')

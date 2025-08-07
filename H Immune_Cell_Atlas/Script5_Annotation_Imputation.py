from anndata import AnnData
import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import csr_matrix, vstack
import gc

adata=sc.read_h5ad('/home/bcd/revision_nature/immune_cell_atlas_xenium/immune_cells_xenium_cosmx_sn_postintegration.h5ad')
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_0_5_postscanvi')

cell_identities = {'0': 'CD8+', '1': 'Macro', '2': 'Macro', '3': 'CD4+', '4': 'Macro', '5': 'Macro', '6': 'CD4+', '7': 'B', 
'8': 'Macro', '9':'Macro','10':'CD8+','11':'Neutrophil','12':'Macro','13':'Plasma','14':'NK','15':'CD4+', 
'16':'Baso_Mast','17':'cDC','18':'pDC','19':'mDC','20':'?','21':'Macro'}

adata.obs["annotation_postscanvi_immune"] = adata.obs['leiden_0_5_postscanvi'].map(cell_identities).astype('category')
sc.pl.umap(adata, color = "annotation_postscanvi_immune", save=f'_annotation_postscanvi_immune.png')


xenium_data_all_cells_all_genes=sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/integration_xenium/all_xenium_raw_combined.h5ad')
# Ensure only cell names present in both datasets are used
filtered_cells_xenium = adata.obs_names.intersection(xenium_data_all_cells_all_genes.obs_names)
# Now subset using the intersection
xenium_data_filtered = xenium_data_all_cells_all_genes[filtered_cells_xenium, :].copy()


cosmx_data_all_cells_all_genes=sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/integration_cosmx/all_cosmx_all_genes_revision_raw.h5ad')
# Ensure only cell names present in both datasets are used
filtered_cells_cosmx = adata.obs_names.intersection(cosmx_data_all_cells_all_genes.obs_names)
# Now subset using the intersection
cosmx_data_filtered = cosmx_data_all_cells_all_genes[filtered_cells_cosmx, :].copy()

adata_concat = ad.concat([xenium_data_filtered, cosmx_data_filtered], join='outer', index_unique=None)


sn=sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/Human_extended_annotated_subset.h5ad')
filtered_cells_sn=adata.obs_names.intersection(sn.obs_names)
sn_data_filtered=sn[filtered_cells_sn, :].copy()
sn_data_filtered.layers['counts']=sn_data_filtered.X

adata_concat = ad.concat([sn_data_filtered, adata_concat], join='outer', index_unique=None)

# Reorder adata_concat so that it matches adata
adata_concat = adata_concat[adata.obs_names.tolist(), :].copy()
# Check if all cells in `adata` are in `adata_concat` and in the same order
same_order = (adata.obs_names == adata_concat.obs_names[:adata.n_obs]).all()
print("Same order as adata:", same_order)

adata_concat.obsm['X_scANVI']=adata.obsm['X_scANVI'].copy()
adata_concat.obsm['X_umap']=adata.obsm['X_umap'].copy()

adata_concat.X=adata_concat.layers['counts'].copy()

Cosmx_cells_mask = (adata_concat.obs['tech'] != 'SN_SEQ')
snRNA_cells_mask = (adata_concat.obs['tech'] == 'SN_SEQ')

CosMx_index = np.where(Cosmx_cells_mask)[0]
snRNA_index = np.where(snRNA_cells_mask)[0]

nn = NearestNeighbors(n_neighbors=15, metric='euclidean')
nn.fit(adata_concat.obsm["X_scANVI"][snRNA_cells_mask])

distances_all_to_non_Cosmx, indices_all_to_non_Cosmx = nn.kneighbors(adata_concat.obsm["X_scANVI"][Cosmx_cells_mask])

print("Min distance:", distances_all_to_non_Cosmx.min())
print("Max distance:", distances_all_to_non_Cosmx.max())
min_nonzero = distances_all_to_non_Cosmx[distances_all_to_non_Cosmx > 0].min()
print(f"Smallest non-zero distance: {min_nonzero}")

adata_concat.X = adata_concat.layers["counts"].copy()
expression_df = adata_concat.to_df()
scRNA_df = expression_df.iloc[snRNA_index].T
expression_df = adata_concat.to_df()
CosMx_df = expression_df.iloc[CosMx_index]
CosMx_df.loc[:, :] = 0  ##for faster operation 0 

scRNA_array = scRNA_df.to_numpy()
CosMx_array = CosMx_df.to_numpy()
print(scRNA_array)
print(scRNA_array.shape)
print(CosMx_array)
print(CosMx_array.shape)

epsilon = 3e-08  # Just above the smallest non-zero distance
print('distances_all:')
print(distances_all_to_non_Cosmx)
distances_safe = distances_all_to_non_Cosmx + epsilon
print('distances_safe:')
print(distances_safe)

adj = (distances_safe ** -2.0).sum(axis=1)
affinity_array = (distances_safe ** -2)
print('adj')
print(adj)
print('affinity_array')
print(affinity_array)
print('indices:')
print(indices_all_to_non_Cosmx)

for i in range(len(indices_all_to_non_Cosmx)):
    CosMx_array[i, :] = (scRNA_array[:, indices_all_to_non_Cosmx[i]] * affinity_array[i]).sum(axis=1) / adj[i]

print('done')
print(CosMx_array)

# Total number of elements
total_values = CosMx_array.size

# Count values by conditions
count_zero = np.sum(CosMx_array == 0)
count_lt_001 = np.sum((CosMx_array > 0) & (CosMx_array < 0.001))
count_lt_01 = np.sum((CosMx_array >= 0.001) & (CosMx_array < 0.01))
count_gte_01 = np.sum(CosMx_array >= 0.01)
count_gte_001 = np.sum(CosMx_array >= 0.001)

# Print counts
print(f"Total values: {total_values:,}")
print(f"== 0: {count_zero:,}")
print(f"> 0 and < 0.001: {count_lt_001:,}")
print(f">= 0.001 and < 0.01: {count_lt_01:,}")
print(f">= 0.01: {count_gte_01:,}")
print(f">= 0.001: {count_gte_001:,}")

CosMx_array[CosMx_array < 0.01] = 0

# Count values by conditions
count_zero = np.sum(CosMx_array == 0)
count_lt_001 = np.sum((CosMx_array > 0) & (CosMx_array < 0.001))
count_lt_01 = np.sum((CosMx_array >= 0.001) & (CosMx_array < 0.01))
count_gte_01 = np.sum(CosMx_array >= 0.01)
count_gte_001 = np.sum(CosMx_array >= 0.001)

# Print counts
print(f"Total values: {total_values:,}")
print(f"== 0: {count_zero:,}")
print(f"> 0 and < 0.001: {count_lt_001:,}")
print(f">= 0.001 and < 0.01: {count_lt_01:,}")
print(f">= 0.01: {count_gte_01:,}")
print(f">= 0.001: {count_gte_001:,}")

chunk_size = 200000
n_rows = CosMx_array.shape[0]
chunks = []

for start in range(0, n_rows, chunk_size):
    end = min(start + chunk_size, n_rows)
    chunk = CosMx_array[start:end, :]
    
    csr_chunk = csr_matrix(chunk, dtype='float32')
    chunks.append(csr_chunk)
    
    print(f"Processed rows {start} to {end}")
    gc.collect()

# Stack all sparse chunks together
final_csr_matrix = vstack(chunks)
print(f"Final shape: {final_csr_matrix.shape}")
print(final_csr_matrix)
adata_CosMx = adata_concat[adata_concat.obs["tech"] != "SN_SEQ"].copy()
adata_CosMx.layers["scanvi_imputed"] = final_csr_matrix
adata_CosMx.X = adata_CosMx.layers["scanvi_imputed"]
adata_CosMx.write('/home/bcd/revision_nature/immune_cell_atlas_xenium/xenium_cosmx_immune_cell_atlas_imputed.h5ad')

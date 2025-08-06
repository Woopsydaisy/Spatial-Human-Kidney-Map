from anndata import AnnData
import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.neighbors import NearestNeighbors
import gc
import csv
from scipy.sparse import csr_matrix
from scipy.sparse import csr_matrix, vstack


adata=sc.read_h5ad('/home/bcd/revision_nature/imputation_new_atlas/xenium_cosmx_sn_preimputation.h5ad')
adata.X=adata.layers['counts'].copy()
print(adata)
print(adata.X)
print(adata.obs['tech'].value_counts())

Cosmx_cells_mask = (adata.obs['tech'] != 'SN_SEQ')
snRNA_cells_mask = (adata.obs['tech'] == 'SN_SEQ')

CosMx_index = np.where(Cosmx_cells_mask)[0]
snRNA_index = np.where(snRNA_cells_mask)[0]

nn = NearestNeighbors(n_neighbors=15, metric='euclidean')
nn.fit(adata.obsm["X_scANVI"][snRNA_cells_mask])

distances_all_to_non_Cosmx, indices_all_to_non_Cosmx = nn.kneighbors(adata.obsm["X_scANVI"][Cosmx_cells_mask])

print("Min distance:", distances_all_to_non_Cosmx.min())
print("Max distance:", distances_all_to_non_Cosmx.max())
min_nonzero = distances_all_to_non_Cosmx[distances_all_to_non_Cosmx > 0].min()
print(f"Smallest non-zero distance: {min_nonzero}")

adata.X = adata.layers["counts"].copy()
expression_df = adata.to_df()
scRNA_df = expression_df.iloc[snRNA_index].T
expression_df = adata.to_df()
CosMx_df = expression_df.iloc[CosMx_index]
CosMx_df.loc[:] = np.nan ##for faster operation 0? 

scRNA_array = scRNA_df.to_numpy()
CosMx_array = CosMx_df.to_numpy()
print(scRNA_array)
print(scRNA_array.shape)
np.save('/home/bcd/revision_nature/imputation_new_atlas/SC_array.npy', scRNA_array)
print(CosMx_array)
print(CosMx_array.shape)
np.save('/home/bcd/revision_nature/imputation_new_atlas/CosMx_array_preimputed.npy', CosMx_array)

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

np.savez(
    '/home/bcd/revision_nature/imputation_new_atlas/imputation_inputs.npz',
    indices_all_to_non_Cosmx=indices_all_to_non_Cosmx,
    affinity_array=affinity_array,
    adj=adj,
    CosMx_index=CosMx_index,
    snRNA_index=snRNA_index
)


# Load arrays
data = np.load('/home/bcd/revision_nature/imputation_new_atlas/imputation_inputs.npz')
indices_all_to_non_Cosmx = data['indices_all_to_non_Cosmx']
affinity_array = data['affinity_array']
adj = data['adj']
CosMx_index = data['CosMx_index']
snRNA_index = data['snRNA_index']

print('adj')
print(adj)
print('affinity_array')
print(affinity_array)
print('indices:')
print(indices_all_to_non_Cosmx)

CosMx_array = np.load('/home/bcd/revision_nature/imputation_new_atlas/CosMx_array_preimputed.npy')
print(CosMx_array)
scRNA_array = np.load('/home/bcd/revision_nature/imputation_new_atlas/SC_array.npy')
print(scRNA_array)


for i in range(len(indices_all_to_non_Cosmx)):
    CosMx_array[i, :] = (scRNA_array[:, indices_all_to_non_Cosmx[i]] * affinity_array[i]).sum(axis=1) / adj[i]

print('done')
print(CosMx_array)

np.save('/home/bcd/revision_nature/imputation_new_atlas/CosMx_array.npy', CosMx_array)

# Load the array
CosMx_array = np.load('/home/bcd/revision_nature/imputation_new_atlas/CosMx_array.npy')

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

np.save('/home/bcd/revision_nature/imputation_new_atlas/CosMx_array.npy', CosMx_array)

# Load the array
CosMx_array = np.load('/home/bcd/revision_nature/imputation_new_atlas/CosMx_array.npy', mmap_mode='r')

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
adata=sc.read_h5ad('/home/bcd/revision_nature/imputation_new_atlas/xenium_cosmx_sn_preimputation.h5ad')
adata_CosMx = adata[adata.obs["tech"] != "SN_SEQ"].copy()
adata_CosMx.layers["scanvi_imputed"] = final_csr_matrix
adata_CosMx.X = adata_CosMx.layers["scanvi_imputed"]
adata_CosMx.write('/home/bcd/revision_nature/imputation_new_atlas/xenium_cosmx_sn_imputed.h5ad')

sc.pl.umap(adata_CosMx, color = "LRP2", color_map='coolwarm', vmax=1, save='_LRP2.png')
sc.pl.umap(adata_CosMx, color = "SLC8A1", color_map='coolwarm', vmax=5, save='_SLC8A1.png')
sc.pl.umap(adata_CosMx, color = "SLC12A3", color_map='coolwarm', vmax=5, save='_SLC12A3.png')
sc.pl.umap(adata_CosMx, color = "MEIS2", color_map='coolwarm', vmax=1, save='_MEIS2.png')
sc.pl.umap(adata_CosMx, color = "EMCN", color_map='coolwarm', vmax=1, save='_EMCN.png')
sc.pl.umap(adata_CosMx, color = "LDB2", color_map='coolwarm', vmax=1, save='_LDB2.png')
sc.pl.umap(adata_CosMx, color = "CALD1", color_map='coolwarm', vmax=1, save='_CALD1.png')
sc.pl.umap(adata_CosMx, color = "SLC26A7", color_map='coolwarm', vmax=3, save='_SLC26A7.png')
sc.pl.umap(adata_CosMx, color = "SLC4A9", color_map='coolwarm', vmax=5, save='_SLC4A9.png')
sc.pl.umap(adata_CosMx, color = "IKZF1", color_map='coolwarm', vmax=1, save='_IKZF1.png')
sc.pl.umap(adata_CosMx, color = "GATA3", color_map='coolwarm', vmax=1, save='_GATA3.png')
sc.pl.umap(adata_CosMx, color = "CFH", color_map='coolwarm', vmax=1, save='_CFH.png')
sc.pl.umap(adata_CosMx, color = "CUBN", color_map='coolwarm', vmax=1, save='_CUBN.png')
sc.pl.umap(adata_CosMx, color = "ACSM2B", color_map='coolwarm', vmax=5, save='_ACSM2B.png')
sc.pl.umap(adata_CosMx, color = "UMOD", color_map='coolwarm', vmax=5, save='_UMOD.png')
sc.pl.umap(adata_CosMx, color = "SLC12A1", color_map='coolwarm', vmax=5, save='_SLC12A1.png')
sc.pl.umap(adata_CosMx, color = "NPHS2", color_map='coolwarm', vmax=5, save='_NPHS2.png')
sc.pl.umap(adata_CosMx, color = "PLA2R1", color_map='coolwarm', vmax=1, save='_PLA2R1.png')
sc.pl.umap(adata_CosMx, color = "HAVCR1", color_map='coolwarm', vmax=1, save='_HAVCR1.png')



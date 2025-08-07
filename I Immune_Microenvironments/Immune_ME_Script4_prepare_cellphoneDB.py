import csv
import gzip
import os
import gc
import numpy as np
import pandas as pd
import scipy.io
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad
import leidenalg as la
from pathlib import Path

# Load the main object and transfer annotations
adata = sc.read_h5ad('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad')

# Load the imputed object and transfer metadata from `adata`
adata_imputed = sc.read_h5ad('/home/bcd/revision_nature/imputation_new_atlas/xenium_cosmx_sn_imputed_updated.h5ad')
adata_imputed.obs['Immune_ME'] = adata.obs['Immune_ME'].copy()
adata_imputed.obs['Immune_ME_20um'] = adata.obs['Immune_ME_20um'].copy()
adata_imputed.obs['immune_cell_annotation_combined'] = adata.obs['immune_cell_annotation_combined'].copy()
sc.pl.umap(adata, color="Immune_ME", save="_Immune_ME_20um_all.png")
adata_imputed.write('/home/bcd/revision_nature/imputation_new_atlas/xenium_cosmx_sn_imputed_updated.h5ad')

# Remove "Unknown" category
mask = adata_imputed.obs['Immune_ME_20um'] != 'Unknown'
adata_imputed.obs['Immune_ME_20um'] = adata_imputed.obs['Immune_ME_20um'].astype('category')
adata_imputed.obs['Immune_ME_20um'].cat.remove_categories(['Unknown'])
adata_subset = adata_imputed[mask]
sc.pl.umap(adata_subset, color="Immune_ME", save="_Immune_ME_20um_no_Unknown.png")
del adata_imputed
gc.collect()

# Subset for cells without immune annotation
mask = adata_subset.obs['immune_cell_annotation_combined'] == 'Unknown'
adata_subset_subset = adata_subset[mask]
sc.pl.umap(adata_subset_subset, color="Immune_ME", save="_Immune_ME_20um_no_Immune.png")
del adata_subset
gc.collect()

# Load imputed immune-only object and prepare it
adata_imputed_immune = sc.read_h5ad('/home/bcd/revision_nature/immune_cell_atlas_xenium/xenium_cosmx_immune_cell_atlas_imputed.h5ad')
adata_imputed_immune.obs['Immune_ME'] = adata_imputed_immune.obs['immune_ME'].copy()
del adata_imputed_immune.obs['immune_ME']
adata_imputed_immune.obs['Immune_ME_20um'] = adata_imputed_immune.obs['Immune_ME'].copy()
adata_imputed_immune.write('/home/bcd/revision_nature/immune_cell_atlas_xenium/xenium_cosmx_immune_cell_atlas_imputed.h5ad')

# Filter unknown immune ME
adata_imputed_immune_subset = adata_imputed_immune[adata_imputed_immune.obs['Immune_ME_20um'] != 'Unknown'].copy()

# Concatenate epithelial and immune subsets
adata_concat = ad.concat([adata_subset_subset, adata_imputed_immune_subset], join='outer', index_unique=None)
adata_concat.write('/home/bcd/revision_nature/immune_ME/immune_ME_20um_object_for_cellphone.h5ad')

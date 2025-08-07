from anndata import AnnData
import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set figure output directory
sc.settings.figdir = "/home/bcd/revision_nature/immune_cell_atlas_xenium/plots_postintegration_xenium"

# Load data
adata = sc.read_h5ad('/home/bcd/revision_nature/immune_cell_atlas_xenium/immune_cells_xenium_sn_scvi_scanvi4.h5ad')

# Run Leiden clustering
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_0_5_postscanvi')
sc.tl.leiden(adata, resolution=1, key_added='leiden_1_postscanvi')

# Plot UMAPs with Leiden clusters
sc.pl.umap(adata, color="leiden_0_5_postscanvi", save="_leiden_0_5_postscanvi.png")
sc.pl.umap(adata, color="leiden_1_postscanvi", save="_leiden_1_postscanvi.png")

# Gene ranking for cluster markers
adata.X = adata.layers['counts'].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.rank_genes_groups(adata, groupby='leiden_0_5_postscanvi', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save="_leiden_0_5_postscanvi.png")

# Subset and plot by modality
adata_sn = adata[adata.obs['tech'] == 'SN_SEQ'].copy()
sc.pl.umap(adata_sn, color="leiden_0_5_postscanvi", save="_leiden_0_5_postscanvi_sn_only.png")

adata_xenium = adata[adata.obs['tech'] != 'SN_SEQ'].copy()
sc.pl.umap(adata_xenium, color="leiden_0_5_postscanvi", save="_leiden_0_5_postscanvi_xenium_only.png")

# Simplified cluster annotation
cell_identities = {
    '0': 'CD8+',
    '1': 'CD4+',
    '2': 'Macro',
    '3': 'Macro',
    '4': 'Macro',
    '5': 'Macro',
    '6': 'B',
    '7': 'Plasma',
    '8': 'NK',
    '9': 'Treg',
    '10': 'CD4+',
    '11': 'Macro',
    '12': 'Baso_Mast',
    '13': 'cDC',
    '14': 'mDC',
    '15': 'pDC'
}

adata.obs["annotation_postscanvi_immune"] = adata.obs['leiden_0_5_postscanvi'].map(cell_identities).astype('category')

sc.settings.figdir = "/home/bcd/revision_nature/immune_cell_atlas_xenium/plots_preintegration_cosmx"
##here are all immune cell annotations stored
adata=sc.read_h5ad('/home/bcd/revision_nature/immune_cell_atlas_xenium/immune_cells_xenium_sn_scvi_scanvi.h5ad')
print(adata)
##this is the integrated CosMx Xenium SN from which we extract all immune cells
adata_all=sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/integration_combined/revision_atlas_cosmx_xenium_ver1_scvi_scanvi.h5ad')
print(adata_all)
immune_cell_subset=adata_all[adata_all.obs['annotation_postscanvi_level2']=='Immune'].copy()
immune_cell_subset.obs['annotation_postscanvi_immune']=adata.obs['annotation_postscanvi_immune'].copy()
immune_cell_subset.obs['annotation_for_scanvi']=immune_cell_subset.obs['annotation_postscanvi_immune'].copy()
immune_cell_subset.obs['annotation_for_scanvi']=immune_cell_subset.obs['annotation_for_scanvi'].fillna('Unknown')
immune_cell_subset.write('/home/bcd/revision_nature/immune_cell_atlas_xenium/xenium_cosmx_sn_prescvi.h5ad')

from anndata import AnnData
import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

adata=sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/integration_combined/revision_atlas_cosmx_xenium_ver1_scvi_scanvi.h5ad')
sc.settings.figdir = "/home/bcd/revision_nature/integration_new_atlas/integration_combined/plots_scanvi1"  # Update this path as needed
sc.settings.set_figure_params(dpi=900)

#sc.tl.leiden(adata, resolution = 1, key_added='leiden_postscanvi_1')

cell_identities = {'0': 'Fibroblast', '1': 'TAL', '2': 'EC_Peritub', '3': 'iTAL', '4': 'TAL', '5': 'CNT', '6': 'PT', '7': 'PT', 
'8': 'PT', '9':'iPT','10':'PC','11':'Immune','12':'PT','13':'PC','14':'IC A','15':'PT','16':'Fibroblast',
'17':'DCT','18':'iPT','19':'PT','20':'iPT','21': 'EC_Peritub', '22': 'Immune', '23': 'PT',
'24': 'PT', '25': 'PT', '26': 'DCT', '27': 'VSMC', '28': 'Podo', '29':'Fibroblast','30':'TAL',
'31': 'EC_glom', '32': 'IC B', '33': 'Fibroblast', '34': 'EC_DVR', '35': 'CNT', '36': 'PEC', '37': 'DTL_ATL', 
'38': 'MC1', '39':'Immune','40':'EC_Lymph'}

adata.obs["annotation_postscanvi_level2"] = adata.obs['leiden_postscanvi_1'].map(cell_identities).astype('category')
sc.pl.umap(adata, color = "annotation_postscanvi_level2", save=f'_annotation_postscanvi_level2.png')

adata_cosmx=adata[adata.obs['tech']=='CosMx'].copy()
adata_subset_cosmx_annotated=adata_cosmx[adata_cosmx.obs['needs_mapping']=='No'].copy()
print(adata_subset_cosmx_annotated)
for annotation in adata_subset_cosmx_annotated.obs['annotation_postscanvi_level2'].unique():
    subset=adata_subset_cosmx_annotated[adata_subset_cosmx_annotated.obs['annotation_postscanvi_level2']==annotation].copy()
    sc.pl.umap(subset, color = "original_cosmx_annotation", title=f'new:{annotation} with old plotted',save=f'_annotation_postscanvi_level2_with_original_cosmx_annotation_group_{annotation}.png')

adata.write('/home/bcd/revision_nature/integration_new_atlas/integration_combined/revision_atlas_cosmx_xenium_ver1_scvi_scanvi.h5ad')

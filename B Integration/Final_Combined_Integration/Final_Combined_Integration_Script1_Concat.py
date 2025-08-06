from anndata import AnnData
import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

xenium_data_all_cells_all_genes=sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/integration_xenium/all_xenium_raw_combined.h5ad')
print(xenium_data_all_cells_all_genes)
xenium_data_integrated=sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/integration_xenium_v3/revision_atlas_xenium_filter_scvi_scanvi.h5ad')
print(xenium_data_integrated)
xenium_data_integrated=xenium_data_integrated[xenium_data_integrated.obs['tech']=='Xenium']
filtered_cells = xenium_data_integrated.obs_names
xenium_data_filtered = xenium_data_all_cells_all_genes[filtered_cells, :].copy()
# Map annotations from the integrated dataset
xenium_data_filtered.obs['annotation_final_level1'] = xenium_data_integrated.obs['annotation_final_level1'].copy()
print(xenium_data_filtered)
print(xenium_data_filtered.obs['annotation_final_level1'].value_counts())

cosmx_data_all_cells_all_genes=sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/integration_cosmx/all_cosmx_all_genes_revision_raw.h5ad')
print(cosmx_data_all_cells_all_genes)
cosmx_data_integrated=sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/integration_cosmx_v3/revision_atlas_cosmx_filter_scvi_scanvi.h5ad')
print(cosmx_data_integrated)
cosmx_data_integrated=cosmx_data_integrated[cosmx_data_integrated.obs['tech']=='CosMx']
filtered_cells = cosmx_data_integrated.obs_names
cosmx_data_filtered = cosmx_data_all_cells_all_genes[filtered_cells, :].copy()
# Map annotations from the integrated dataset
cosmx_data_filtered.obs['annotation_final_level1'] = cosmx_data_integrated.obs['annotation_final_level1'].copy()
print(cosmx_data_filtered)
print(cosmx_data_filtered.obs['annotation_final_level1'].value_counts())

adata_concat = ad.concat([xenium_data_filtered, cosmx_data_filtered], join='outer', index_unique=None)
print(adata_concat)
print(adata_concat.obs['annotation_final_level1'].value_counts())

sn=sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/Human_extended_annotated.h5ad')
sn.obs['annotation_final_level1']=sn.obs['annotation'].copy()
sn.obs['annotation_final_level1']=sn.obs['annotation_final_level1'].replace('DCT_CNT_PC', 'DCT_CNT_CD')
sn.obs['annotation_final_level1']=sn.obs['annotation_final_level1'].replace('iTAL', 'TAL')
sn.obs['annotation_final_level1']=sn.obs['annotation_final_level1'].replace('iPT', 'PT')
sn.obs['tech']='SN_SEQ'
print(sn.X)
sn.layers['counts']=sn.X.copy()
sn.obs['needs_mapping']='No'
sn.obs['original_cosmx_annotation']="Unknown"

adata_concat = ad.concat([sn, adata_concat], join='inner', index_unique=None)
print(adata_concat)
print(adata_concat.obs['tech'].value_counts())
print(adata_concat.obs['annotation_final_level1'].value_counts())
adata_concat.write('/home/bcd/revision_nature/integration_new_atlas/integration_combined/xenium_cosmx_sn_prescvi.h5ad')

from anndata import AnnData
import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# Load individual Xenium AnnData objects
adata1 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/HK3626_Xenium_raw.h5ad')
adata2 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/1001_Xenium_raw.h5ad')
adata3 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/1003_Xenium_raw.h5ad')
adata4 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/1004_Xenium_raw.h5ad')
adata5 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/1005_Xenium_raw.h5ad')
adata6 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/1007_Xenium_raw.h5ad')
adata7 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/1008_Xenium_raw.h5ad')
adata8 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/1009_Xenium_raw.h5ad')
adata9 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/1006_Xenium_raw.h5ad')
adata10 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/1013_Xenium_raw.h5ad')
adata11 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/1012_Xenium_raw.h5ad')
adata12 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/1011_Xenium_raw.h5ad')
adata13 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/1010_Xenium_raw.h5ad')
adata14 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/HK2695_Xenium_raw.h5ad')
adata15 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/HK2753_Xenium_raw.h5ad')
adata16 = sc.read_h5ad('/home/bcd/Xenium_Rawdata/adata_objects/HK3106_Xenium_raw.h5ad')

adata_list = [
    adata1, adata2, adata3, adata4, adata5, adata6, adata7, adata8,
    adata9, adata10, adata11, adata12, adata13, adata14, adata15, adata16
]

# Annotate and format all individual datasets
for adata in adata_list:
    print(adata)
    adata.obs['orig_ident'] = adata.obs['sample'].copy()
    print(adata.X)
    adata.layers['counts'] = adata.X.copy()
    adata.obs['tech'] = 'Xenium'
    adata.obs['nCount_RNA'] = adata.obs['transcript_counts'].copy()
    adata.obs['nFeature_RNA'] = adata.obs['n_genes_by_counts'].copy()
    adata.obs['annotation_final_level1'] = "Unknown"
    adata.obs['annotation_final_level1B'] = "Unknown"
    adata.obs['annotation'] = 'Unknown'
    adata.obs['proj'] = 'Xenium_Susztak'
    adata.obs['percent_mt'] = 0
    adata.obs['needs_mapping'] = "Yes"
    adata.obs['original_cosmx_annotation'] = "Unknown"
    adata.obs.index = adata.obs['sample'].astype(str) + '_' + adata.obs.index.astype(str)
    print(adata.obs['proj'])

# Concatenate all Xenium datasets
adata_concat = ad.concat([adata1, adata2], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata3], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata4], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata5], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata6], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata7], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata8], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata9], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata10], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata11], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata12], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata13], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata14], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata15], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata16], join='inner', index_unique=None)

print(adata_concat)

# Save concatenated Xenium data
adata_concat.write('/home/bcd/revision_nature/integration_new_atlas/integration_xenium/all_xenium_raw_combined.h5ad')

# Basic filtering
sc.pp.filter_cells(adata_concat, min_counts=30)
print(adata_concat)

# Load annotated SN-seq reference dataset
sn = sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/Human_extended_annotated.h5ad')
print(sn)

# Annotate and standardize metadata for integration
sn.obs['annotation_final_level1'] = sn.obs['annotation'].copy()
sn.obs['annotation_final_level1'] = sn.obs['annotation_final_level1'].replace('DCT_CNT_PC', 'DCT_CNT_CD')
sn.obs['annotation_final_level1'] = sn.obs['annotation_final_level1'].replace('iTAL', 'TAL')
sn.obs['annotation_final_level1'] = sn.obs['annotation_final_level1'].replace('iPT', 'PT')
sn.obs['tech'] = 'SN_SEQ'
print(sn.X)
sn.layers['counts'] = sn.X.copy()
sn.obs['needs_mapping'] = 'No'
sn.obs['original_cosmx_annotation'] = "Unknown"

# Concatenate Xenium and SN-seq data for scVI
scvi_object = ad.concat([adata_concat, sn], join='inner', index_unique=None)
print(scvi_object)

# Save the combined object for scVI processing
scvi_object.write('/home/bcd/revision_nature/integration_new_atlas/integration_xenium/xenium_atlas_prescvi.h5ad')
print(scvi_object.obs['tech'])
print(scvi_object.obs['proj'])


from anndata import AnnData
import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

adata1=sc.read_h5ad('/home/bcd/revision_nature/new_samples/1061_raw_noImages.h5ad')
adata2=sc.read_h5ad('/home/bcd/revision_nature/new_samples/1062_raw_noImages.h5ad')
adata3=sc.read_h5ad('/home/bcd/revision_nature/new_samples/1063_raw_noImages.h5ad')
adata4=sc.read_h5ad('/home/bcd/revision_nature/new_samples/1064_raw_noImages.h5ad')
adata5=sc.read_h5ad('/home/bcd/revision_nature/new_samples/1068_raw_noImages.h5ad')
adata6=sc.read_h5ad('/home/bcd/revision_nature/new_samples/1070_raw_noImages.h5ad')
adata7=sc.read_h5ad('/home/bcd/revision_nature/new_samples/1071_raw_noImages.h5ad')
adata8=sc.read_h5ad('/home/bcd/revision_nature/new_samples/1072_raw_noImages.h5ad')
adata9=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK2990_raw_noImages.h5ad')
adata10=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3068_raw_noImages.h5ad')
adata11=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3069_raw_noImages.h5ad')
adata12=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3421_raw_noImages.h5ad')
adata13=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3588_raw_noImages.h5ad')
adata14=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3594_raw_noImages.h5ad')
adata15=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3606_raw_noImages.h5ad')
adata16=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3614_raw_noImages.h5ad')
adata17=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK2924_raw_noImages.h5ad')
adata18=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK2989_raw_noImages.h5ad')
adata19=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3591_raw_noImages.h5ad')
adata20=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3612_raw_noImages.h5ad')
adata21=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3647_raw_noImages.h5ad')
adata22=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3616_raw_noImages.h5ad')
adata23=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3623_raw_noImages.h5ad')
adata24=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3469_raw_noImages.h5ad')
adata25=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3070_raw_noImages.h5ad')
adata26=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3474_raw_noImages.h5ad')
adata27=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3066_raw_noImages.h5ad')
adata28=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3624_raw_noImages.h5ad')
adata29=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3626_raw_noImages.h5ad')
adata30=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3631_raw_noImages.h5ad')
adata31=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3535_raw_noImages.h5ad')
adata32=sc.read_h5ad('/home/bcd/revision_nature/new_samples/HK3565_raw_noImages.h5ad')

adata_list=[adata1, adata2, adata3, adata4, adata5, adata6, adata7, adata8, adata9, adata10, adata11, adata12, 
adata13, adata14, adata15, adata16, adata17, adata18, adata19, adata20, adata21, adata22, adata23, adata24, adata25, 
adata26, adata27, adata28, adata29, adata30, adata31, adata32]

for adata in adata_list:
    print(adata)
    adata.obs['orig_ident']=adata.obs['sample'].copy()
    print(adata.X)
    adata.layers['counts']=adata.X.copy()
    adata.obs['tech']='CosMx'
    adata.obs['annotation_final_level1']="Unknown"
    adata.obs['annotation_final_level1B']="Unknown"
    adata.obs['annotation']='Unknown'
    adata.obs['proj']='CosMx_Susztak'
    adata.obs['percent_mt']=0
    adata.obs['needs_mapping']="Yes"
    adata.obs['original_cosmx_annotation']="Unknown"
    adata.obs.index = adata.obs['sample'].astype(str) + '_' + adata.obs.index.astype(str)
    print(adata.obs['proj'])



adata_annotated=sc.read_h5ad('/home/bcd/integration_all_cosmx/allCosmx_allGenes.h5ad') ##these are the first 16 samples
genes_to_remove = ['XIST', 'SRY', 'MALAT1']
# Filter out these genes from the AnnData object
adata_annotated = adata_annotated[:, ~adata_annotated.var_names.isin(genes_to_remove)].copy()
adata_annotated.obs['original_cosmx_annotation']=adata_annotated.obs["annotation_post_scanvi70_broad"].copy()
adata_annotated.obs["annotation_post_scanvi70_broad"].value_counts()
print(adata_annotated.obs['annotation_post_scanvi70_broad'])
cell_identities = {'TAL': 'TAL', 'iTAL': 'TAL', 'DTL_ATL': 'DTL_ATL',
                   'iPT': 'PT', 'PT': 'PT',
                   'IC A': 'IC', 'IC B': 'IC',
                   'Immune':'Immune',
                   'PEC':'PEC',
                   'Podo':'Podo',
                   'EC_glom': 'EC', 'EC_Peritub': 'EC', 'EC_DVR': 'EC',
                   'CNT': 'DCT_CNT_CD', 'DCT': 'DCT_CNT_CD', 'PC': 'DCT_CNT_CD',
                   'Fibroblast':'Stroma','VSMC':'Stroma','MC1':'Stroma'}

adata_annotated.obs["annotation_final_level1"] = adata_annotated.obs['annotation_post_scanvi70_broad'].map(cell_identities).astype('category')
adata_annotated.obs['needs_mapping']='No'
adata_annotated.obs['tech']='CosMx'
adata_annotated.obs["annotation_final_level1B"]='Unknown'
adata_annotated.obs["percent_mt"]=0
adata_annotated.obs["proj"]='CosMx_Susztak'
adata_annotated.obs['orig_ident']=adata_annotated.obs['sample'].copy()

sn=sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/Human_extended_annotated_subset.h5ad')
print(sn)
sn.obs['annotation_final_level1']=sn.obs['annotation'].copy()
sn.obs['annotation_final_level1']=sn.obs['annotation_final_level1'].replace('DCT_CNT_PC', 'DCT_CNT_CD')
sn.obs['annotation_final_level1']=sn.obs['annotation_final_level1'].replace('iTAL', 'TAL')
sn.obs['annotation_final_level1']=sn.obs['annotation_final_level1'].replace('iPT', 'PT')
sn.obs['tech']='SN_SEQ'
print(sn.X)
sn.layers['counts']=sn.X.copy()
sn.obs['needs_mapping']='No'
sn.obs['original_cosmx_annotation']="Unknown"


adata_concat = ad.concat([adata_annotated, adata1], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata2], join='inner', index_unique=None)
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
adata_concat = ad.concat([adata_concat, adata17], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata18], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata19], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata20], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata21], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata22], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata23], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata24], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata25], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata26], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata27], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata28], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata29], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata30], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata31], join='inner', index_unique=None)
adata_concat = ad.concat([adata_concat, adata32], join='inner', index_unique=None)

adata_concat = ad.concat([sn, adata_concat], join='inner', index_unique=None)
print(adata_concat)
print(adata_concat.obs['annotation_final_level1'].value_counts())
print(adata_concat.obs['tech'].value_counts())
sc.pp.filter_cells(adata_concat, min_counts=30)
print(adata_concat)
adata_concat.write('/home/bcd/revision_nature/integration_new_atlas/integration_cosmx/revision_atlas_cosmx_prescvi.h5ad')
print(adata_concat.obs['orig_ident'].value_counts())


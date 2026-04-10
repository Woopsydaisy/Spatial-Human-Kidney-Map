from anndata import AnnData
import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.sparse as sp
import scipy.io
import scipy.sparse as sp

adata=sc.read_h5ad('/home/bcd/revision_nature/B_Cell_Atlas/b_atlas_full_genes_annotated.h5ad')
adata.obs['tech']='SC_SEQ'
adata.obs['proj']=adata.obs['dataid'].copy()
print(adata)

xenium_ckd=sc.read_h5ad('/home/bcd/revision_nature/immune_cell_atlas_xenium/extracted_b_plasma_cells_xenium.h5ad')
print(xenium_ckd)
xenium_ckd.obs['annotation_postscvi']='Unknown'

adata_concat = ad.concat([adata, xenium_ckd], join='inner', index_unique=None)
print('after_concat xenium ckd with SC')
print(adata_concat)

adata_concat.obs['tech']=adata_concat.obs['tech'].replace('Xenium_5k', 'Xenium')
adata_concat.write('/home/bcd/revision_nature/B_Cell_Atlas/xenium_sc_b_cells_prescvi.h5ad')
print(adata_concat.obs['tech'].value_counts())

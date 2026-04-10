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
import squidpy as sq

adata=sc.read_h5ad('/home/bcd/revision_nature/B_Cell_Atlas/xenium_sc_b_cells_scvi_scanvi4.h5ad')
print(adata)

sc.tl.leiden(adata, resolution=1, key_added='leiden_postscvi_1')
adata.write('/home/bcd/revision_nature/B_Cell_Atlas/xenium_sc_b_cells_scvi_scanvi4.h5ad')

adata.obs["leiden_postscvi_1"].value_counts()
cell_identities = {'0': 'Plasma', '1': 'Bmem', '2': 'Bmem', '3': 'Bnaive', '4': 'AtM', '5': 'Bnaive', '6': 'Bmem', '7': 'Bmem', 
'8': 'ACB', '9': 'Bnaive', '10':'Bnaive', '11':'ACB','12':'ACB','13':'Bnaive','14':'GCB','15':'Bnaive','16':'AtM','17':'Bnaive','18':'Bnaive','19':'Plasma'}
adata.obs["annotation_postscanvi"] = adata.obs['leiden_postscvi_1'].map(cell_identities).astype('category')

sc.pl.umap(adata, color = "annotation_postscanvi", legend_loc='on data', legend_fontsize=10, legend_fontoutline=2, frameon = False,save='annotation_postscanvi2_v2.png')

adata.write('/home/bcd/revision_nature/B_Cell_Atlas/xenium_sc_b_cells_scvi_scanvi4.h5ad')

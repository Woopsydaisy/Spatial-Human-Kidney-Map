!pip --quiet install scanpy

import os
import random
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns

!pip uninstall dask -y
!pip uninstall datashader -y

!pip install dask==2022.10.0
!pip install datashader==0.14.4

!pip install pyliger

import pyliger
import scanpy as sc
import numpy as np

adata=sc.read_h5ad('/content/drive/MyDrive/Bernhard/revision_nature/integration_benchmarking_object.h5ad')

batch_cats = adata.obs['tech'].cat.categories

bdata = adata.copy()
# Pyliger normalizes by library size with a size factor of 1
# So here we give it the count data
bdata.X = bdata.layers["counts"]
# List of adata per batch
adata_list = [bdata[bdata.obs['tech'] == b].copy() for b in batch_cats]
for i, ad in enumerate(adata_list):
    ad.uns["sample_name"] = batch_cats[i]
    # Hack to make sure each method uses the same genes
    ad.uns["var_gene_idx"] = np.arange(bdata.n_vars)

liger_data = pyliger.create_liger(adata_list, remove_missing=False, make_sparse=False)

# Hack to make sure each method uses the same genes
liger_data.var_genes = bdata.var_names

pyliger.normalize(liger_data)

pyliger.scale_not_center(liger_data)

pyliger.optimize_ALS(liger_data, k=30)

pyliger.quantile_norm(liger_data)

adata.obsm["LIGER"] = np.zeros((adata.shape[0], liger_data.adata_list[0].obsm["H_norm"].shape[1]))
for i, b in enumerate(batch_cats):
    adata.obsm["LIGER"][adata.obs['tech'] == b] = liger_data.adata_list[i].obsm["H_norm"]

adata.write('/content/drive/MyDrive/Bernhard/revision_nature/integration_benchmarking_object.h5ad')

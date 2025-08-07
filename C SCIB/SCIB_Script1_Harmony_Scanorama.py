!pip install --quiet scanorama
!pip install --quiet harmonypy

import os
import random
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns


adata=sc.read_h5ad('/content/drive/MyDrive/Bernhard/revision_nature/integration_benchmarking_object.h5ad')

adata.X=adata.layers['counts'].copy()
print(adata.X)
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.pca(adata)
adata.obsm['Unintegrated']=adata.obsm['X_pca'].copy()

##first scanorama
sc.external.pp.scanorama_integrate(adata, key='tech')
adata.obsm['Scanorama']=adata.obsm['X_scanorama'].copy()

sc.external.pp.harmony_integrate(adata, key='tech', max_iter_harmony=2000)
adata.obsm['Harmony']=adata.obsm['X_pca_harmony'].copy()
adata.write('/content/drive/MyDrive/Bernhard/revision_nature/integration_benchmarking_object.h5ad')


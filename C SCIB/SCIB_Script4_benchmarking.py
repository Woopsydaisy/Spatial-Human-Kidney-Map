!pip install --quiet scib-metrics
!pip install --quiet scib

import scib
import os
import random
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns

adata.obsm['Harmony']=adata.obsm['X_pca_harmony'].copy()
adata.obsm['scVI']=adata.obsm['X_scVI'].copy()
adata.obsm['scANVI']=adata.obsm['X_scANVI'].copy()
adata.obsm['Scanorama']=adata.obsm['X_scanorama'].copy()
adata.obsm['Unintegrated']=adata.obsm['X_pca'].copy()

bm = Benchmarker(
    adata,
    batch_key="tech",
    label_key="annotation",
    embedding_obsm_keys=["Unintegrated", "Scanorama", "LIGER", "Harmony", "scVI", "scANVI"],
    n_jobs=6,
)
bm.benchmark()

bm.plot_results_table()

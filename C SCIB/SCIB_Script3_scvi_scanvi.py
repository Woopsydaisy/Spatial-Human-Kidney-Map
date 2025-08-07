!pip install scanpy
!pip install --quiet scvi-colab
from scvi_colab import install

install()

import os
import random
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
import scvi

adata=sc.read_h5ad('/content/drive/MyDrive/Bernhard/revision_nature/integration_benchmarking_object.h5ad')

scvi.model.SCVI.setup_anndata(
    adata,
    batch_key="tech",
    layer="counts", categorical_covariate_keys=["sample"],
    continuous_covariate_keys=["nCount_RNA", "percent_mt"])
model = scvi.model.SCVI(adata)
model
vae = scvi.model.SCVI(adata, n_layers=4, n_latent=30, gene_likelihood="nb", dropout_rate=0.1)
vae.train(max_epochs = 600, plan_kwargs={"lr":0.001}, early_stopping = True, early_stopping_patience = 25, batch_size=1024)
model = vae
model.save("/content/drive/MyDrive/Bernhard/revision_nature/benchmark_integration/models/models_scvi", overwrite = True)

latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent
adata.layers["scvi_normalized"] = model.get_normalized_expression(
    library_size=10e4)

adata.write('/content/drive/MyDrive/Bernhard/revision_nature/integration_benchmarking_object.h5ad')

#model = scvi.model.SCVI.load("/content/drive/MyDrive/Bernhard/revision_nature/benchmark_integration/models/models_scvi", adata=adata)

seed = 10
scvi.settings.seed = 10

model_path = "/content/drive/MyDrive/Bernhard/revision_nature/benchmark_integration/models/models_scanvi"

scanvi_model = scvi.model.SCANVI.from_scvi_model(model, labels_key="annotation", unlabeled_category = "Unknown")

scanvi_model.train(max_epochs = 600, plan_kwargs={"lr":0.001}, early_stopping = True, early_stopping_patience = 15, batch_size=512)

scanvi_model.save(model_path, overwrite=True)

adata.obs["C_scANVI"] = scanvi_model.predict(adata)
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)

adata.write('/content/drive/MyDrive/Bernhard/revision_nature/integration_benchmarking_object.h5ad')

import os
from anndata import AnnData
import scanpy as sc
import pandas as pd
import anndata as ad
import scvi
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import torch
from sklearn.neighbors import NearestNeighbors
import gc

os.environ["CUDA_VISIBLE_DEVICES"] = "0"
torch.set_float32_matmul_precision('high')
seed=10

sc.settings.figdir = "/home/bcd/revision_nature/integration_new_atlas/integration_combined/plots_scvi1"  # Update this path as needed

adata=sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/integration_combined/xenium_cosmx_sn_prescvi.h5ad')

print(adata)
print(adata.obs['tech'].value_counts())
print(adata.obs['proj'].value_counts())
print(adata.obs['needs_mapping'].value_counts())
print(adata.layers['counts'])
adata=adata.copy()
scvi.model.SCVI.setup_anndata(adata, layer="counts",
        batch_key="tech",
        categorical_covariate_keys=["proj","orig_ident"],
        continuous_covariate_keys=["nCount_RNA", "percent_mt"])

# Initialize and train the SCVI model
vae = scvi.model.SCVI(adata, n_layers=4, n_latent=30, gene_likelihood="nb", dropout_rate=0.1)
vae.train(max_epochs=1000, plan_kwargs={"lr": 0.001}, early_stopping=True, 
              batch_size=2048, early_stopping_patience=25)

# Save the model with a unique name for each HVG count
model_name = f"/home/bcd/revision_nature/integration_new_atlas/integration_combined/models/scvi/model1"
vae.save(model_name, overwrite=True)

# Get latent representation and normalized expression
latent = vae.get_latent_representation()
adata.obsm["X_scVI"] = latent
adata.layers["scvi_normalized"] = vae.get_normalized_expression(library_size=10e4)

# Compute neighbors and UMAP
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)

adata.write("/home/bcd/revision_nature/integration_new_atlas/integration_combined/revision_atlas_cosmx_xenium_ver1_scvi.h5ad")

sc.pl.umap(adata, color="tech", save="_tech.png")
sc.pl.umap(adata, color="proj", save="_proj.png")
sc.pl.umap(adata, color="annotation_final_level1", save="_annotation_final_level1.png")
sc.pl.umap(adata, color="original_cosmx_annotation", save="_original_cosmx_annotation.png")

for annotation in adata.obs['annotation_final_level1'].unique():
        sc.pl.umap(adata, color="annotation_final_level1", groups=annotation, save=f"_annotation_final_level1_{annotation}.png")
        subset=adata[adata.obs['annotation_final_level1']==annotation]
        sc.pl.umap(subset, color="tech", save=f"_tech_{annotation}.png")

#adata=sc.read_h5ad("/home/bcd/revision_nature/integration_new_atlas/integration_combined/revision_atlas_cosmx_xenium_ver1_scvi.h5ad")
# Load pre-trained SCVI model
model = scvi.model.SCVI.load("/home/bcd/revision_nature/integration_new_atlas/integration_combined/models/scvi/model1", adata=adata)

# Define the model save path
model_path = "/home/bcd/revision_nature/integration_new_atlas/integration_combined/models/scanvi/model1"
history_plot_path = "/home/bcd/revision_nature/integration_new_atlas/integration_combined/models/scanvi/training_history_scanvi_1.png"  # Define plot save path

# Initialize SCANVI model
scanvi_model = scvi.model.SCANVI.from_scvi_model(model, labels_key="annotation_final_level1", unlabeled_category="Unknown")

# Train SCANVI model
scanvi_model.train(
    max_epochs=50, 
    plan_kwargs={"lr": 0.001}, 
    early_stopping=True, 
    early_stopping_patience=25,
    n_samples_per_label=100, 
    batch_size=512
)

# Save trained SCANVI model
scanvi_model.save(model_path, overwrite=True)

adata.obs["C_scANVI"] = scanvi_model.predict(adata)
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)

sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)

adata.write('/home/bcd/revision_nature/integration_new_atlas/integration_combined/revision_atlas_cosmx_xenium_ver1_scvi_scanvi.h5ad')

sc.settings.figdir = "/home/bcd/revision_nature/integration_new_atlas/integration_combined/plots_scanvi1"  # Update this path as needed

print(adata)
# Save each UMAP plot with specific file names
sc.pl.umap(adata, color="tech", save="_tech.png")
sc.pl.umap(adata, color="proj", save="_proj.png")
sc.pl.umap(adata, color="annotation_final_level1", save="_filter1_annotation_final_level1.png")
sc.pl.umap(adata, color="original_cosmx_annotation", save="_original_cosmx_annotation.png")


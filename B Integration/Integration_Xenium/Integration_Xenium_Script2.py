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

sc.settings.figdir = "/home/bcd/revision_nature/integration_new_atlas/integration_xenium/plots_scvi1"

adata = sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/integration_xenium/xenium_atlas_prescvi.h5ad')
yes_cells = adata[adata.obs['needs_mapping'] == 'Yes']
no_cells = adata[adata.obs['needs_mapping'] != 'Yes']

adata_reordered = yes_cells.concatenate(no_cells, join='inner', batch_key=None, index_unique=None)
adata = adata_reordered.copy()

print(adata)
print(adata.obs['tech'].value_counts())
print(adata.obs['proj'].value_counts())
print(adata.obs['needs_mapping'].value_counts())
print(adata.layers['counts'])

adata = adata.copy()
scvi.model.SCVI.setup_anndata(adata, layer="counts",
        batch_key="tech",
        categorical_covariate_keys=["proj","orig_ident"],
        continuous_covariate_keys=["nCount_RNA", "percent_mt"])

vae = scvi.model.SCVI(adata, n_layers=4, n_latent=30, gene_likelihood="nb", dropout_rate=0.1)
vae.train(max_epochs=1000, plan_kwargs={"lr": 0.001}, early_stopping=True, 
              batch_size=2048, early_stopping_patience=25)

model_name = f"/home/bcd/revision_nature/integration_new_atlas/integration_xenium/models/scvi/model1"
vae.save(model_name, overwrite=True)

latent = vae.get_latent_representation()
adata.obsm["X_scVI"] = latent
adata.layers["scvi_normalized"] = vae.get_normalized_expression(library_size=10e4)

sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)

adata.write("/home/bcd/revision_nature/integration_new_atlas/integration_xenium/revision_atlas_xenium_ver1_scvi.h5ad")

sc.pl.umap(adata, color="tech", save="_tech.png")
sc.pl.umap(adata, color="proj", save="_proj.png")
sc.pl.umap(adata, color="annotation_final_level1", save="_annotation_final_level1.png")
sc.pl.umap(adata, color="original_cosmx_annotation", save="_original_cosmx_annotation.png")

adata = sc.read_h5ad("/home/bcd/revision_nature/integration_new_atlas/integration_xenium/revision_atlas_xenium_ver1_scvi.h5ad")
Cosmx_cells_mask = (adata.obs["needs_mapping"] == "Yes")
umap_coords_Cosmx = adata.obsm["X_scVI"][Cosmx_cells_mask]
SN_cells_mask = (adata.obs["tech"] == "SN_SEQ")
umap_coords_non_Cosmx = adata.obsm["X_scVI"][SN_cells_mask]

nn = NearestNeighbors(n_neighbors=10, metric="euclidean")
nn.fit(umap_coords_non_Cosmx)
distances_all_to_non_Cosmx, indices_all_to_non_Cosmx = nn.kneighbors(umap_coords_Cosmx)

Cosmx_cell_ids = adata.obs_names[Cosmx_cells_mask]
sn_cell_ids = adata.obs_names[SN_cells_mask]

data = []
for j, cosmx_cell_id in enumerate(Cosmx_cell_ids):
    neighbor_indices = indices_all_to_non_Cosmx[j]
    neighbor_distances = distances_all_to_non_Cosmx[j]
    neighbor_sn_cell_ids = sn_cell_ids[neighbor_indices]
    rows = [cosmx_cell_id] + neighbor_sn_cell_ids.tolist() + neighbor_distances.tolist()
    data.append(rows)

df = pd.DataFrame(
    data,
    columns=["Cosmx_Cell_ID"] + [f"SN_Neighbor_{k}" for k in range(1, 11)] + [f"Distance_to_SN_{k}" for k in range(1, 11)]
)
df["Mean_Distance"] = df.iloc[:, -10:].mean(axis=1)

distance_dict = dict(zip(df["Cosmx_Cell_ID"], df["Mean_Distance"]))
adata.obs["mean_distance_10_scvi"] = adata.obs.index.map(distance_dict)

Cosmx_cells_mask = (adata.obs['needs_mapping'] == 'Yes')
snRNA_cells_mask = (adata.obs['tech'] == 'SN_SEQ')

CosMx_index = np.where(Cosmx_cells_mask)[0]
snRNA_index = np.where(snRNA_cells_mask)[0]

percentile_threshold = np.percentile(adata[CosMx_index, :].obs['mean_distance_10_scvi'], 10)
snRNA_cell_types = adata[snRNA_index, :].obs.annotation_final_level1.values

adata.obs.loc[adata.obs["needs_mapping"] == "Yes", "annotation_final_level1"] = "Unknown"
print(percentile_threshold)
for i in CosMx_index:
    distance = adata.obs.iloc[i, 13]
    if distance < percentile_threshold:
        neighbor_index = indices_all_to_non_Cosmx[i]
        unique, counts = np.unique((snRNA_cell_types[neighbor_index]), return_counts=True)
        adata.obs.iloc[i, 4] = unique[np.argmax(counts)]

print(adata.obs['annotation_final_level1'].value_counts())
adata.write('/home/bcd/revision_nature/integration_new_atlas/integration_xenium/revision_atlas_xenium_ver1_scvi.h5ad')

model = scvi.model.SCVI.load("/home/bcd/revision_nature/integration_new_atlas/integration_xenium/models/scvi/model1", adata=adata)

model_path = "/home/bcd/revision_nature/integration_new_atlas/integration_xenium/models/scanvi/model1"
history_plot_path = "/home/bcd/revision_nature/integration_new_atlas/integration_xenium/models/scanvi/training_history_scanvi_1.png"

scanvi_model = scvi.model.SCANVI.from_scvi_model(model, labels_key="annotation_final_level1", unlabeled_category="Unknown")

scanvi_model.train(
    max_epochs=50, 
    plan_kwargs={"lr": 0.001}, 
    early_stopping=True, 
    early_stopping_patience=25,
    n_samples_per_label=200, 
    batch_size=2048
)

scanvi_model.save(model_path, overwrite=True)

adata.obs["C_scANVI"] = scanvi_model.predict(adata)
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)

sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)

adata.write('/home/bcd/revision_nature/integration_new_atlas/integration_xenium/revision_atlas_xenium_scvi_scanvi.h5ad')

sc.settings.figdir = "/home/bcd/revision_nature/integration_new_atlas/integration_xenium/plots_scvi_scanvi1"

print(adata)
sc.pl.umap(adata, color="tech", save="_tech.png")
sc.pl.umap(adata, color="proj", save="_proj.png")
sc.pl.umap(adata, color="annotation_final_level1", save="_annotation_final_level1.png")
sc.pl.umap(adata, color="original_cosmx_annotation", save="_original_cosmx_annotation.png")

Cosmx_cells_mask = (adata.obs["needs_mapping"] == "Yes")
umap_coords_Cosmx = adata.obsm["X_scANVI"][Cosmx_cells_mask]
SN_cells_mask = (adata.obs["tech"] == "SN_SEQ")
umap_coords_non_Cosmx = adata.obsm["X_scANVI"][SN_cells_mask]

nn = NearestNeighbors(n_neighbors=10, metric="euclidean")
nn.fit(umap_coords_non_Cosmx)
distances_all_to_non_Cosmx, indices_all_to_non_Cosmx = nn.kneighbors(umap_coords_Cosmx)

Cosmx_cell_ids = adata.obs_names[Cosmx_cells_mask]
sn_cell_ids = adata.obs_names[SN_cells_mask]

data = []
for j, cosmx_cell_id in enumerate(Cosmx_cell_ids):
    neighbor_indices = indices_all_to_non_Cosmx[j]
    neighbor_distances = distances_all_to_non_Cosmx[j]
    neighbor_sn_cell_ids = sn_cell_ids[neighbor_indices]
    rows = [cosmx_cell_id] + neighbor_sn_cell_ids.tolist() + neighbor_distances.tolist()
    data.append(rows)

df = pd.DataFrame(
    data,
    columns=["Cosmx_Cell_ID"] + [f"SN_Neighbor_{k}" for k in range(1, 11)] + [f"Distance_to_SN_{k}" for k in range(1, 11)]
)
df["Mean_Distance"] = df.iloc[:, -10:].mean(axis=1)

distance_dict = dict(zip(df["Cosmx_Cell_ID"], df["Mean_Distance"]))
adata.obs["mean_distance_10_scvi"] = adata.obs.index.map(distance_dict)

Cosmx_cells_mask = (adata.obs['needs_mapping'] == 'Yes')
snRNA_cells_mask = (adata.obs['tech'] == 'SN_SEQ')

CosMx_index = np.where(Cosmx_cells_mask)[0]
snRNA_index = np.where(snRNA_cells_mask)[0]

snRNA_cell_types = adata[snRNA_index, :].obs.annotation_final_level1.values
percentile_threshold = np.percentile(adata[CosMx_index, :].obs['mean_distance_10_scvi'], 95)

adata.obs['mean_distance_10_scvi'] = adata.obs['mean_distance_10_scvi'].fillna(0)
adata_subset = adata[adata.obs['mean_distance_10_scvi'] < percentile_threshold].copy()
adata_subset.write('/home/bcd/revision_nature/integration_new_atlas/integration_xenium/revision_atlas_xenium_filtered.h5ad')

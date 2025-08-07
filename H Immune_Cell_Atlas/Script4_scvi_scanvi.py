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

# Load AnnData
adata = sc.read_h5ad('/home/bcd/revision_nature/immune_cell_atlas_xenium/xenium_cosmx_sn_prescvi.h5ad')
print(adata)
print(adata.obs['orig_ident'].value_counts())

# Filter to patients with at least 50 cells
patient_counts = adata.obs['orig_ident'].value_counts()
patients_to_keep = patient_counts[patient_counts >= 50].index
adata = adata[adata.obs['orig_ident'].isin(patients_to_keep)].copy()

# Setup SCVI
scvi.model.SCVI.setup_anndata(
    adata, 
    layer="counts",
    batch_key="tech",
    categorical_covariate_keys=["proj", "orig_ident"],
    continuous_covariate_keys=["nCount_RNA", "percent_mt"]
)

# Initialize and train SCVI
vae = scvi.model.SCVI(adata, n_layers=4, n_latent=30, gene_likelihood="nb", dropout_rate=0.1)
vae.train(max_epochs=1000, plan_kwargs={"lr": 0.001}, early_stopping=True, early_stopping_patience=25)

# Save model
model_name = "/home/bcd/revision_nature/immune_cell_atlas_xenium/models/scvi/model1_cosmx_xenium"
vae.save(model_name, overwrite=True)

# Get latent representation and normalized expression
latent = vae.get_latent_representation()
adata.obsm["X_scVI"] = latent
adata.layers["scvi_normalized"] = vae.get_normalized_expression(library_size=1e5)

# Compute neighbors and UMAP
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)

# Save final AnnData
adata.write("/home/bcd/revision_nature/immune_cell_atlas_xenium/immune_cells_xenium_cosmx_sn_scvi.h5ad")



adata=sc.read_h5ad('/home/bcd/revision_nature/immune_cell_atlas_xenium/immune_cells_xenium_cosmx_sn_scvi.h5ad')

model = scvi.model.SCVI.load("/home/bcd/revision_nature/immune_cell_atlas_xenium/models/scvi/model1_cosmx_xenium", adata=adata)
model_path = "/home/bcd/revision_nature/immune_cell_atlas_xenium/models/scanvi"

scanvi_model = scvi.model.SCANVI.from_scvi_model(model, labels_key="annotation_for_scanvi", unlabeled_category = "Unknown")

scanvi_model.train(max_epochs = 50, plan_kwargs={"lr":0.001}, early_stopping = True, early_stopping_patience = 25, n_samples_per_label=10)
scanvi_model.save(model_path, overwrite=True)

adata.obs["C_scANVI"] = scanvi_model.predict(adata)
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)

sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)

adata.write('/home/bcd/revision_nature/immune_cell_atlas_xenium/immune_cells_xenium_cosmx_sn_scvi_scanvi.h5ad')



adata=sc.read_h5ad('/home/bcd/revision_nature/immune_cell_atlas_xenium/immune_cells_xenium_cosmx_sn_scvi_scanvi.h5ad')

yes_cells = adata[adata.obs['tech'] == 'CosMx']
no_cells = adata[adata.obs['tech'] != 'CosMx']

# Concatenate the subsets
adata_reordered = yes_cells.concatenate(no_cells, join='inner', batch_key=None, index_unique=None)
adata = adata_reordered.copy()

Cosmx_cells_mask = (adata.obs["tech"] == "CosMx")
umap_coords_Cosmx = adata.obsm["X_scANVI"][Cosmx_cells_mask]
SN_cells_mask=(adata.obs["tech"] == "SN_SEQ")
umap_coords_non_Cosmx = adata.obsm["X_scANVI"][SN_cells_mask]

nn = NearestNeighbors(n_neighbors=10, metric="euclidean")
nn.fit(umap_coords_non_Cosmx)
distances_all_to_non_Cosmx, indices_all_to_non_Cosmx = nn.kneighbors(umap_coords_Cosmx)

Cosmx_cell_ids = adata.obs_names[Cosmx_cells_mask]
sn_cell_ids = adata.obs_names[SN_cells_mask]

# Create dataframe for distances
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

# Map mean distances back to AnnData
distance_dict = dict(zip(df["Cosmx_Cell_ID"], df["Mean_Distance"]))
adata.obs["mean_distance_10_scvi"] = adata.obs.index.map(distance_dict)

Cosmx_cells_mask = (adata.obs['tech'] == 'CosMx')
snRNA_cells_mask = (adata.obs['tech'] == 'SN_SEQ')

CosMx_index = np.where(Cosmx_cells_mask)[0]
snRNA_index = np.where(snRNA_cells_mask)[0]

snRNA_cell_types = adata[snRNA_index, :].obs.annotation_for_scanvi.values
percentile_threshold = np.percentile(adata[CosMx_index, :].obs['mean_distance_10_scvi'], 50)

adata.obs.loc[adata.obs["tech"] == "CosMx", "annotation_for_scanvi"] = "Unknown"
print(adata)
print(percentile_threshold)
for i in CosMx_index:
            distance = adata.obs.iloc[i, 20]  #column 20 of adata.obs should be distance
            if distance < percentile_threshold:
                neighbor_index = indices_all_to_non_Cosmx[i]
                unique, counts = np.unique((snRNA_cell_types[neighbor_index]), return_counts=True)
                adata.obs.iloc[i, 19] = unique[np.argmax(counts)]


model = scvi.model.SCVI.load("/home/bcd/revision_nature/immune_cell_atlas_xenium/models/scvi/model1_cosmx_xenium", adata=adata)
model_path = "/home/bcd/revision_nature/immune_cell_atlas_xenium/models/scanvi"

scanvi_model = scvi.model.SCANVI.from_scvi_model(model, labels_key="annotation_for_scanvi", unlabeled_category = "Unknown")

scanvi_model.train(max_epochs = 100, plan_kwargs={"lr":0.001}, early_stopping = True, early_stopping_patience = 25, n_samples_per_label=10)
scanvi_model.save(model_path, overwrite=True)

adata.obs["C_scANVI"] = scanvi_model.predict(adata)
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)

sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)

adata.write('/home/bcd/revision_nature/immune_cell_atlas_xenium/immune_cells_xenium_cosmx_sn_scvi_scanvi2.h5ad')

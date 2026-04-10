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

adata=sc.read_h5ad('/home/bcd/revision_nature/B_Cell_Atlas/b_cell_atlas_processed.h5ad')
adata.obs.index=adata.obs['orig_ident'].astype(str)+'_'+adata.obs.index
print(adata.var)
print(adata)
df=pd.read_csv('/home/bcd/revision_nature/B_Cell_Atlas/genes_black_WYC.csv')

genes_to_exclude = df['Genes'].drop_duplicates().tolist()  # Assuming the genes are in the first column

# Subset adata to exclude these genes
adata_subset = adata[:, ~adata.var_names.isin(genes_to_exclude)].copy()
print(adata_subset)

sc.pp.normalize_total(adata_subset, inplace=True)
sc.pp.log1p(adata_subset)
sc.pp.highly_variable_genes(adata_subset, n_top_genes=3000, batch_key='dataid', flavor='seurat')

adata_hvg = adata_subset[:, adata_subset.var['highly_variable']].copy()
print(adata_hvg)
sc.pp.filter_cells(adata_hvg, min_counts=200)
adata_hvg.write('/home/bcd/revision_nature/B_Cell_Atlas/b_atlas_3000hvg_prescvi.h5ad')

##keeping only samples with enough cells to remove more batch effect
adata=adata[adata.obs['dataid']!='OA020']
adata=adata[adata.obs['dataid']!='IH002']
adata=adata[adata.obs['dataid']!='OA044']
adata=adata[adata.obs['dataid']!='OA013']
adata=adata[adata.obs['dataid']!='OA025']
adata=adata[adata.obs['dataid']!='OA004']
adata=adata[adata.obs['dataid']!='OA011']
adata=adata[adata.obs['dataid']!='OA070']
adata=adata[adata.obs['dataid']!='OA047']
adata=adata[adata.obs['dataid']!='OA069']
adata=adata[adata.obs['dataid']!='OA028']
adata=adata[adata.obs['dataid']!='OA041']
# Count the number of cells per patient
patient_counts = adata.obs['patient'].value_counts()
# Identify patients with at least 500 cells
patients_to_keep = patient_counts[patient_counts >= 500].index
# Filter the AnnData object to keep only those patients
adata = adata[adata.obs['patient'].isin(patients_to_keep)].copy()

# Print the updated counts
print(adata.obs['dataid'].value_counts())
print(adata.obs['patient'].value_counts())

adata.obs["celltype"].value_counts()
cell_identities = {'B.09.DUSP4+AtM': 'AtM', 'B.08.ITGB1+SwBm': 'SwBm', 'B.07.CCR7+ACB3': 'ACB', 'B.06.NR4A2+ACB2': 'ACB', 
                   'B.10.ENO1+Pre_GCB': 'Unknown', 'B.03.HSP+B': 'Unknown', 'B.14.Plasmablast':'Plasma','B.15.Plasma cell':'Plasma',
                   'B.02.IFIT3+B':'', 'B.01.TCL1A+naiveB': 'naiveB', 
                   'B.04.MT1X+B': 'Unknown', 'B.05.EGR1+ACB': 'ACB', 'B.11.SUGCT+DZ_GCB': 'GCB', 'B.12.LMO2+LZ_GCB': 'GCB', 'B.13.Cycling_GCB': 'GCB'}
adata.obs["celltype_coarse"] = adata.obs['celltype'].map(cell_identities).astype('category')
sc.pl.umap(adata, color = "celltype_coarse", save='_celltype_coarse.png')

adata.write('/home/bcd/revision_nature/B_Cell_Atlas/b_atlas_3000hvg_filter1.h5ad')

print(adata)
adata=adata.copy()

scvi.model.SCVI.setup_anndata(adata, 
        layer="counts",
        batch_key="dataid",
        categorical_covariate_keys=["orig_ident"],
        continuous_covariate_keys=["nCount_RNA", "percent_mt"])

# Initialize and train the SCVI model
vae = scvi.model.SCVI(adata, n_layers=3, n_latent=30, gene_likelihood="nb", dropout_rate=0.1)
vae.train(max_epochs=1000, plan_kwargs={"lr": 0.001}, early_stopping=True, 
              batch_size=2048, early_stopping_patience=25)

# Save the model with a unique name for each HVG count
model_name = f"/home/bcd/revision_nature/B_Cell_Atlas/models/scvi/model2"
vae.save(model_name, overwrite=True)

# Get latent representation and normalized expression
latent = vae.get_latent_representation()
adata.obsm["X_scVI"] = latent
adata.layers["scvi_normalized"] = vae.get_normalized_expression(library_size=10e4)

# Compute neighbors and UMAP
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)

adata.write("/home/bcd/revision_nature/B_Cell_Atlas/b_atlas_3000hvg_filter1_scvi.h5ad")
sc.settings.figdir = "/home/bcd/revision_nature/B_Cell_Atlas/plots_scvi2" 


model = scvi.model.SCVI.load("/home/bcd/revision_nature/B_Cell_Atlas/models/scvi/model2", adata=adata)
model_path = "/home/bcd/revision_nature/B_Cell_Atlas/models/scanvi"

scanvi_model = scvi.model.SCANVI.from_scvi_model(model, labels_key="celltype_coarse", unlabeled_category = "Unknown")

scanvi_model.train(max_epochs = 100, plan_kwargs={"lr":0.001}, early_stopping = True, early_stopping_patience = 25, batch_size=2048)
scanvi_model.save(model_path, overwrite=True)

adata.obs["C_scANVI"] = scanvi_model.predict(adata)
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)

sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)

adata.write('/home/bcd/revision_nature/B_Cell_Atlas/b_atlas_3000hvg_filter1_scvi_scanvi.h5ad')

model = scvi.model.SCVI.load("/home/bcd/revision_nature/B_Cell_Atlas/models/scvi/model2", adata=adata)
model_path = "/home/bcd/revision_nature/B_Cell_Atlas/models/scanvi"

scanvi_model = scvi.model.SCANVI.from_scvi_model(model, labels_key="celltype_coarse", unlabeled_category = "Unknown")

scanvi_model.train(max_epochs = 100, plan_kwargs={"lr":0.001}, early_stopping = True, early_stopping_patience = 25, n_samples_per_label=50, batch_size=256)
scanvi_model.save(model_path, overwrite=True)

adata.obs["C_scANVI"] = scanvi_model.predict(adata)
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)

sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)

adata.write('/home/bcd/revision_nature/B_Cell_Atlas/b_atlas_3000hvg_filter1_scvi_scanvi_v2.h5ad')

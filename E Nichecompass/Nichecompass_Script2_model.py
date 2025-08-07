import os
import random
import warnings
from datetime import datetime
import json
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
import squidpy as sq
from matplotlib import gridspec
from sklearn.preprocessing import MinMaxScaler
from nichecompass.models import NicheCompass
from nichecompass.utils import (add_gps_from_gp_dict_to_adata,
                                create_new_color_dict,
                                compute_communication_gp_network,
                                visualize_communication_gp_network,
                                extract_gp_dict_from_mebocost_ms_interactions,
                                extract_gp_dict_from_nichenet_lrt_interactions,
                                extract_gp_dict_from_omnipath_lr_interactions,
                                filter_and_combine_gp_dict_gps_v2,
                                generate_enriched_gp_info_plots)

os.environ["CUDA_VISIBLE_DEVICES"] = "0"
sc.settings.figdir = "/home/bcd/revision_nature/nichecompass/revision_atlas/plots_script1"
adata=sc.read_h5ad('/home/bcd/revision_nature/nichecompass/revision_atlas/prenichecompass_adata.h5ad')

print(adata)
adata.obsp['spatial_connectivities']=adata.obsp['20_micron_connectivities'].copy() ##for some reason there is a bug 
print(adata)
species = "human"
batches = [ "1061_CosMx", "1062_CosMx", "1063_CosMx", "1064_CosMx", "1068_CosMx", "1070_CosMx", "1071_CosMx", "1072_CosMx",  "HK2695_CosMx", 
 "HK2753_CosMx", "HK2841_CosMx", "HK2844_CosMx", "HK2844_2_CosMx", "HK2873_CosMx", "HK2874_CosMx", "HK2924_CosMx", "HK2989_CosMx", "HK2990_CosMx", 
"HK3035_CosMx", "HK3035_2_CosMx", "HK3039_CosMx", "HK3063_CosMx", "HK3066_CosMx", "HK3068_CosMx", "HK3069_CosMx", "HK3070_CosMx",  "HK3106_CosMx", 
"HK3421_CosMx", "HK3469_CosMx", "HK3474_CosMx", "HK3531_CosMx", "HK3531_2_CosMx", "HK3535_CosMx", "HK3542_CosMx", "HK3565_CosMx", "HK3588_CosMx", "HK3591_CosMx", 
"HK3594_CosMx", "HK3606_CosMx", "HK3612_CosMx", "HK3614_CosMx", "HK3616_CosMx", "HK3623_CosMx", "HK3624_CosMx",  "HK3626_CosMx", "HK3631_CosMx", 
"HK3647_CosMx", "Pediatric1_CosMx"]

### Model ###
# AnnData keys
counts_key = "counts"
adj_key = "spatial_connectivities" #20_micron_connectivities

cat_covariates_keys = ["unique_sample_identifier"]
original_niche_key= "niches_annotation_based"
gp_names_key = "nichecompass_gp_names"
active_gp_names_key = "nichecompass_active_gp_names"
gp_targets_mask_key = "nichecompass_gp_targets"
gp_targets_categories_mask_key = "nichecompass_gp_targets_categories"
gp_sources_mask_key = "nichecompass_gp_sources"
gp_sources_categories_mask_key = "nichecompass_gp_sources_categories"
latent_key = "nichecompass_latent"

# Architecture
cat_covariates_embeds_injection = ["gene_expr_decoder"]
cat_covariates_embeds_nums = [48]
cat_covariates_no_edges = [True]
conv_layer_encoder = "gatv2conv" # change to "gatv2conv" if enough compute and memory
active_gp_thresh_ratio = 0.01

# Trainer
n_epochs = 1000
n_epochs_all_gps = 25
lr = 0.001
lambda_edge_recon = 500000.
lambda_gene_expr_recon = 300.
lambda_l1_masked = 0. # prior GP  regularization
lambda_l1_addon = 30. # de novo GP regularization
edge_batch_size = 12000 # increase if more memory available or decrease to save memory
n_sampled_neighbors = 10
use_cuda_if_available = True

### Analysis ###
cell_type_key = "annotation_updated"
latent_leiden_resolution = 0.6
latent_cluster_key = f"latent_leiden_{str(latent_leiden_resolution)}"
sample_key = "unique_sample_identifier"
spot_size = 0.2
differential_gp_test_results_key = "nichecompass_differential_gp_test_results"
warnings.filterwarnings("ignore")
##

# Define paths
ga_data_folder_path = "/home/bcd/revision_nature/nichecompass/data/gene_annotations"
gp_data_folder_path = "/home/bcd/revision_nature/nichecompass/data/gene_programs"
omnipath_lr_network_file_path = f"{gp_data_folder_path}/omnipath_lr_network.csv"
collectri_tf_network_file_path = f"{gp_data_folder_path}/collectri_tf_network_{species}.csv"
nichenet_lr_network_file_path = f"{gp_data_folder_path}/nichenet_lr_network_v2_{species}.csv"
nichenet_ligand_target_matrix_file_path = f"{gp_data_folder_path}/nichenet_ligand_target_matrix_v2_{species}.csv"
mebocost_enzyme_sensor_interactions_folder_path = f"{gp_data_folder_path}/metabolite_enzyme_sensor_gps"
gene_orthologs_mapping_file_path = f"{ga_data_folder_path}/human_mouse_gene_orthologs.csv"

artifacts_folder_path = f"/home/bcd/revision_nature/nichecompass/data/artifacts"
model_folder_path = "/home/bcd/revision_nature/nichecompass/revision_atlas/model"
figure_folder_path = "/home/bcd/revision_nature/nichecompass/revision_atlas/plots_script1"

# Initialize model
model = NicheCompass(adata,
                     counts_key=counts_key,
                     adj_key=adj_key,
                     cat_covariates_embeds_injection=cat_covariates_embeds_injection,
                     cat_covariates_keys=cat_covariates_keys,
                     cat_covariates_no_edges=cat_covariates_no_edges,
                     cat_covariates_embeds_nums=cat_covariates_embeds_nums,
                     gp_names_key=gp_names_key,
                     active_gp_names_key=active_gp_names_key,
                     gp_targets_mask_key=gp_targets_mask_key,
                     gp_targets_categories_mask_key=gp_targets_categories_mask_key,
                     gp_sources_mask_key=gp_sources_mask_key,
                     gp_sources_categories_mask_key=gp_sources_categories_mask_key,
                     latent_key=latent_key,
                     conv_layer_encoder=conv_layer_encoder,
                     active_gp_thresh_ratio=active_gp_thresh_ratio)

# Train model
model.train(n_epochs=n_epochs,
            n_epochs_all_gps=n_epochs_all_gps,
            lr=lr,
            lambda_edge_recon=lambda_edge_recon,
            lambda_gene_expr_recon=lambda_gene_expr_recon,
            lambda_l1_masked=lambda_l1_masked,
            edge_batch_size=edge_batch_size,
            n_sampled_neighbors=n_sampled_neighbors,
            use_cuda_if_available=use_cuda_if_available,
            verbose=False)

# Save trained model
model.save(dir_path=model_folder_path,
           overwrite=True,
           save_adata=True,
           adata_file_name="postnichecompass_adata.h5ad")

# Compute latent neighbor graph
sc.pp.neighbors(model.adata,
                use_rep=latent_key,
                key_added=latent_key)

# Compute UMAP embedding
sc.tl.umap(model.adata,
           neighbors_key=latent_key)

# Save trained model
model.save(dir_path=model_folder_path,
           overwrite=True,
           save_adata=True,
           adata_file_name="postnichecompass_adata.h5ad")

# Load trained model
model = NicheCompass.load(dir_path=model_folder_path,
                          adata=None,
                          adata_file_name="postnichecompass_adata.h5ad",
                          gp_names_key=gp_names_key)

samples = model.adata.obs[sample_key].unique().tolist()

batch_colors = create_new_color_dict(
    adata=model.adata,
    cat_key=cat_covariates_keys[0])

cell_type_colors = create_new_color_dict(
    adata=model.adata,
    cat_key=cell_type_key)

sc.tl.leiden(adata=model.adata,
             resolution=latent_leiden_resolution,
             key_added=latent_cluster_key,
             neighbors_key=latent_key)

latent_cluster_colors = create_new_color_dict(
    adata=model.adata,
    cat_key=latent_cluster_key)

sc.pl.umap(adata=model.adata,
           color=[latent_cluster_key],
           palette=latent_cluster_colors,
           title=f"Niches in Latent Space",
           save='_niches_in_latent_space_0_6.png')

sc.pl.umap(adata=model.adata,
           color=[original_niche_key],
           title=f"Original niche annotation in Latent Space",
           save='_original_niches_in_latent_space_0_6.png')

adata.write('/home/bcd/revision_nature/nichecompass/revision_atlas/postnichecompass_adata2.h5ad')

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
adata=sc.read_h5ad('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad') #update


print(adata)
print(adata.obs['unique_sample_identifier'].value_counts())
adata=adata[adata.obs['niches_annotation_based']!='Unknown'].copy() #update
adata=adata[adata.obs['tech']=='CosMx'].copy()
adata.X=adata.layers['counts'].copy()
print(adata)

cosmx_raw=sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/integration_cosmx/all_cosmx_all_genes_revision_raw.h5ad')
##subet adata to only include CosMx genes
common_genes = adata.var_names.intersection(cosmx_raw.var_names)
# Subset adata to only keep those genes
adata = adata[:, common_genes].copy()
print(adata)
print(adata.obs['unique_sample_identifier'].value_counts())
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
adj_key = "20_micron_connectivities"
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
edge_batch_size = 2048 # increase if more memory available or decrease to save memory
n_sampled_neighbors = 10
use_cuda_if_available = True

### Analysis ###
cell_type_key = "annotation_updated"
latent_leiden_resolution = 0.2
latent_cluster_key = f"latent_leiden_{str(latent_leiden_resolution)}"
sample_key = "unique_sample_identifier"
spot_size = 0.2
differential_gp_test_results_key = "nichecompass_differential_gp_test_results"
warnings.filterwarnings("ignore")
##
# Get time of notebook execution for timestamping saved artifacts
now = datetime.now()
current_timestamp = now.strftime("%d%m%Y_%H%M%S")

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
##
os.makedirs(model_folder_path, exist_ok=True)
os.makedirs(figure_folder_path, exist_ok=True)

with open('/home/bcd/revision_nature/nichecompass/combined_gp_dict.json', 'r') as file:
    loaded_gp_dict = json.load(file)

#print(loaded_gp_dict)
combined_gp_dict=loaded_gp_dict.copy()
#print(combined_gp_dict)

print(f"Number of gene programs after filtering and combining: "
      f"{len(combined_gp_dict)}.")


adata_batch_list = []

for batch in batches:
    print(f"Processing batch {batch}...")
    print("Loading data...")
    adata_batch = adata[adata.obs['unique_sample_identifier']==batch]
    print("Computing spatial neighborhood graph...\n") 
    # Make adjacency matrix symmetric
    adata_batch.obsp[adj_key] = (
        adata_batch.obsp[adj_key].maximum(
            adata_batch.obsp[adj_key].T))
    adata_batch_list.append(adata_batch)
adata = ad.concat(adata_batch_list, join="inner")
print(adata)

# Combine spatial neighborhood graphs as disconnected components
batch_connectivities = []
len_before_batch = 0
for i in range(len(adata_batch_list)):
    if i == 0: # first batch
        after_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[0].shape[0],
            (adata.shape[0] -
            adata_batch_list[0].shape[0])))
        batch_connectivities.append(sp.hstack(
            (adata_batch_list[0].obsp[adj_key],
            after_batch_connectivities_extension)))
    elif i == (len(adata_batch_list) - 1): # last batch
        before_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0],
            (adata.shape[0] -
            adata_batch_list[i].shape[0])))
        batch_connectivities.append(sp.hstack(
            (before_batch_connectivities_extension,
            adata_batch_list[i].obsp[adj_key])))
    else: # middle batches
        before_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0], len_before_batch))
        after_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0],
            (adata.shape[0] -
            adata_batch_list[i].shape[0] -
            len_before_batch)))
        batch_connectivities.append(sp.hstack(
            (before_batch_connectivities_extension,
            adata_batch_list[i].obsp[adj_key],
            after_batch_connectivities_extension)))
    len_before_batch += adata_batch_list[i].shape[0]
adata.obsp[adj_key] = sp.vstack(batch_connectivities)
print(adata)
add_gps_from_gp_dict_to_adata(
    gp_dict=combined_gp_dict,
    adata=adata,
    gp_targets_mask_key=gp_targets_mask_key,
    gp_targets_categories_mask_key=gp_targets_categories_mask_key,
    gp_sources_mask_key=gp_sources_mask_key,
    gp_sources_categories_mask_key=gp_sources_categories_mask_key,
    gp_names_key=gp_names_key,
    min_genes_per_gp=2,
    min_source_genes_per_gp=1,
    min_target_genes_per_gp=1,
    max_genes_per_gp=None,
    max_source_genes_per_gp=None,
    max_target_genes_per_gp=None)

cell_type_colors = create_new_color_dict(
    adata=adata,
    cat_key=cell_type_key)

samples = adata.obs[sample_key].unique().tolist()

for sample in samples:
    adata_batch = adata[adata.obs[sample_key] == sample]
    
    print(f"Summary of sample {sample}:")
    print(f"Number of nodes (observations): {adata_batch.layers[counts_key].shape[0]}")
    print(f"Number of node features (genes): {adata_batch.layers[counts_key].shape[1]}")

    # Visualize cell-level annotated data in physical space
    sc.pl.spatial(adata_batch,color=cell_type_key,palette=cell_type_colors,spot_size=spot_size,save='_cell_types.png')  

adata.write('/home/bcd/revision_nature/nichecompass/revision_atlas/prenichecompass_adata.h5ad')


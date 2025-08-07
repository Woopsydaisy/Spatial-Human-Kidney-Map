import csv
import anndata as ad
import gzip
import os
import scipy.io
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gc
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

adata = sc.read_h5ad('/home/bcd/revision_nature/imputation_new_atlas/xenium_cosmx_sn_imputed_neighborhood.h5ad')

adata.obs['niches_annotation_based'].value_counts()

# Initialize a mask with all False
keep_mask = pd.Series(False, index=adata.obs_names)

# Group by niche
for niche, niche_df in adata.obs.groupby("niches_annotation_based"):
    total_cells = len(niche_df)
    counts = niche_df['annotation_updated'].value_counts()
    # Default: keep annotations >1% within this niche
    valid_annotations = counts[counts > 0.01 * total_cells].index.tolist()
    # Exception: always keep EC_Lymph in the Immune niche
    if niche == "Immune niche":
        if "EC_Lymph" not in valid_annotations:
            valid_annotations.append("EC_Lymph")
    
    # Mark those cells to keep
    niche_keep = niche_df[niche_df['annotation_updated'].isin(valid_annotations)].index
    keep_mask[niche_keep] = True

# Subset AnnData object
adata_filtered = adata[keep_mask]

print(adata_filtered)

adata_filtered.obs['microenvironment'] = adata_filtered.obs['niches_annotation_based']
adata_filtered.obs['microenvironment_cell_types'] = (
    adata_filtered.obs['microenvironment'].astype(str) + '_' + adata_filtered.obs['annotation_updated'].astype(str)
)

counts = adata_filtered.obs['microenvironment'].value_counts()
niches_to_keep = counts[counts >= 50].index

# Create a boolean mask
mask = adata_filtered.obs['microenvironment'].isin(niches_to_keep)
# Assign a filtered view for downstream use
adata_filtered_view = adata_filtered[mask]  # This is a view, not a full copy

counts = adata_filtered_view.obs['microenvironment_cell_types'].value_counts()
# Identify microenvironments with at least 50 cells
niches_to_keep = counts[counts >= 50].index
mask = adata_filtered_view.obs['microenvironment_cell_types'].isin(niches_to_keep)
adata_filtered_view2 = adata_filtered_view[mask]

print(adata_filtered_view2.obs['microenvironment_cell_types'].value_counts())

metaData = adata_filtered_view2.obs[["microenvironment_cell_types"]]
metaData["barcode_sample"] = adata_filtered_view2.obs.index
metaData = metaData.rename(columns={'microenvironment_cell_types': 'cell_type'})
metaData = metaData[['barcode_sample', 'cell_type']]
print(metaData)

metaData.to_csv("/home/bcd/revision_nature/neighbors_neighborhood/cellphonedb/metadata_niches.csv", index=False)

microenvironment = adata_filtered_view2.obs[["microenvironment"]]
microenvironment["microenvironment_cell_types"] = adata_filtered_view2.obs["microenvironment_cell_types"]
microenvironment = microenvironment.rename(columns={'microenvironment_cell_types': 'cell_type'})
microenvironment = microenvironment[['cell_type', 'microenvironment']]

microenvironment.to_csv("/home/bcd/revision_nature/neighbors_neighborhood/cellphonedb/microenvironment_niches.csv", index=False)
adata_filtered_view2.write('/home/bcd/revision_nature/neighbors_neighborhood/cellphonedb/niches_for_cellphone.h5ad')

counts_file_path = '/home/bcd/revision_nature/neighbors_neighborhood/cellphonedb/niches_for_cellphone.h5ad'
meta_file_path = "/home/bcd/revision_nature/neighbors_neighborhood/cellphonedb/metadata_niches.csv"
metadata = pd.read_csv(meta_file_path, sep = ',')
microenvs_file_path = "/home/bcd/revision_nature/neighbors_neighborhood/cellphonedb/microenvironment_niches.csv"
microenvironment = pd.read_csv(microenvs_file_path, sep = ',')
out_path = '/home/bcd/revision_nature/neighbors_neighborhood/cellphonedb/'
cpdb_file_path = '/home/bcd/revision_nature/neighbors_neighborhood/cellphonedb/cellphonedb.zip'


cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
    microenvs_file_path = microenvs_file_path,       # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                       # optional: whether to score interactions or not.
    iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
    threshold = 0.01,                             # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 4,                                     # number of threads to use in the analysis.
    debug_seed = 42,                                 # debug randome seed. To disable >=0.
    result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    #subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
    #subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = out_path,                          # Path to save results.
    output_suffix = None                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )

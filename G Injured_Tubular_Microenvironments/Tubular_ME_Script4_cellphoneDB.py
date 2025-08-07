import csv
import anndata as ad
import gzip
import os
import scipy.io
import gc
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import leidenalg as la
from pathlib import Path
import umap
import seaborn as sns
from scipy.sparse import csr_matrix
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# Load the AnnData object
adata = sc.read_h5ad('/home/bcd/revision_nature/iPT_iLOH_ME/iPT_iLOH_ME_20um_imputed_object.h5ad')
print(adata.obs['iPT_iLOH_ME_20um'].value_counts())

# Initialize a keep mask
keep_mask = pd.Series(False, index=adata.obs_names)

# Filter annotations that make up more than 1% of each microenvironment
for niche, niche_df in adata.obs.groupby("iPT_iLOH_ME_20um"):
    total_cells = len(niche_df)
    counts = niche_df['annotation_updated'].value_counts()
    valid_annotations = counts[counts > 0.01 * total_cells].index.tolist()
    niche_keep = niche_df[niche_df['annotation_updated'].isin(valid_annotations)].index
    keep_mask[niche_keep] = True

# Subset the data
adata_filtered = adata[keep_mask]
print(adata_filtered.obs['iPT_iLOH_ME_20um'].value_counts())

# Clean up memory
del adata
gc.collect()

# Create metadata columns for downstream analysis
adata_filtered.obs['microenvironment'] = (
    adata_filtered.obs['iPT_iLOH_ME_20um'].astype(str) + '_' +
    adata_filtered.obs['unique_sample_identifier'].astype(str)
)
adata_filtered.obs['microenvironment_cell_types'] = (
    adata_filtered.obs['microenvironment'].astype(str) + '_' +
    adata_filtered.obs['annotation_updated'].astype(str)
)

# Save filtered AnnData object
adata_filtered.write('/home/bcd/revision_nature/iPT_iLOH_ME/iPT_iLOH_ME_20um_imputed_object_subset.h5ad')
print(adata_filtered.obs['microenvironment_cell_types'].value_counts())

# Prepare metadata for CellPhoneDB
metaData = adata_filtered.obs[["microenvironment_cell_types"]].copy()
metaData["barcode_sample"] = adata_filtered.obs.index
metaData = metaData.rename(columns={'microenvironment_cell_types': 'cell_type'})
metaData = metaData[['barcode_sample', 'cell_type']]
metaData.to_csv("/home/bcd/revision_nature/iPT_iLOH_ME/cellphonedb_output/metadata_niches.csv", index=False)

# Create and save microenvironment-cell type mapping
microenvironment = adata_filtered.obs[["microenvironment"]].copy()
microenvironment["microenvironment_cell_types"] = adata_filtered.obs["microenvironment_cell_types"]
microenvironment = microenvironment.rename(columns={'microenvironment_cell_types': 'cell_type'})
microenvironment = microenvironment[['cell_type', 'microenvironment']]
microenvironment.to_csv("/home/bcd/revision_nature/iPT_iLOH_ME/cellphonedb_output/microenvironment_niches.csv", index=False)

counts_file_path = '/home/bcd/revision_nature/iPT_iLOH_ME/iPT_iLOH_ME_20um_imputed_object_subset.h5ad'
meta_file_path = "/home/bcd/revision_nature/iPT_iLOH_ME/cellphonedb_output/metadata_niches.csv"
metadata = pd.read_csv(meta_file_path, sep = ',')
microenvs_file_path = "/home/bcd/revision_nature/iPT_iLOH_ME/cellphonedb_output/microenvironment_niches.csv"
microenvironment = pd.read_csv(microenvs_file_path, sep = ',')
out_path = '/home/bcd/revision_nature/iPT_iLOH_ME/cellphonedb_output/'
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
    threads = 8,                                     # number of threads to use in the analysis.
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


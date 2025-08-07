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

adata=sc.read_h5ad('/home/bcd/revision_nature/immune_ME/immune_ME_20um_object_for_cellphone.h5ad')
print(adata)
print(adata.obs['Immune_ME_20um'].value_counts())
print(adata.obs['annotation_updated'].value_counts())
print(adata.obs['immune_cell_annotation_combined'].value_counts())

adata.obs['annotation_updated']=adata.obs['annotation_updated'].fillna('Immune')
print(adata.obs['annotation_updated'].value_counts())

##one annotation
adata.obs['immune_annotation_with_annotation_updated'] = np.where(
    adata.obs['annotation_updated'] != 'Immune',
    adata.obs['annotation_updated'],
    adata.obs['immune_cell_annotation_combined']
)
print(adata.obs['immune_annotation_with_annotation_updated'].value_counts())

adata=adata[adata.obs['immune_annotation_with_annotation_updated']!='Unknown'].copy()
adata=adata[adata.obs['immune_annotation_with_annotation_updated']!='mDC'].copy()

adata.obs['microenvironment']=adata.obs['Immune_ME_20um'].astype(str)+'_'+adata.obs['unique_sample_identifier'].astype(str)
print(adata.obs['microenvironment'].value_counts())
cell_type_counts = adata.obs['microenvironment'].value_counts()
valid_types = cell_type_counts[cell_type_counts >= 50].index
adata = adata[adata.obs['microenvironment'].isin(valid_types)].copy()
print(adata.obs['microenvironment'].value_counts())

adata.obs['microenvironment_cell_types']=adata.obs['microenvironment'].astype(str)+'_'+adata.obs['immune_annotation_with_annotation_updated'].astype(str)
adata.write('/home/bcd/revision_nature/immune_ME/immune_ME_20um_object_for_cellphone_subset.h5ad')
print(adata.obs['microenvironment_cell_types'].value_counts())
cell_type_counts = adata.obs['microenvironment_cell_types'].value_counts()
valid_types = cell_type_counts[cell_type_counts >= 10].index
adata = adata[adata.obs['microenvironment_cell_types'].isin(valid_types)].copy()
print(adata.obs['microenvironment_cell_types'].value_counts())

metaData = adata.obs[["microenvironment_cell_types"]]
metaData["barcode_sample"] = adata.obs.index
metaData = metaData.rename(columns={'microenvironment_cell_types': 'cell_type'})
metaData = metaData[['barcode_sample', 'cell_type']]
print(metaData)

metaData.to_csv("/home/bcd/revision_nature/immune_ME/cellphonedb_output/metadata_niches.csv", index = False)

microenvironment=adata.obs[["microenvironment"]]
microenvironment["microenvironment_cell_types"] = adata.obs["microenvironment_cell_types"]
microenvironment = microenvironment.rename(columns={'microenvironment_cell_types': 'cell_type'})
microenvironment = microenvironment[['cell_type','microenvironment']]
print(microenvironment)
microenvironment.to_csv("/home/bcd/revision_nature/immune_ME/cellphonedb_output/microenvironment_niches.csv", index = False)

counts_file_path = '/home/bcd/revision_nature/immune_ME/immune_ME_20um_object_for_cellphone_subset.h5ad'
meta_file_path = "/home/bcd/revision_nature/immune_ME/cellphonedb_output/metadata_niches.csv"
metadata = pd.read_csv(meta_file_path, sep = ',')
microenvs_file_path = "/home/bcd/revision_nature/immune_ME/cellphonedb_output/microenvironment_niches.csv"
microenvironment = pd.read_csv(microenvs_file_path, sep = ',')
out_path = '/home/bcd/revision_nature/immune_ME/cellphonedb_output/'
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
    threads = 32,                                     # number of threads to use in the analysis.
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

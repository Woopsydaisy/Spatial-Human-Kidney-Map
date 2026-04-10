from anndata import AnnData
import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.sparse as sp
import scipy.io
import scipy.sparse as sp

# Define input file paths
counts_file = "/home/bcd/revision_nature/B_Cell_Atlas/expression_matrix_b_cells.mtx"
metadata_file = "/home/bcd/revision_nature/B_Cell_Atlas/metadata_b_cells.csv"
genes_file = "/home/bcd/revision_nature/B_Cell_Atlas/gene_names_b_cell.csv"
# Load sparse counts matrix from MTX file
print("Loading counts matrix from MTX file...")
counts_sparse = scipy.io.mmread(counts_file).tocsr()

print("Counts matrix loaded.")
print(f"Counts matrix shape: {counts_sparse.shape}")

# Load metadata
print("Loading metadata...")
metadata = pd.read_csv(metadata_file, index_col=0)

print("Metadata loaded.")
print(metadata.head(5))
print(f"Metadata shape: {metadata.shape}")

# Load gene names
print("Loading gene names...")
gene_names = pd.read_csv(genes_file)
print("Gene names loaded.")
print(gene_names.head(5))
print(gene_names.shape)
assert counts_sparse.shape[0] == len(gene_names), "Number of genes does not match rows in counts matrix!"
# Ensure indices match between counts and metadata
assert counts_sparse.shape[1] == metadata.shape[0], "Counts columns and metadata rows do not match!"

# Transpose counts matrix for AnnData (cells as rows, genes as columns)
counts_sparse = counts_sparse.T

print("Creating AnnData object...")
# Create AnnData object
adata = sc.AnnData(
    X=counts_sparse,              # Expression matrix (sparse)
    obs=metadata,                 # Cell metadata
    var=pd.DataFrame(index=gene_names.iloc[:, 0])  # Gene metadata
)

# Save the AnnData object
output_file = "/home/bcd/revision_nature/B_Cell_Atlas/b_cell_atlas.h5ad"
adata.write(output_file)
print(f"AnnData object saved to: {output_file}")


adata.layers['counts']=adata.X.copy()

adata.obs['orig_ident']=adata.obs['orig.ident'].copy()
adata.obs['orig_ident']=adata.obs['orig_ident'].astype('category')
del adata.obs['orig.ident']

adata.obs['percent_mt']=adata.obs['percent.mt'].copy()
adata.obs['percent_mt']=adata.obs['percent_mt'].astype(float)
del adata.obs['percent.mt']

adata.obs['nCount_RNA']=adata.obs['nCount_RNA'].astype(int)
adata.obs['nFeature_RNA']=adata.obs['nFeature_RNA'].astype(int)

print(adata.obs['site'].value_counts())
print(adata.obs['type'].value_counts())
print(adata.obs['cancer'].value_counts())

adata.obs['DIG_Score1']=adata.obs['DIG.Score1'].copy()
adata.obs['DIG_Score1']=adata.obs['DIG_Score1'].astype(float)
del adata.obs['DIG.Score1']

adata.obs['S_Score']=adata.obs['S.Score'].copy()
adata.obs['S_Score']=adata.obs['S_Score'].astype(float)
del adata.obs['S.Score']

adata.obs['G2M_Score']=adata.obs['G2M.Score'].copy()
adata.obs['G2M_Score']=adata.obs['G2M_Score'].astype(float)
del adata.obs['G2M.Score']

print(adata)
adata.write('/home/bcd/revision_nature/B_Cell_Atlas/b_cell_atlas_processed.h5ad')

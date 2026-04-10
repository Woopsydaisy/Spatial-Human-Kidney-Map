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
from google.colab import drive
import leidenalg as la
from pathlib import Path

##this is not imputed but normal measured expession!

adata=sc.read_h5ad('/content/drive/MyDrive/Bernhard/revision_nature/revision_object_for_colab.h5ad')

# Assuming 'adata' is your AnnData object
markers = ['CALB1','SLC8A1','SPP1','SLC12A3','GSTP1','PECAM1','VIM','TM4SF1', 'KLF2','RGS5','CCL21','PLVAP', 'IGFBP5', 'TGFBR2',  'COL1A1', 'COL1A2', 'SLC26A7',
           'KRT7','SLC4A9', 'HLA-DRB', 'HLA-DRA', 'CXCR4','PDGFRB','ENG','GATA3', 'AQP2', 'AQP3', 'DUSP1', 'IGFBP7',  'ACSM2B', 'LRP2', 'CUBN',
           'NPHS2', 'VEGFA','SPOCK2', 'UMOD', 'SLC12A1', 'TAGLN', 'ACTA2','NOTCH3', 'CXCL14', 'IL32','MMP7','ITGB6', 'S100A6']

# Ensure markers list only includes valid gene names
markers = [gene for gene in markers if gene in adata.var_names]

##now the cosmx

adata_cosmx=adata[adata.obs['tech']=='CosMx'].copy()
adata_cosmx

adata_cosmx.X=adata_cosmx.layers['counts'].copy()

adata_cosmx=adata_cosmx[adata_cosmx.obs['annotation_updated']!='EC_Lymph'].copy()
adata_cosmx

sc.pp.log1p(adata_cosmx)

# Assuming 'adata' is your AnnData object
markers = ['CALB1','SLC8A1','SPP1','SLC12A3','GSTP1','PECAM1','VIM','TM4SF1', 'KLF2','RGS5','PLVAP', 'IGFBP5', 'TGFBR2',  'COL1A1', 'COL1A2', 'SLC26A7',
           'KRT7','SLC4A9', 'C1QA','C1QB', 'CXCR4','PDGFRB','GATA3', 'AQP2', 'AQP3', 'DUSP1', 'IGFBP7',  'ACSM2B', 'LRP2', 'CUBN',
           'NPHS2', 'VEGFA','SPOCK2', 'UMOD', 'SLC12A1', 'TAGLN', 'ACTA2','NOTCH3', 'CXCL14', 'IL32','MMP7', 'S100A6']

# Ensure markers list only includes valid gene names
markers = [gene for gene in markers if gene in adata_cosmx.var_names]

# Subset the .X matrix for the markers and convert to a DataFrame
gene_expression_df_imputed = pd.DataFrame(
    adata_cosmx[:, markers].X.toarray(),
    index=adata_cosmx.obs_names,
    columns=markers
)

gene_expression_df_imputed['annotation'] = adata_cosmx.obs['annotation_updated'].values
mean_expression_per_cluster_unimputed = gene_expression_df_imputed.groupby('annotation').mean()

from sklearn.preprocessing import MinMaxScaler

# Create a MinMaxScaler object
scaler = MinMaxScaler()

# Scale each column individually
gene_expression_scaled_imputed = scaler.fit_transform(mean_expression_per_cluster_unimputed)

# Convert the scaled array back to a DataFrame
gene_expression_scaled_df_imputed = pd.DataFrame(gene_expression_scaled_imputed, index=mean_expression_per_cluster_unimputed.index, columns=mean_expression_per_cluster_unimputed.columns)


plt.figure(figsize=(20, 7.5))  # Adjust the figure size as needed
heatmap = sns.heatmap(gene_expression_scaled_df_imputed, cmap='plasma', vmax=1, cbar_kws={"shrink": .3})

# Adjust the size of the colorbar
heatmap.collections[0].colorbar.ax.set_ylabel('Expression level')
heatmap.collections[0].colorbar.ax.set_xlabel('')
# Removing the tick marks as requested (optional)
plt.tick_params(axis='both', which='both', length=0,labelsize=12)  # Ensures no tick marks on both x and y axes

# Save the figure
plt.savefig('suppl_heatmap_cosmx_expression.png', dpi=450, bbox_inches='tight')

# Show the plot
plt.show()

adata_xenium=adata[adata.obs['tech']=='Xenium'].copy()
adata_xenium

adata_xenium.X=adata_xenium.layers['counts'].copy()

sc.pp.log1p(adata_xenium)

# Assuming 'adata' is your AnnData object
markers = ['CALB1','SLC8A1','SPP1','SLC12A3','TACTSD2','PFKFB3','PAX8','PECAM1','TM4SF1','PDPN','FLT4','TFPI','PLVAP','TGFBR2',  'COL1A1', 'COL1A2','C7', 'SLC26A7',
           'ATP6AP2','ATP6V1B1', 'C1QA','C1QB', 'CXCR4','PDGFRB','IL1RL1','GATA3', 'AQP2', 'AQP3', 'ITGB3',   'ACSM2B', 'LRP2', 'CUBN','ALDOB',
           'VEGFA','PLA2R1','SPOCK2','UMOD', 'SLC12A1', 'TAGLN', 'ACTA2','NOTCH3', 'IL32','MMP7','ITGB6', 'SLPI']

# Ensure markers list only includes valid gene names
markers = [gene for gene in markers if gene in adata_xenium.var_names]

# Subset the .X matrix for the markers and convert to a DataFrame
gene_expression_df_imputed = pd.DataFrame(
    adata_xenium[:, markers].X.toarray(),
    index=adata_xenium.obs_names,
    columns=markers
)

gene_expression_df_imputed['annotation'] = adata_xenium.obs['annotation_updated'].values
mean_expression_per_cluster_unimputed = gene_expression_df_imputed.groupby('annotation').mean()

from sklearn.preprocessing import MinMaxScaler

# Create a MinMaxScaler object
scaler = MinMaxScaler()

# Scale each column individually
gene_expression_scaled_imputed = scaler.fit_transform(mean_expression_per_cluster_unimputed)

# Convert the scaled array back to a DataFrame
gene_expression_scaled_df_imputed = pd.DataFrame(gene_expression_scaled_imputed, index=mean_expression_per_cluster_unimputed.index, columns=mean_expression_per_cluster_unimputed.columns)


plt.figure(figsize=(20, 7.5))  # Adjust the figure size as needed
heatmap = sns.heatmap(gene_expression_scaled_df_imputed, cmap='plasma', vmax=1, cbar_kws={"shrink": .3})

# Adjust the size of the colorbar
heatmap.collections[0].colorbar.ax.set_ylabel('Expression level')
heatmap.collections[0].colorbar.ax.set_xlabel('')
# Removing the tick marks as requested (optional)
plt.tick_params(axis='both', which='both', length=0,labelsize=12)  # Ensures no tick marks on both x and y axes

# Save the figure
plt.savefig('suppl_heatmap_xenium_expression.png', dpi=450, bbox_inches='tight')

# Show the plot
plt.show()

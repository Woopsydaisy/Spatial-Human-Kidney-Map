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
import scipy

adata=sc.read_h5ad('/content/drive/MyDrive/Bernhard/revision_nature/revision_object.h5ad')
adata_cosmx=adata[adata.obs['tech']=='CosMx']
adata_cosmx
adata_cosmx=adata_cosmx[adata_cosmx.obs['niches_annotation_based']!='Unknown']

adata_cosmx.X = adata_cosmx.layers['counts'].copy()  

# 2. Extract relevant identifiers
sample_ids = adata_cosmx.obs['unique_sample_identifier'].astype(str)
group_ids = adata_cosmx.obs['niches_annotation_based'].astype(str)

# 3. Create a dataframe of expression values
X = adata_cosmx.X.toarray()
expr_df = pd.DataFrame(X, columns=adata_cosmx.var_names)
expr_df = expr_df.loc[:, expr_df.sum(axis=0) != 0]
expr_df['sample'] = sample_ids.values
expr_df['group'] = group_ids.values

# 4. Aggregate (sum) by sample and niche (group) — pseudobulk
pseudobulk = expr_df.groupby(['sample', 'group']).sum()

# 5. Reset index and create a combined index of "sample_group"
pseudobulk.reset_index(inplace=True)
pseudobulk['sample_group'] = pseudobulk['sample'] + '_' + pseudobulk['group']
pseudobulk_matrix = pseudobulk.drop(columns=['sample', 'group']).set_index('sample_group')

# 6. Create metadata using the same combined index
sample_metadata = pseudobulk[['sample_group']].copy()
sample_metadata['group'] = pseudobulk['group']
sample_metadata = sample_metadata.set_index('sample_group')

pseudobulk_matrix.to_csv('pseudobulk_counts_niches.csv')
sample_metadata.to_csv('pseudobulk_metadata_niches.csv')

##after running DESeq2:
CNT_niche=pd.read_csv('/content/DE_niche_vs_rest_CNT niche.csv')
CNT_niche['padj'] = CNT_niche['padj'].replace(0, 1e-300)
CNT_niche=CNT_niche[CNT_niche['padj']<0.05]
CNT_niche=CNT_niche[CNT_niche['log2FoldChange']>0.5]
CNT_niche['pi_score'] = np.sign(CNT_niche['log2FoldChange']) * -np.log10(CNT_niche['padj'])
CNT_niche=CNT_niche.sort_values(by='pi_score', ascending=False)

Fibro_niche=pd.read_csv('/content/DE_niche_vs_rest_Fibroblast niche.csv')
Fibro_niche['padj'] = Fibro_niche['padj'].replace(0, 1e-300)
Fibro_niche=Fibro_niche[Fibro_niche['padj']<0.05]
Fibro_niche=Fibro_niche[Fibro_niche['log2FoldChange']>0.5]
Fibro_niche['pi_score'] = np.sign(Fibro_niche['log2FoldChange']) * -np.log10(Fibro_niche['padj'])
Fibro_niche=Fibro_niche.sort_values(by='pi_score', ascending=False)

Glom_niche=pd.read_csv('/content/DE_niche_vs_rest_Glomerular niche.csv')
Glom_niche['padj'] = Glom_niche['padj'].replace(0, 1e-300)
Glom_niche=Glom_niche[Glom_niche['padj']<0.05]
Glom_niche=Glom_niche[Glom_niche['log2FoldChange']>0.5]
Glom_niche['pi_score'] = np.sign(Glom_niche['log2FoldChange']) * -np.log10(Glom_niche['padj'])
Glom_niche=Glom_niche.sort_values(by='pi_score', ascending=False)

Immune_niche=pd.read_csv('/content/DE_niche_vs_rest_Immune niche.csv')
Immune_niche['padj'] = Immune_niche['padj'].replace(0, 1e-300)
Immune_niche=Immune_niche[Immune_niche['padj']<0.05]
Immune_niche=Immune_niche[Immune_niche['log2FoldChange']>0.5]
Immune_niche['pi_score'] = np.sign(Immune_niche['log2FoldChange']) * -np.log10(Immune_niche['padj'])
Immune_niche=Immune_niche.sort_values(by='pi_score', ascending=False)

LOH_niche=pd.read_csv('/content/DE_niche_vs_rest_LOH niche.csv')
LOH_niche['padj'] = LOH_niche['padj'].replace(0, 1e-300)
LOH_niche=LOH_niche[LOH_niche['padj']<0.05]
LOH_niche=LOH_niche[LOH_niche['log2FoldChange']>0.5]
LOH_niche['pi_score'] = np.sign(LOH_niche['log2FoldChange']) * -np.log10(LOH_niche['padj'])
LOH_niche=LOH_niche.sort_values(by='pi_score', ascending=False)

PC_niche=pd.read_csv('/content/DE_niche_vs_rest_PC niche.csv')
PC_niche['padj'] = PC_niche['padj'].replace(0, 1e-300)
PC_niche=PC_niche[PC_niche['padj']<0.05]
PC_niche=PC_niche[PC_niche['log2FoldChange']>0.5]
PC_niche['pi_score'] = np.sign(PC_niche['log2FoldChange']) * -np.log10(PC_niche['padj'])
PC_niche=PC_niche.sort_values(by='pi_score', ascending=False)

PT_niche=pd.read_csv('/content/DE_niche_vs_rest_PT niche.csv')
PT_niche['padj'] = PT_niche['padj'].replace(0, 1e-300)
PT_niche=PT_niche[PT_niche['padj']<0.05]
PT_niche=PT_niche[PT_niche['log2FoldChange']>0.5]
PT_niche['pi_score'] = np.sign(PT_niche['log2FoldChange']) * -np.log10(PT_niche['padj'])
PT_niche=PT_niche.sort_values(by='pi_score', ascending=False)

Vascular_niche=pd.read_csv('/content/DE_niche_vs_rest_Vascular niche.csv')
Vascular_niche['padj'] = Vascular_niche['padj'].replace(0, 1e-300)
Vascular_niche=Vascular_niche[Vascular_niche['padj']<0.05]
Vascular_niche=Vascular_niche[Vascular_niche['log2FoldChange']>0.5]
Vascular_niche['pi_score'] = np.sign(Vascular_niche['log2FoldChange']) * -np.log10(Vascular_niche['padj'])
Vascular_niche=Vascular_niche.sort_values(by='pi_score', ascending=False)

iLOH_niche=pd.read_csv('/content/DE_niche_vs_rest_iLOH niche.csv')
iLOH_niche['padj'] = iLOH_niche['padj'].replace(0, 1e-300)
iLOH_niche=iLOH_niche[iLOH_niche['padj']<0.05]
iLOH_niche=iLOH_niche[iLOH_niche['log2FoldChange']>0.5]
iLOH_niche['pi_score'] = np.sign(iLOH_niche['log2FoldChange']) * -np.log10(iLOH_niche['padj'])
iLOH_niche=iLOH_niche.sort_values(by='pi_score', ascending=False)

iPT_niche=pd.read_csv('/content/DE_niche_vs_rest_iPT niche.csv')
iPT_niche['padj'] = iPT_niche['padj'].replace(0, 1e-300)
iPT_niche=iPT_niche[iPT_niche['padj']<0.05]
iPT_niche=iPT_niche[iPT_niche['log2FoldChange']>0.5]
iPT_niche['pi_score'] = np.sign(iPT_niche['log2FoldChange']) * -np.log10(iPT_niche['padj'])
iPT_niche=iPT_niche.sort_values(by='pi_score', ascending=False)

DCT_niche=pd.read_csv('/content/DE_niche_vs_rest_DCT niche.csv')
DCT_niche['padj'] = DCT_niche['padj'].replace(0, 1e-300)
DCT_niche=DCT_niche[DCT_niche['padj']<0.05]
DCT_niche=DCT_niche[DCT_niche['log2FoldChange']>0.5]
DCT_niche['pi_score'] = np.sign(DCT_niche['log2FoldChange']) * -np.log10(DCT_niche['padj'])
DCT_niche=DCT_niche.sort_values(by='pi_score', ascending=False)


# Dictionary of your loaded dataframes with their labels
niche_dfs = {
    "CNT niche": CNT_niche,
    "Fibroblast niche": Fibro_niche,
    "Glomerular niche": Glom_niche,
    "Immune niche": Immune_niche,
    "LOH niche": LOH_niche,
    "PC niche": PC_niche,
    "PT niche": PT_niche,
    "Vascular niche": Vascular_niche,
    "iLOH niche": iLOH_niche,
    "iPT niche": iPT_niche,
    "DCT niche": DCT_niche
}

output_file = "/content/spatial_genesignatures.csv"

with open(output_file, "w") as f_out:
    for label, df in niche_dfs.items():
        f_out.write(f"Top 10 DEGs - {label}\n")
        df.head(10).to_csv(f_out, index=False)  # index=False if 'gene' is a column
        f_out.write("\n")  # blank line between sections

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
import leidenalg as la
from pathlib import Path
import gc

sc.settings.figdir = "/home/bcd/revision_nature/neighbors_neighborhood/plots_script8"  # Update this path as needed
sc.settings.set_figure_params(dpi=900)

adata = sc.read_h5ad('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad')
print(adata)
print(adata.obs['niches_annotation_based'].value_counts())
print(adata.obs['unique_sample_identifier'].value_counts())
print(adata.obs['Condition'].value_counts())

adata_imputed=sc.read_h5ad('/home/bcd/revision_nature/imputation_new_atlas/xenium_cosmx_sn_imputed.h5ad')
print(adata_imputed)
adata_imputed.obs['Condition']=adata.obs['Condition'].copy()
adata_imputed.obs['niches_annotation_based']=adata.obs['niches_annotation_based'].copy()

import gc
import scanpy as sc
import pandas as pd

niches = adata_imputed.obs['niches_annotation_based'].unique()
conditions_to_compare = ['DKD', 'Control']

for niche in niches:
    print(f"Processing niche: {niche}")
    
    # Subset for current niche
    adata_niche = adata_imputed[adata_imputed.obs['niches_annotation_based'] == niche].copy()
    
    # Subset for DKD and Control only
    adata_niche = adata_niche[adata_niche.obs['Condition'].isin(conditions_to_compare)].copy()
    
    if adata_niche.n_obs < 10:
        print(f"Skipping {niche} due to low cell count ({adata_niche.n_obs})")
        continue

    # Perform DE analysis
    sc.tl.rank_genes_groups(adata_niche, groupby='Condition', method='wilcoxon', pts=True)
    
    # Save rank genes plot
    sc.pl.rank_genes_groups(
        adata_niche,
        n_genes=25,
        sharey=False,
        save=f'_{niche.replace(" ", "_")}_control_vs_DKD_imputed.png'
    )

    # Extract DEGs and save CSV for each comparison
    result = adata_niche.uns['rank_genes_groups']
    categories = result['names'].dtype.names
    for category in categories:
        df_group = pd.DataFrame({
            'gene': result['names'][category],
            'logfoldchanges': result['logfoldchanges'][category],
            'pvals': result['pvals'][category],
            'pvals_adj': result['pvals_adj'][category],
            'scores': result['scores'][category],
            'pts': result['pts'][category],
            'pts_rest': result['pts_rest'][category]
        })
        df_group.to_csv(
            f'/home/bcd/revision_nature/neighbors_neighborhood/DEGS/niche_{niche.replace(" ", "_")}_degs_{category}.csv',
            index=False
        )
    
    del adata_niche
    gc.collect()

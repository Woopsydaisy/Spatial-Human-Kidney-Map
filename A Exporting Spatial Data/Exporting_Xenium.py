from anndata import AnnData
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import spatialdata as sd
from spatialdata_io import xenium
import spatialdata_plot
import os

neighbors = 15
nPCs = 30
leiden_resolution = 0.5
output_dir = "/home/bcd/Xenium_Rawdata/export_plots"
os.makedirs(output_dir, exist_ok=True)


xenium_path = ""
sdata = xenium(xenium_path)
print(sdata)
adata = sdata.tables["table"]
print(adata)

adata.obs['n_genes_by_counts'] = (adata.X > 0).sum(axis=1).A1

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata.obs["transcript_counts"],
    kde=False,
    ax=axs[0],
)

axs[1].set_title("Unique transcripts per cell")
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1],
)


axs[2].set_title("Area of segmented cells")
sns.histplot(
    adata.obs["cell_area"],
    kde=False,
    ax=axs[2],
)

axs[3].set_title("Nucleus ratio")
sns.histplot(
    adata.obs["nucleus_area"] / adata.obs["cell_area"],
    kde=False,
    ax=axs[3],
)
output_path = os.path.join(output_dir, "QC_metrics_1.png")
fig.savefig(output_path, dpi=450, bbox_inches="tight")

adata.obs['sample']='1'
print(adata)
adata.write('/home/bcd/Xenium_Rawdata/adata_objects/Xenium1_raw.h5ad')

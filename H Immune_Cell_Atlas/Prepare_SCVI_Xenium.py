from anndata import AnnData
import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

xenium_data_annotated=sc.read_h5ad('/home/bcd/revision_nature/integration_new_atlas/integration_combined/revision_atlas_cosmx_xenium_ver1_scvi_scanvi.h5ad')
print(xenium_data_annotated)
xenium_data_annotated=xenium_data_annotated[xenium_data_annotated.obs['annotation_postscanvi_level2']=='Immune'].copy()
xenium_data_annotated=xenium_data_annotated[xenium_data_annotated.obs['tech']!='CosMx'].copy()
print(xenium_data_annotated)

# Map annotations from the dataframe
df=pd.read_csv('/home/bcd/revision_nature/immune_cell_atlas_xenium/annotation_SN_SEQ.csv')
df=df.set_index('Unnamed: 0')
print(df)
# Assuming `df` is your DataFrame with annotations and its index matches cell IDs
annotation_dict = df['annotations_after_scanvi_simple'].to_dict()

xenium_data_filtered.obs['annotations_after_scanvi_simple']='Unknown'
xenium_data_filtered.obs['annotations_after_scanvi_simple'] = xenium_data_filtered.obs_names.map(annotation_dict)
xenium_data_filtered.obs['annotations_after_scanvi_simple']=xenium_data_filtered.obs['annotations_after_scanvi_simple'].astype('category')

xenium_data_filtered.write('/home/bcd/revision_nature/immune_cell_atlas_xenium/immune_cells_xenium_sn_prescvi.h5ad')

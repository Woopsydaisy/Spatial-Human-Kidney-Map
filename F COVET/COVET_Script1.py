import os

os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "0"  

import jax

print("JAX devices:", jax.devices())

import importlib.util
import sys

import warnings
import numpy as np
import pandas as pd
import scipy.sparse
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.neighbors
import scanpy as sc
import umap.umap_ as umap
from fa2 import ForceAtlas2
import importlib.util
import sys
import os

import sys
print("message")
sys.stdout.flush()


# Set the path to the scenvi folder
scenvi_path = "/home/liranmao/katalin/envi/0_revision/envi_new_object/ENVI/scenvi"  # update this path
module_name = "scenvi"

# Load the module from the specific path
spec = importlib.util.spec_from_file_location(module_name, os.path.join(scenvi_path, "__init__.py"))
scenvi = importlib.util.module_from_spec(spec)
sys.modules[module_name] = scenvi
spec.loader.exec_module(scenvi)

print(scenvi.__file__)  

st_data = sc.read_h5ad('/home/liranmao/katalin/envi/0_revision/envi_new_object/revision_object.h5ad')

st_data = st_data[st_data.obs['tech'] == 'CosMx']
st_data=st_data[st_data.obs['niches_annotation_based']!='Unknown']  

import scipy.sparse

# Only apply log1p to sparse matrix without converting to dense
st_data.layers["log"] = st_data.X.copy()
st_data.layers["log"].data = np.log1p(st_data.layers["log"].data)



###### get the ENVI matrix
# Assuming st_data is an AnnData object and st_data.X is a sparse matrix
st_data.obsm['spatial'] = np.array(st_data.obsm['spatial_fov'])


# Convert sparse matrix to dense before addition
if scipy.sparse.issparse(st_data.X):
    st_data.X = st_data.X.toarray().astype(np.float16)


print('starting computing COVET')
# Compute the covariance matrix using ENVI
# st_data.obsm['COVET'], st_data.obsm['COVET_SQRT'], st_data.uns['CovGenes'] = scenvi.compute_covet(st_data, batch_key = 'unique_sample_identifier')

st_data.obsm['COVET'], st_data.obsm['COVET_SQRT'], st_data.uns['CovGenes'] = scenvi.compute_covet(st_data,  batch_key='unique_sample_identifier', batch_size=50000)

st_data.write_h5ad('/home/liranmao/katalin/envi/0_revision/envi_new_object/revision_object_for_colab_ENVI_envi_all.h5ad')

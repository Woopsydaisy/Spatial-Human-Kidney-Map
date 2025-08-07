import os
import sys
import warnings

import numpy as np
import pandas as pd
import scipy.sparse
import cupy as cp

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import sklearn.neighbors
from sklearn.decomposition import PCA

import umap
import umap.umap_ as umap 
import pynndescent
from fa2 import ForceAtlas2

print("message")
sys.stdout.flush()

st_data = sc.read_h5ad('/home/liranmao/katalin/envi/0_revision/envi_new_object/revision_object_for_colab_ENVI_envi_all.h5ad')
st_data
# Define the mapping from Condition to condition_plot

condition_map = {
    'DKD': 'Disease',
    'Control': 'Healthy'
}

# Apply the mapping to create the new obs column
st_data.obs['condition_plot'] = st_data.obs['Condition'].map(condition_map)

# Optional: check for unmapped conditions
unmapped = st_data.obs.loc[st_data.obs['condition_plot'].isna(), 'Condition'].unique()
if len(unmapped) > 0:
    print("Warning: These conditions were not mapped:", unmapped)

### Utility Functions
def flatten(arr):
    return np.asarray(arr).reshape(arr.shape[0], -1)

def fast_knn_pynndescent(data, k):
    index = pynndescent.NNDescent(data, n_neighbors=k, metric='euclidean', n_jobs=-1)
    neighbors, distances = index.neighbor_graph
    N = data.shape[0]
    rows = np.repeat(np.arange(N), k)
    cols = neighbors.flatten()
    dists = distances.flatten()
    kNN = scipy.sparse.csr_matrix((dists, (rows, cols)), shape=(N, N))
    return kNN

def run_diffusion_maps_gpu(data, n_components=10, knn=30, alpha=0):
    if isinstance(data, np.ndarray):
        data = pd.DataFrame(data)

    if not scipy.sparse.issparse(data):
        print("Building fast kNN graph with pynndescent...")
        kNN = fast_knn_pynndescent(data.values, k=knn)

        adaptive_k = int(np.floor(knn / 3))
        kNN_adaptive = fast_knn_pynndescent(data.values, k=adaptive_k)
        adaptive_std = kNN_adaptive.max(axis=1).toarray().ravel()

        x, y, dists = scipy.sparse.find(kNN)
        dists = dists / adaptive_std[x]
        N = data.shape[0]
        W = scipy.sparse.csr_matrix((np.exp(-dists), (x, y)), shape=(N, N))

        kernel = W + W.T
    else:
        kernel = data

    N = kernel.shape[0]
    D = np.ravel(kernel.sum(axis=1))
    if alpha > 0:
        D[D != 0] = D[D != 0] ** (-alpha)
        mat = scipy.sparse.csr_matrix((D, (range(N), range(N))), shape=[N, N])
        kernel = mat.dot(kernel).dot(mat)
        D = np.ravel(kernel.sum(axis=1))

    D[D != 0] = 1 / D[D != 0]
    T = scipy.sparse.csr_matrix((D, (range(N), range(N))), shape=[N, N]).dot(kernel)

    print("Running GPU-accelerated SVD with CuPy...")
    T_dense = T.toarray().astype(np.float32)
    T_gpu = cp.asarray(T_dense)

    U, S, VT = cp.linalg.svd(T_gpu, full_matrices=False)

    U = U[:, :n_components]
    S = S[:n_components]

    res = {
        "T": T,
        "EigenVectors": pd.DataFrame(cp.asnumpy(U), index=data.index if isinstance(data, pd.DataFrame) else None),
        "EigenValues": pd.Series(cp.asnumpy(S))
    }
    return res

def compute_umap_layout_after_pca(data, n_components_pca=50, n_neighbors_umap=30):
    """Apply PCA before UMAP to avoid memory explosion."""
    print(f"Running PCA to reduce to {n_components_pca} components...")

    # Convert to float32 to reduce memory usage
    data = data.astype(np.float32)

    pca = PCA(n_components=n_components_pca, random_state=42)
    data_pca = pca.fit_transform(data)

    print(f"Running UMAP on reduced PCA data...")
    reducer = umap.UMAP(
        n_neighbors=n_neighbors_umap,
        min_dist=0.1,
        spread=1.0,
        metric='euclidean',
        n_components=2,
        random_state=42
    )
    layout = reducer.fit_transform(data_pca)
    return pd.DataFrame(layout, columns=['x', 'y']), data_pca


### Main Code

print('Starting UMAP layout and Diffusion Maps computation (safe PCA+UMAP)...')

flattened_covet = flatten(st_data.obsm['COVET_SQRT'])

# PCA + UMAP Layout
UMAP_COVET, flattened_covet_small = compute_umap_layout_after_pca(flattened_covet, n_components_pca=50, n_neighbors_umap=30)
st_data.obsm['UMAP_COVET'] = np.asarray(UMAP_COVET)

print('Finished computing UMAP layout.')
# save umap
import numpy as np
import pandas as pd
import pickle

# Method 1: Save as NumPy array (recommended for simple reloading)
def save_umap_numpy(st_data, filename="umap_covet_coordinates.npy"):
    """Save UMAP coordinates as NumPy array"""
    np.save(filename, st_data.obsm['UMAP_COVET'])
    print(f"UMAP coordinates saved to {filename}")

def load_umap_numpy(st_data, filename="umap_covet_coordinates.npy"):
    """Load UMAP coordinates from NumPy array"""
    umap_coords = np.load(filename)
    st_data.obsm['UMAP_COVET'] = umap_coords
    print(f"UMAP coordinates loaded from {filename}")
    return st_data

# Method 2: Save as CSV (human-readable, good for sharing)
def save_umap_csv(st_data, filename="umap_covet_coordinates.csv"):
    """Save UMAP coordinates as CSV with cell IDs"""
    umap_df = pd.DataFrame(
        st_data.obsm['UMAP_COVET'], 
        columns=['UMAP1', 'UMAP2'],
        index=st_data.obs_names
    )
    umap_df.to_csv(filename)
    print(f"UMAP coordinates saved to {filename}")

def load_umap_csv(st_data, filename="umap_covet_coordinates.csv"):
    """Load UMAP coordinates from CSV"""
    umap_df = pd.read_csv(filename, index_col=0)
    # Ensure the order matches st_data.obs_names
    umap_df = umap_df.reindex(st_data.obs_names)
    st_data.obsm['UMAP_COVET'] = umap_df.values
    print(f"UMAP coordinates loaded from {filename}")
    return st_data

# Method 3: Save as pickle (preserves exact format)
def save_umap_pickle(st_data, filename="umap_covet_coordinates.pkl"):
    """Save UMAP coordinates as pickle"""
    with open(filename, 'wb') as f:
        pickle.dump(st_data.obsm['UMAP_COVET'], f)
    print(f"UMAP coordinates saved to {filename}")

def load_umap_pickle(st_data, filename="umap_covet_coordinates.pkl"):
    """Load UMAP coordinates from pickle"""
    with open(filename, 'rb') as f:
        umap_coords = pickle.load(f)
    st_data.obsm['UMAP_COVET'] = umap_coords
    print(f"UMAP coordinates loaded from {filename}")
    return st_data

# Method 4: Save entire AnnData object (most comprehensive)
def save_full_anndata(st_data, filename="st_data_with_umap.h5ad"):
    """Save the entire AnnData object including UMAP"""
    st_data.write(filename)
    print(f"Full AnnData object saved to {filename}")

def load_full_anndata(filename="st_data_with_umap.h5ad"):
    """Load the entire AnnData object"""
    import scanpy as sc
    st_data = sc.read(filename)
    print(f"Full AnnData object loaded from {filename}")
    return st_data


# Save using your preferred method
save_umap_numpy(st_data, filename = './0_envi_16reso_results/umap_covet_coordinates_min_dist_0.1.npy')           # Simplest, fastest



def plot_umap_layout_by_condition(st_data, color_by="niches_annotation_based", condition_col="Condition", save_prefix="niche_umap"):
    """
    Plot the UMAP layout stored in st_data.obsm['UMAP_COVET'] separately for each condition,
    and also generate a global plot for all cells.

    Args:
        st_data: AnnData object
        color_by: Column name in st_data.obs to color points by
        condition_col: Column name for conditions to separate plots by
        save_prefix: Prefix for saved plot files
    """
    assert "UMAP_COVET" in st_data.obsm, "UMAP_COVET not found in st_data.obsm!"
    assert color_by in st_data.obs.columns, f"{color_by} not found in st_data.obs!"
    assert condition_col in st_data.obs.columns, f"{condition_col} not found in st_data.obs!"
    
    # Get overall UMAP coordinates and set expanded axis limits
    layout = pd.DataFrame(st_data.obsm['UMAP_COVET'], columns=["UMAP1", "UMAP2"], index=st_data.obs_names)
    
    # Calculate expanded axis limits (add 20% padding to overall range)
    x_min, x_max = layout["UMAP1"].min(), layout["UMAP1"].max()
    y_min, y_max = layout["UMAP2"].min(), layout["UMAP2"].max()
    
    x_padding = (x_max - x_min) * 0.2
    y_padding = (y_max - y_min) * 0.2
    
    x_lim = (x_min - x_padding, x_max + x_padding)
    y_lim = (y_min - y_padding, y_max + y_padding)
    
    # === GLOBAL PLOT ===
    plt.figure(figsize=(10, 8))
    
    metadata_all = st_data.obs[color_by]
    sns.scatterplot(
        x=layout["UMAP1"],
        y=layout["UMAP2"],
        hue=metadata_all,
        s=1,
        linewidth=0,
        palette="tab20",
        alpha=0.8,
        edgecolor=None,
        legend='full'
    )
    
    plt.title(f"Global UMAP Layout colored by {color_by}", fontsize=14)
    plt.xlabel("UMAP 1", fontsize=12)
    plt.ylabel("UMAP 2", fontsize=12)
    plt.xlim(x_lim)
    plt.ylim(y_lim)
    plt.legend(markerscale=4, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., ncol=1)
    plt.tight_layout()
    
    save_path_global = f"{save_prefix}_all.pdf"
    plt.savefig(save_path_global, dpi=300, bbox_inches='tight')
    print(f"Global plot saved to {save_path_global}")
    plt.show()
    
    # === CONDITION-WISE PLOTS ===
    for condition in st_data.obs[condition_col].unique():
        plt.figure(figsize=(10, 8))
        
        condition_mask = st_data.obs[condition_col] == condition
        condition_data = st_data[condition_mask]
        
        condition_layout = pd.DataFrame(
            condition_data.obsm['UMAP_COVET'], 
            columns=["UMAP1", "UMAP2"], 
            index=condition_data.obs_names
        )
        condition_metadata = condition_data.obs[color_by]
        
        sns.scatterplot(
            x=condition_layout["UMAP1"], 
            y=condition_layout["UMAP2"],
            hue=condition_metadata,
            s=1,
            linewidth=0,
            palette="tab20",
            alpha=0.8,
            edgecolor=None,
            legend='full'
        )
        
        plt.title(f"UMAP Layout colored by {color_by} - Condition: {condition}", fontsize=14)
        plt.xlabel("UMAP 1", fontsize=12)
        plt.ylabel("UMAP 2", fontsize=12)
        plt.xlim(x_lim)
        plt.ylim(y_lim)
        plt.legend(markerscale=4, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., ncol=1)
        plt.tight_layout()
        
        save_path = f"{save_prefix}_{condition}.pdf"
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
        plt.show()

plot_umap_layout_by_condition(st_data, color_by="niches_annotation_based", condition_col="condition_plot", save_prefix="niche_umap")
st_data_sst = st_data_sst[st_data_sst.obs['Condition'].isin(['DKD', 'Control'])]
plot_umap_layout_by_condition(st_data_sst, color_by="niches_annotation_based", condition_col="Condition", save_prefix="niche_umap_DKD_healthy")

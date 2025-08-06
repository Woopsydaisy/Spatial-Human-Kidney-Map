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
import squidpy as sq

adata = sc.read_h5ad('/home/bcd/revision_nature/imputation_new_atlas/xenium_cosmx_sn_imputed.h5ad')
output_dir = "/home/bcd/revision_nature/neighbors_neighborhood/neighbor_histograms_v2"
os.makedirs(output_dir, exist_ok=True)

adata_list = []

for tech in adata.obs['tech'].unique():
    if tech=='CosMx':
        adata_subset=adata[adata.obs['tech']=='CosMx'].copy()
        for i, j in enumerate(adata_subset.obs["orig_ident"].unique()):
            adata_sample = adata_subset[adata_subset.obs["orig_ident"] == j]

            sq.gr.spatial_neighbors(
                adata_sample,
                radius=20/0.12,
                coord_type="generic",
                key_added="20_micron",
                spatial_key = "spatial_fov",
            )
            print(j)
            neighbor_check = np.array(adata_sample.obsp["20_micron_connectivities"].sum(axis=1)).flatten()
            adata_sample.obs["n_neighbors_20um"] = neighbor_check
            plt.hist(neighbor_check, bins=np.arange(0, 30, 1))
            plt.title(f"Neighbors in {j}")
            plt.xlabel("Number of neighbors")
            plt.ylabel("Frequency")
            plt.tight_layout()
            filename = f"{output_dir}/neighbors_{j}_{tech}_20um.png"
            plt.savefig(filename, dpi=300)
            plt.close()


            adata_list.append(adata_sample)
    if tech=='Xenium':
        adata_subset=adata[adata.obs['tech']=='Xenium'].copy()
        for i, j in enumerate(adata_subset.obs["orig_ident"].unique()):
            adata_sample = adata_subset[adata_subset.obs["orig_ident"] == j]

            sq.gr.spatial_neighbors(
                adata_sample,
                radius=20,
                coord_type="generic",
                key_added="20_micron",
                spatial_key = "spatial",
            )
            print(j)
            neighbor_check = np.array(adata_sample.obsp["20_micron_connectivities"].sum(axis=1)).flatten()
            adata_sample.obs["n_neighbors_20um"] = neighbor_check
            plt.hist(neighbor_check, bins=np.arange(0, 30, 1))
            plt.title(f"Neighbors in {j}")
            plt.xlabel("Number of neighbors")
            plt.ylabel("Frequency")
            plt.tight_layout()
            filename = f"{output_dir}/neighbors_{j}_{tech}_20um.png"
            plt.savefig(filename, dpi=300)
            plt.close()
            adata_list.append(adata_sample)

adata_neighbored = ad.concat(adata_list, pairwise = True)
print(adata_neighbored)
print(adata_neighbored.obs["n_neighbors_20um"].describe())

# Transfer '20_micron_connectivities'
adata.obsp['20_micron_connectivities'] = adata_neighbored.obsp['20_micron_connectivities'].copy()

# Transfer '20_micron_distances'
adata.obsp['20_micron_distances'] = adata_neighbored.obsp['20_micron_distances'].copy()

adata.obs['n_neighbors_20um']=adata_neighbored.obs["n_neighbors_20um"].copy()

adata.write('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad')
print('20_done')


#now 40um
adata_list = []

for tech in adata.obs['tech'].unique():
    if tech=='CosMx':
        adata_subset=adata[adata.obs['tech']=='CosMx'].copy()
        for i, j in enumerate(adata_subset.obs["orig_ident"].unique()):
            adata_sample = adata_subset[adata_subset.obs["orig_ident"] == j]

            sq.gr.spatial_neighbors(
                adata_sample,
                radius=40/0.12,
                coord_type="generic",
                key_added="40_micron",
                spatial_key = "spatial_fov",
            )
            print(j)
            neighbor_check = np.array(adata_sample.obsp["40_micron_connectivities"].sum(axis=1)).flatten()
            adata_sample.obs["n_neighbors_40um"] = neighbor_check
            plt.hist(neighbor_check, bins=np.arange(0, 50, 1))
            plt.title(f"Neighbors in {j}")
            plt.xlabel("Number of neighbors")
            plt.ylabel("Frequency")
            plt.tight_layout()
            filename = f"{output_dir}/neighbors_{j}_{tech}_40um.png"
            plt.savefig(filename, dpi=300)
            plt.close()


            adata_list.append(adata_sample)
    if tech=='Xenium':
        adata_subset=adata[adata.obs['tech']=='Xenium'].copy()
        for i, j in enumerate(adata_subset.obs["orig_ident"].unique()):
            adata_sample = adata_subset[adata_subset.obs["orig_ident"] == j]

            sq.gr.spatial_neighbors(
                adata_sample,
                radius=40,
                coord_type="generic",
                key_added="40_micron",
                spatial_key = "spatial",
            )
            print(j)
            neighbor_check = np.array(adata_sample.obsp["40_micron_connectivities"].sum(axis=1)).flatten()
            adata_sample.obs["n_neighbors_40um"] = neighbor_check
            plt.hist(neighbor_check, bins=np.arange(0, 50, 1))
            plt.title(f"Neighbors in {j}")
            plt.xlabel("Number of neighbors")
            plt.ylabel("Frequency")
            plt.tight_layout()
            filename = f"{output_dir}/neighbors_{j}_{tech}_40um.png"
            plt.savefig(filename, dpi=300)
            plt.close()
            adata_list.append(adata_sample)

adata_neighbored = ad.concat(adata_list, pairwise = True)
print(adata_neighbored)
print(adata_neighbored.obs["n_neighbors_40um"].describe())

# Transfer '40_micron_connectivities'
adata.obsp['40_micron_connectivities'] = adata_neighbored.obsp['40_micron_connectivities'].copy()

# Transfer '40_micron_distances'
adata.obsp['40_micron_distances'] = adata_neighbored.obsp['40_micron_distances'].copy()

adata.obs['n_neighbors_40um']=adata_neighbored.obs["n_neighbors_40um"].copy()

adata.write('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad')
print('40_done')



#now 60um
adata_list = []

for tech in adata.obs['tech'].unique():
    if tech=='CosMx':
        adata_subset=adata[adata.obs['tech']=='CosMx'].copy()
        for i, j in enumerate(adata_subset.obs["orig_ident"].unique()):
            adata_sample = adata_subset[adata_subset.obs["orig_ident"] == j]

            sq.gr.spatial_neighbors(
                adata_sample,
                radius=60/0.12,
                coord_type="generic",
                key_added="60_micron",
                spatial_key = "spatial_fov",
            )
            print(j)
            neighbor_check = np.array(adata_sample.obsp["60_micron_connectivities"].sum(axis=1)).flatten()
            adata_sample.obs["n_neighbors_60um"] = neighbor_check
            plt.hist(neighbor_check, bins=np.arange(0, 100, 1))
            plt.title(f"Neighbors in {j}")
            plt.xlabel("Number of neighbors")
            plt.ylabel("Frequency")
            plt.tight_layout()
            filename = f"{output_dir}/neighbors_{j}_{tech}_60um.png"
            plt.savefig(filename, dpi=300)
            plt.close()


            adata_list.append(adata_sample)
    if tech=='Xenium':
        adata_subset=adata[adata.obs['tech']=='Xenium'].copy()
        for i, j in enumerate(adata_subset.obs["orig_ident"].unique()):
            adata_sample = adata_subset[adata_subset.obs["orig_ident"] == j]

            sq.gr.spatial_neighbors(
                adata_sample,
                radius=60,
                coord_type="generic",
                key_added="60_micron",
                spatial_key = "spatial",
            )
            print(j)
            neighbor_check = np.array(adata_sample.obsp["60_micron_connectivities"].sum(axis=1)).flatten()
            adata_sample.obs["n_neighbors_60um"] = neighbor_check
            plt.hist(neighbor_check, bins=np.arange(0, 100, 1))
            plt.title(f"Neighbors in {j}")
            plt.xlabel("Number of neighbors")
            plt.ylabel("Frequency")
            plt.tight_layout()
            filename = f"{output_dir}/neighbors_{j}_{tech}_60um.png"
            plt.savefig(filename, dpi=300)
            plt.close()
            adata_list.append(adata_sample)

adata_neighbored = ad.concat(adata_list, pairwise = True)
print(adata_neighbored)
print(adata_neighbored.obs["n_neighbors_60um"].describe())

# Transfer '60_micron_connectivities'
adata.obsp['60_micron_connectivities'] = adata_neighbored.obsp['60_micron_connectivities'].copy()

# Transfer '60_micron_distances'
adata.obsp['60_micron_distances'] = adata_neighbored.obsp['60_micron_distances'].copy()

adata.obs['n_neighbors_60um']=adata_neighbored.obs["n_neighbors_60um"].copy()

adata.write('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad')
print('60_done')


#now 80um
adata_list = []

for tech in adata.obs['tech'].unique():
    if tech=='CosMx':
        adata_subset=adata[adata.obs['tech']=='CosMx'].copy()
        for i, j in enumerate(adata_subset.obs["orig_ident"].unique()):
            adata_sample = adata_subset[adata_subset.obs["orig_ident"] == j]

            sq.gr.spatial_neighbors(
                adata_sample,
                radius=80/0.12,
                coord_type="generic",
                key_added="80_micron",
                spatial_key = "spatial_fov",
            )
            print(j)
            neighbor_check = np.array(adata_sample.obsp["80_micron_connectivities"].sum(axis=1)).flatten()
            adata_sample.obs["n_neighbors_80um"] = neighbor_check
            plt.hist(neighbor_check, bins=np.arange(0, 150, 1))
            plt.title(f"Neighbors in {j}")
            plt.xlabel("Number of neighbors")
            plt.ylabel("Frequency")
            plt.tight_layout()
            filename = f"{output_dir}/neighbors_{j}_{tech}_80um.png"
            plt.savefig(filename, dpi=300)
            plt.close()


            adata_list.append(adata_sample)
    if tech=='Xenium':
        adata_subset=adata[adata.obs['tech']=='Xenium'].copy()
        for i, j in enumerate(adata_subset.obs["orig_ident"].unique()):
            adata_sample = adata_subset[adata_subset.obs["orig_ident"] == j]

            sq.gr.spatial_neighbors(
                adata_sample,
                radius=80,
                coord_type="generic",
                key_added="80_micron",
                spatial_key = "spatial",
            )
            print(j)
            neighbor_check = np.array(adata_sample.obsp["80_micron_connectivities"].sum(axis=1)).flatten()
            adata_sample.obs["n_neighbors_80um"] = neighbor_check
            plt.hist(neighbor_check, bins=np.arange(0, 150, 1))
            plt.title(f"Neighbors in {j}")
            plt.xlabel("Number of neighbors")
            plt.ylabel("Frequency")
            plt.tight_layout()
            filename = f"{output_dir}/neighbors_{j}_{tech}_80um.png"
            plt.savefig(filename, dpi=300)
            plt.close()
            adata_list.append(adata_sample)

adata_neighbored = ad.concat(adata_list, pairwise = True)
print(adata_neighbored)
print(adata_neighbored.obs["n_neighbors_80um"].describe())

# Transfer '80_micron_connectivities'
adata.obsp['80_micron_connectivities'] = adata_neighbored.obsp['80_micron_connectivities'].copy()

# Transfer '80_micron_distances'
adata.obsp['80_micron_distances'] = adata_neighbored.obsp['80_micron_distances'].copy()

adata.obs['n_neighbors_80um']=adata_neighbored.obs["n_neighbors_80um"].copy()

adata.write('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad')
print('80_done')


#now 40um
adata_list = []

for tech in adata.obs['tech'].unique():
    if tech=='CosMx':
        adata_subset=adata[adata.obs['tech']=='CosMx'].copy()
        for i, j in enumerate(adata_subset.obs["orig_ident"].unique()):
            adata_sample = adata_subset[adata_subset.obs["orig_ident"] == j]

            sq.gr.spatial_neighbors(
                adata_sample,
                radius=200/0.12,
                coord_type="generic",
                key_added="200_micron",
                spatial_key = "spatial_fov",
            )
            print(j)
            neighbor_check = np.array(adata_sample.obsp["200_micron_connectivities"].sum(axis=1)).flatten()
            adata_sample.obs["n_neighbors_200um"] = neighbor_check
            plt.hist(neighbor_check, bins=np.arange(0, 200, 1))
            plt.title(f"Neighbors in {j}")
            plt.xlabel("Number of neighbors")
            plt.ylabel("Frequency")
            plt.tight_layout()
            filename = f"{output_dir}/neighbors_{j}_{tech}_200um.png"
            plt.savefig(filename, dpi=300)
            plt.close()


            adata_list.append(adata_sample)
    if tech=='Xenium':
        adata_subset=adata[adata.obs['tech']=='Xenium'].copy()
        for i, j in enumerate(adata_subset.obs["orig_ident"].unique()):
            adata_sample = adata_subset[adata_subset.obs["orig_ident"] == j]

            sq.gr.spatial_neighbors(
                adata_sample,
                radius=200,
                coord_type="generic",
                key_added="200_micron",
                spatial_key = "spatial",
            )
            print(j)
            neighbor_check = np.array(adata_sample.obsp["200_micron_connectivities"].sum(axis=1)).flatten()
            adata_sample.obs["n_neighbors_200um"] = neighbor_check
            plt.hist(neighbor_check, bins=np.arange(0, 200, 1))
            plt.title(f"Neighbors in {j}")
            plt.xlabel("Number of neighbors")
            plt.ylabel("Frequency")
            plt.tight_layout()
            filename = f"{output_dir}/neighbors_{j}_{tech}_200um.png"
            plt.savefig(filename, dpi=300)
            plt.close()
            adata_list.append(adata_sample)

adata_neighbored = ad.concat(adata_list, pairwise = True)
print(adata_neighbored)
print(adata_neighbored.obs["n_neighbors_200um"].describe())

# Transfer '200_micron_connectivities'
#adata.obsp['200_micron_connectivities'] = adata_neighbored.obsp['200_micron_connectivities'].copy()

# Transfer '200_micron_distances'
#adata.obsp['200_micron_distances'] = adata_neighbored.obsp['200_micron_distances'].copy()

adata.obs['n_neighbors_200um']=adata_neighbored.obs["n_neighbors_200um"].copy()

adata.write('/home/bcd/revision_nature/neighbors_neighborhood/xenium_cosmx_imputed_neighbored.h5ad')

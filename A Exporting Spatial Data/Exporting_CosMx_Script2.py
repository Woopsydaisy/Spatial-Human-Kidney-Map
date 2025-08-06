#!/usr/bin/env python
# coding: utf-8

# In[1]:


from anndata import AnnData
import scanpy as sc
import squidpy as sq
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# ## Set parameters

# In[2]:


adata = sc.read_h5ad("/home/levinsj/spatial/adata/slides/run1_S1_rawExport.h5ad")


# In[3]:


adata.layers["All_counts"] = adata.X.copy()


# In[4]:


sc.pl.scatter(
    adata,
    x="CenterX_global_px",
    y="CenterY_global_px",
    color="Run_Tissue_name",
    size=1,
)


# In[5]:


print(adata)


# partition slide into samples.
# This will require you to augment the following script to ensure that sample name, sample type and location of tissue are correct.

# In[6]:


result = []
type_1 = []
for i in adata.obs.index:
    if (adata.obs.loc[i, "CenterY_global_px"] > 100000):
        result.append("HK3524")
        type_1.append("Fetal")       
    else:
        result.append("HK2753")
        type_1.append("Healthy")
adata.obs["sample"] = result
adata.obs["type"] = type_1


# In[7]:


sc.pl.scatter(
    adata,
    x="CenterX_global_px",
    y="CenterY_global_px",
    color="sample",
    size=1,
)


# In[8]:


sc.pl.scatter(
    adata,
    x="CenterX_global_px",
    y="CenterY_global_px",
    color="type",
    size=1,
)


# Save each sample with and without images

# In[9]:


sample1_adata = adata[adata.obs["sample"] == "HK3524"]
sample2_adata = adata[adata.obs["sample"] == "HK2753"]


# # QC of each sample

# In[10]:


samples = [sample1_adata, sample2_adata]
for i in samples:
    print(i)
    # calculate mean counts per FOV
    fov_mean = []
    fov_id = []

    for j in i.obs["fov"].unique():
        fov_mean.append(np.mean(i[i.obs["fov"] == j].obs["nCount_RNA"]))
        fov_id.append(j)

    fov_dict = {k: v for k, v in zip(fov_id, fov_mean)}

    i.obs["meanCounts_perFOV"] = i.obs['fov'].map(fov_dict)
    print(i.obs["meanCounts_perFOV"])
    
    # plot QC for each sample in space
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12.5, 7.5))
    axes[0,0].scatter(i.obs["CenterX_global_px"], i.obs["CenterY_global_px"], c=(i.obs['nCount_RNA']),vmax = 200,vmin = 60, s = 0.1)
    axes[0,1].scatter(i.obs["CenterX_global_px"], i.obs["CenterY_global_px"], c=(i.obs['nFeature_RNA']),vmax = 200,vmin = 60, s = 0.1)
    axes[0,2].scatter(i.obs["CenterX_global_px"], i.obs["CenterY_global_px"], c=(i.obs['nCount_negprobes']),vmax = 5,vmin = 0, s = 0.1)
    axes[1,0].scatter(i.obs["CenterX_global_px"], i.obs["CenterY_global_px"], c=(i.obs['nCount_falsecode']),vmax = 40,vmin = 0, s = 0.1)
    axes[1,1].scatter(i.obs["CenterX_global_px"], i.obs["CenterY_global_px"], c=(i.obs['meanCounts_perFOV']),vmax = 250,vmin = 60, s = 0.1)
    plt.show()
    
    # plot QC on histograms for each sample
    fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(30, 5))
    axes[0].hist(i.obs["nCount_RNA"], bins = np.arange(0, 1000, 50))
    axes[0].axvline(np.mean(i.obs["nCount_RNA"]), color='black', linestyle='dashed', linewidth=1)
    axes[1].hist(i.obs["nFeature_RNA"], bins = np.arange(0, 400, 10))
    axes[1].axvline(np.mean(i.obs["nFeature_RNA"]), color='black', linestyle='dashed', linewidth=1)
    axes[2].hist(i.obs["nCount_falsecode"], bins = np.arange(0, 120, 2))
    axes[2].axvline(np.mean(i.obs["nCount_falsecode"]), color='black', linestyle='dashed', linewidth=1)
    axes[3].hist(i.obs["nCount_negprobes"], bins = np.arange(0, 20, 0.5))
    axes[3].axvline(np.mean(i.obs["nCount_negprobes"]), color='black', linestyle='dashed', linewidth=1)
    axes[4].hist(i.obs["meanCounts_perFOV"], bins = np.arange(50, 250, 10))
    plt.show()
    
    print("Average counts per cell = " + str(np.mean(i.obs["nCount_RNA"])))
    print("Average features per cell = " + str(np.mean(i.obs["nFeature_RNA"])))
    print("Average Negative Probes per cell = " + str(np.mean(i.obs["nCount_negprobes"])))
    print("Average False Codes per cell = " + str(np.mean(i.obs["nCount_falsecode"])))


# In[11]:


for i in samples:
    sq.gr.spatial_neighbors(
    i,
    radius=0.02,
    coord_type="generic",
    key_added = "20_micron")
    
    sq.gr.spatial_neighbors(
    i,
    radius=0.05,
    coord_type="generic",
    key_added = "50_micron")


# In[12]:


sample1_adata.write_h5ad("/home/levinsj/spatial/adata/individual_Samples/HK3524_raw.h5ad")
sample1_adata.uns["spatial"] = []
sample1_adata.write_h5ad("/home/levinsj/spatial/adata/individual_Samples/HK3524_raw_noImages.h5ad")

sample2_adata.write_h5ad("/home/levinsj/spatial/adata/individual_Samples/HK2753_raw.h5ad")
sample2_adata.uns["spatial"] = []
sample2_adata.write_h5ad("/home/levinsj/spatial/adata/individual_Samples/HK2753_raw_noImages.h5ad")

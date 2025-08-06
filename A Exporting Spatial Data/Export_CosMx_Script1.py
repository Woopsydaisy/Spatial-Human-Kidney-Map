#!/usr/bin/env python
# coding: utf-8

# # s01 imports to create a squidpy object with (and without) linked images

# ### import packages, currently running in apptainer container spatial_version2.sif

# In[1]:


import tiledb
from anndata import AnnData
import scanpy as sc
import pandas as pd
import squidpy as sq
import scipy

import os


# ### AWS S3 pull command for images:
# 
# cd ~/spatial/CosMx_raw_export/run5_11_8_23/S1/
# 
# /home/levinsj/Applications/aws/bin/aws s3 cp s3://ksusztak.cosmx/Run5_20231108/S1/20231108_234000_S4/CellStatsDir/ . --recursive

# # Organization of files to run, load files in following format (using export_V1.2.2)
#     
# ### —> Folder for run/
# ### —————>S0/ (Inside this file you need to add the flat files and contents of cell stats directory— can remove RnD and morphology folders)
# ### ————————-> S0_tx_file.csv
# ### ————————-> S0_exprMat_file.csv
# ### ————————-> S0_fov_positions_file.csv
# ### ————————-> S0_metadata_file.csv
# ### ————————-> CellComposite/
# ### ————————-> CellOvelay/
# ### ————————-> FOV001/
# ### ————————-> FOV002/
# ### ————————-> FOV00…/
# 
# ## This script will automatically move and create the cell Labels and cell compartments

# # Load files
# Edit file location here

# In[2]:


sample_dir = "/home/levinsj/spatial/CosMx_raw_export/run1/S1/"
exp = pd.read_csv("/home/levinsj/spatial/CosMx_raw_export/run1/S1/S1_exprMat_file.csv")
tx = pd.read_csv("/home/levinsj/spatial/CosMx_raw_export/run1/S1/S1_tx_file.csv")
fov = pd.read_csv("/home/levinsj/spatial/CosMx_raw_export/run1/S1/S1_fov_positions_file.csv")


# ## Reformat the files and file structure to work with squidpy's function. 

# In[3]:


exp.drop('cell', axis=1).to_csv('/home/levinsj/spatial/CosMx_raw_export/run1/S1/tmp.csv', index=False)

fov = fov.drop('Slide', axis=1)
fov = fov.rename(columns={'FOV': 'fov'})
fov.to_csv('/home/levinsj/spatial/CosMx_raw_export/run1/S1/tmp2.csv', index=False)

print(tx.CellComp.value_counts())


# ## Now moving the images into thee format where we can use them and re-format the 

# In[4]:


os.system(f"mkdir -p {sample_dir}/CellLabels")
os.system(f"mkdir -p {sample_dir}/CompartmentLabels")

for file in os.listdir(sample_dir):
    if "FOV" in file:
        os.system(f"cp {sample_dir}/{file}/CellLabels* {sample_dir}/CellLabels/")
        os.system(f"cp {sample_dir}/{file}/CompartmentLabels* {sample_dir}/CompartmentLabels/")


# ## Load data into an object

# In[5]:


adata = sq.read.nanostring(
    path=sample_dir,
    counts_file="tmp.csv",
    meta_file="S1_metadata_file.csv",
    fov_file="tmp2.csv",
)

print(adata)


# In[6]:


probeNames = pd.read_csv('/home/levinsj/spatial/adata/probeNames.csv')
print(probeNames["x"])


# add spike in to adata.var object

# In[7]:


spike_in = ["ESRRB", "SLC12A1", "UMOD", "CD247", "SLC8A1", "SNTG1", "SLC12A3", "TRPM6", "ACSL4", "SCN2A",
          "SATB2", "STOX2", "EMCN", "MEIS2", "SEMA3A", "PLVAP", "NEGR1", "SERPINE1", "CSMD1", "SLC26A7",
          "SLC22A7", "SLC4A9", "SLC26A4", "CREB5", "HAVCR1", "REN", "AP1S3", "LAMA3", "NOS1", "PAPPA2",
          "SYNPO2", "RET", "LHX1", "SIX2", "CITED1", "WNT9B", "AQP2", "SCNN1G", "ALDH1A2", "CFH", "NTRK3",
          "WT1", "NPHS2", "PTPRQ", "CUBN", "LRP2", "SLC13A3", "ACSM2B", "SLC4A4", "PARD3", "XIST","UTY"]


# In[8]:


adata.var["custom_probes"] = adata.var_names.isin(spike_in)
adata.var["orig_probes"] = ~adata.var_names.isin(spike_in + ["NegativeAdd"])


# # create layers with location of transcripts
# each layer has counts parsed by location of read. Counts are raw counts.

# In[9]:


exp_cyt = exp

for col in exp_cyt.columns[3:]:
    exp_cyt[col].values[:] = 0

txCyt = tx[tx.CellComp == "Cytoplasm"]

exp_cyt.index = exp_cyt.cell

for i in txCyt.index:
  exp_cyt.loc[txCyt.loc[i, "cell"],txCyt.loc[i, "target"]] += 1

adata.layers["Cytoplasm"] = scipy.sparse.csr_matrix(exp_cyt.loc[:, exp_cyt.columns.isin(adata.var_names)])


# In[10]:


exp_cyt = exp

for col in exp_cyt.columns[3:]:
    exp_cyt[col].values[:] = 0

txCyt = tx[tx.CellComp == "Membrane"]

exp_cyt.index = exp_cyt.cell

for i in txCyt.index:
  exp_cyt.loc[txCyt.loc[i, "cell"],txCyt.loc[i, "target"]] += 1

adata.layers["Membrane"] = scipy.sparse.csr_matrix(exp_cyt.loc[:, exp_cyt.columns.isin(adata.var_names)])


# In[11]:


exp_cyt = exp

for col in exp_cyt.columns[3:]:
    exp_cyt[col].values[:] = 0

txCyt = tx[tx.CellComp == "Nuclear"]

exp_cyt.index = exp_cyt.cell

for i in txCyt.index:
  exp_cyt.loc[txCyt.loc[i, "cell"],txCyt.loc[i, "target"]] += 1

adata.layers["Nuclear"] = scipy.sparse.csr_matrix(exp_cyt.loc[:, exp_cyt.columns.isin(adata.var_names)])


# given concerns with autoflorescence, let's annotate probe by fluor, and add this information to the .var object

# In[12]:


fluor_codes = pd.read_csv("/home/levinsj/spatial/referenceFiles/probe_fluorophores.csv", index_col=0, header=0)


# In[13]:


adata.var["green_counts"] = 4
adata.var["blue_counts"] = 4
adata.var["red_counts"] = 4
adata.var["yellow_counts"] = 4

g_count_3 = list(fluor_codes[fluor_codes["Green Spots"] == 3].index)
g_count_2 = list(fluor_codes[fluor_codes["Green Spots"] == 2].index)
g_count_1 = list(fluor_codes[fluor_codes["Green Spots"] == 1].index)
g_count_0 = list(fluor_codes[fluor_codes["Green Spots"] == 0].index)

b_count_3 = list(fluor_codes[fluor_codes["Blue spots"] == 3].index)
b_count_2 = list(fluor_codes[fluor_codes["Blue spots"] == 2].index)
b_count_1 = list(fluor_codes[fluor_codes["Blue spots"] == 1].index)
b_count_0 = list(fluor_codes[fluor_codes["Blue spots"] == 0].index)

r_count_3 = list(fluor_codes[fluor_codes["Red spots"] == 3].index)
r_count_2 = list(fluor_codes[fluor_codes["Red spots"] == 2].index)
r_count_1 = list(fluor_codes[fluor_codes["Red spots"] == 1].index)
r_count_0 = list(fluor_codes[fluor_codes["Red spots"] == 0].index)

y_count_3 = list(fluor_codes[fluor_codes["Yellow spots"] == 3].index)
y_count_2 = list(fluor_codes[fluor_codes["Yellow spots"] == 2].index)
y_count_1 = list(fluor_codes[fluor_codes["Yellow spots"] == 1].index)
y_count_0 = list(fluor_codes[fluor_codes["Yellow spots"] == 0].index)


# In[14]:


for i in adata.var["green_counts"].index:
    if i in g_count_3:
        adata.var["green_counts"][i] = 3
    if i in g_count_2:
        adata.var["green_counts"][i] = 2
    if i in g_count_1:
        adata.var["green_counts"][i] = 1
    if i in g_count_0:
        adata.var["green_counts"][i] = 0 

for i in adata.var["blue_counts"].index:
    if i in b_count_3:
        adata.var["blue_counts"][i] = 3
    if i in b_count_2:
        adata.var["blue_counts"][i] = 2
    if i in b_count_1:
        adata.var["blue_counts"][i] = 1
    if i in b_count_0:
        adata.var["blue_counts"][i] = 0 

for i in adata.var["red_counts"].index:
    if i in r_count_3:
        adata.var["red_counts"][i] = 3
    if i in r_count_2:
        adata.var["red_counts"][i] = 2
    if i in r_count_1:
        adata.var["red_counts"][i] = 1
    if i in r_count_0:
        adata.var["red_counts"][i] = 0 

for i in adata.var["yellow_counts"].index:
    if i in y_count_3:
        adata.var["yellow_counts"][i] = 3
    if i in y_count_2:
        adata.var["yellow_counts"][i] = 2
    if i in y_count_1:
        adata.var["yellow_counts"][i] = 1
    if i in y_count_0:
        adata.var["yellow_counts"][i] = 0 


# # Save objects

# In[15]:


var_subset = probeNames["x"]
adata = adata[:, var_subset]


# In[16]:


adata.layers["All_counts"] = adata.X.copy()


# In[17]:


print(adata)


# In[18]:


adata.write_h5ad("/home/levinsj/spatial/adata/slides/run1_S1_rawExport.h5ad")

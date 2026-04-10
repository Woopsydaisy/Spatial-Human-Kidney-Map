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

adata=sc.read_h5ad('/content/drive/MyDrive/Bernhard/revision_nature/revision_object_for_colab.h5ad')

# Extract unique samples
unique_samples = adata.obs['unique_sample_identifier'].unique()

# Define the number of rows, each row having 4 UMAP plots
num_rows = (len(unique_samples) + 6) // 7  # Ensures an extra row if number of samples is not multiple of 4

# Create a figure with subplots
fig, axes = plt.subplots(num_rows, 7, figsize=(20, 4 * num_rows))  # Adjust subplot size as needed

# Ensure axes is iterable in case of a single row or column
if num_rows == 1:
    axes = axes[np.newaxis, :]  # Adds an extra dimension to make indexing consistent
if len(unique_samples) < 4 and num_rows > 1:
    axes = axes[:, np.newaxis]

# Flatten axes for easier handling
axes = axes.flatten()

# Iterate through each sample and plot the UMAP
for i, sample in enumerate(unique_samples):
    ax = axes[i]
    # Subset the data to include only the current sample
    sc.pl.umap(adata[adata.obs['unique_sample_identifier'] == sample, :], color='unique_sample_identifier', legend_loc='none',show=False, ax=ax, title=f'{sample}', frameon=False)

# Hide any unused subplots
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

# Adjust layout and spacing
plt.tight_layout()

# Save the figure
plt.savefig('UMAP_per_sample.png', dpi=450)  # You can adjust the DPI as needed

# Show the figure
plt.show()


import matplotlib.pyplot as plt

# Extract unique samples
unique_samples = adata_cosmx.obs['unique_sample_identifier'].unique()

# Define the number of rows, with 3 sample pairs per figure row
num_rows = (len(unique_samples) + 2) // 3  # Ensures an extra row if the number of samples is not a multiple of 3

# Create a figure with subplots, 6 columns per row for nCount_RNA and nFeature_RNA for three samples per figure row
fig, axes = plt.subplots(num_rows, 6, figsize=(18, 3 * num_rows))

# Ensure axes is iterable in case of a single row
if num_rows == 1:
    axes = [axes]

# Flatten the axes array for easy indexing if it's multidimensional
if num_rows > 1:
    axes = axes.flatten()

# Iterate through each sample and plot the histograms
for i, sample in enumerate(unique_samples):
    row = i // 3  # Determine which row in the figure
    col = (i % 3) * 2  # Determine start column in the figure, 0, 2, or 4

    # Indices for subplots
    ax1 = axes[row * 6 + col]      # nCount_RNA plot position
    ax2 = axes[row * 6 + col + 1]  # nFeature_RNA plot position

    # Plot nCount_RNA histogram
    sample_data_count = adata_cosmx.obs[adata_cosmx.obs['unique_sample_identifier'] == sample]['nCount_RNA']
    ax1.hist(sample_data_count, bins=30, color=plt.get_cmap('tab10')(0), edgecolor='black')
    ax1.set_title(f'Counts per Cell {sample}')
    ax1.set_ylabel('Number of Cells')
    ax1.set_xlabel('Counts per Cell')
    ax1.set_xlim(0, 1000)
    ax1.grid(False)

    # Plot nFeature_RNA histogram
    sample_data_feature = adata_cosmx.obs[adata_cosmx.obs['unique_sample_identifier'] == sample]['nFeature_RNA']
    ax2.hist(sample_data_feature, bins=30, color=plt.get_cmap('tab10')(1), edgecolor='black')
    ax2.set_title(f'Features per Cell {sample}')
    ax2.set_ylabel('Number of Cells')
    ax2.set_xlabel('Features per Cell')
    ax2.set_xlim(0, 400)
    ax2.grid(False)

# Adjust layout
plt.tight_layout()

# Save the figure with 900 DPI
plt.savefig('nCount_nFeature_RNA_histograms_grouped_cosmx.png', dpi=450)

# Show the figure
plt.show()


import matplotlib.pyplot as plt

# Extract unique samples
unique_samples = adata_xenium.obs['unique_sample_identifier'].unique()

# Define the number of rows, with 3 sample pairs per figure row
num_rows = (len(unique_samples) + 2) // 3  # Ensures an extra row if the number of samples is not a multiple of 3

# Create a figure with subplots, 6 columns per row for nCount_RNA and nFeature_RNA for three samples per figure row
fig, axes = plt.subplots(num_rows, 6, figsize=(18, 3 * num_rows))

# Ensure axes is iterable in case of a single row
if num_rows == 1:
    axes = [axes]

# Flatten the axes array for easy indexing if it's multidimensional
if num_rows > 1:
    axes = axes.flatten()

# Iterate through each sample and plot the histograms
for i, sample in enumerate(unique_samples):
    row = i // 3  # Determine which row in the figure
    col = (i % 3) * 2  # Determine start column in the figure, 0, 2, or 4

    # Indices for subplots
    ax1 = axes[row * 6 + col]      # nCount_RNA plot position
    ax2 = axes[row * 6 + col + 1]  # nFeature_RNA plot position

    # Plot nCount_RNA histogram
    sample_data_count = adata_xenium.obs[adata_xenium.obs['unique_sample_identifier'] == sample]['nCount_RNA']
    ax1.hist(sample_data_count, bins=30, color=plt.get_cmap('tab10')(0), edgecolor='black')
    ax1.set_title(f'Counts per Cell {sample}')
    ax1.set_ylabel('Number of Cells')
    ax1.set_xlabel('Counts per Cell')
    ax1.set_xlim(0, 2000)
    ax1.grid(False)

    # Plot nFeature_RNA histogram
    sample_data_feature = adata_xenium.obs[adata_xenium.obs['unique_sample_identifier'] == sample]['nFeature_RNA']
    ax2.hist(sample_data_feature, bins=30, color=plt.get_cmap('tab10')(1), edgecolor='black')
    ax2.set_title(f'Features per Cell {sample}')
    ax2.set_ylabel('Number of Cells')
    ax2.set_xlabel('Features per Cell')
    ax2.set_xlim(0, 1000)
    ax2.grid(False)

# Adjust layout
plt.tight_layout()

# Save the figure with 900 DPI
plt.savefig('nCount_nFeature_RNA_histograms_grouped_xenium.png', dpi=450)

# Show the figure
plt.show()




#supplementary figure for SN QC stuff

adata=sc.read_h5ad('/content/drive/MyDrive/Bernhard/SN_SEQ_ATLAS/Human_extended_annotated.h5ad')

adata

import matplotlib.pyplot as plt
import seaborn as sns
import math
import pandas as pd

# Extract data
df = adata.obs[['orig_ident', 'nCount_RNA']].copy()

# Sort or order by group if desired
ordered_groups = df['orig_ident'].unique()  # or sorted(df['orig_ident'].unique())

# Plot settings
n_per_row = 30
total_groups = len(ordered_groups)
n_rows = math.ceil(total_groups / n_per_row)

# Set figure size
fig_width = 16  # adjust as needed
fig_height = 4.8 * n_rows  # increase height per row
fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, fig_height), squeeze=False)

# Flatten axes for easier indexing
axes = axes.flatten()

# Plot in chunks
for i in range(n_rows):
    start = i * n_per_row
    end = min((i + 1) * n_per_row, total_groups)
    groups_subset = ordered_groups[start:end]

    ax = axes[i]
    sns.violinplot(
        data=df[df['orig_ident'].isin(groups_subset)],
        x='orig_ident',
        y='nCount_RNA',
        order=groups_subset,
        ax=ax,
        inner='quartile',
        scale='width'
    )
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel('nCount_RNA')
    ax.tick_params(axis='x', rotation=90)

plt.tight_layout()
plt.savefig('nCount_RNA_violinplot_by_orig_ident.png', dpi=900)
plt.show()


adata

import matplotlib.pyplot as plt
import seaborn as sns
import math
import pandas as pd

# Extract data
df = adata.obs[['orig_ident', 'nFeature_RNA']].copy()

# Sort or order by group if desired
ordered_groups = df['orig_ident'].unique()  # or sorted(df['orig_ident'].unique())

# Plot settings
n_per_row = 30
total_groups = len(ordered_groups)
n_rows = math.ceil(total_groups / n_per_row)

# Set figure size
fig_width = 16  # adjust as needed
fig_height = 4.3 * n_rows  # increase height per row
fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, fig_height), squeeze=False)

# Flatten axes for easier indexing
axes = axes.flatten()

# Plot in chunks
for i in range(n_rows):
    start = i * n_per_row
    end = min((i + 1) * n_per_row, total_groups)
    groups_subset = ordered_groups[start:end]

    ax = axes[i]
    sns.violinplot(
        data=df[df['orig_ident'].isin(groups_subset)],
        x='orig_ident',
        y='nFeature_RNA',
        order=groups_subset,
        ax=ax,
        inner='quartile',
        scale='width'
    )
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel('nFeature_RNA')
    ax.tick_params(axis='x', rotation=90)

plt.tight_layout()
plt.savefig('nFeature_violinplot_by_orig_ident.png', dpi=900)
plt.show()


adata

import matplotlib.pyplot as plt
import seaborn as sns
import math
import pandas as pd

# Extract data
df = adata.obs[['orig_ident', 'percent_mt']].copy()

# Sort or order by group if desired
ordered_groups = df['orig_ident'].unique()  # or sorted(df['orig_ident'].unique())

# Plot settings
n_per_row = 30
total_groups = len(ordered_groups)
n_rows = math.ceil(total_groups / n_per_row)

# Set figure size
fig_width = 16  # adjust as needed
fig_height = 4.3 * n_rows  # increase height per row
fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, fig_height), squeeze=False)

# Flatten axes for easier indexing
axes = axes.flatten()

# Plot in chunks
for i in range(n_rows):
    start = i * n_per_row
    end = min((i + 1) * n_per_row, total_groups)
    groups_subset = ordered_groups[start:end]

    ax = axes[i]
    sns.violinplot(
        data=df[df['orig_ident'].isin(groups_subset)],
        x='orig_ident',
        y='percent_mt',
        order=groups_subset,
        ax=ax,
        inner='quartile',
        scale='width'
    )
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel('percent_mt')
    ax.tick_params(axis='x', rotation=90)

plt.tight_layout()
plt.savefig('percent_mt_violinplot_by_orig_ident.png', dpi=900)
plt.show()


import numpy as np

# 1. Identify ribosomal genes
ribo_genes = [gene for gene in adata.var_names if gene.startswith('RPS') or gene.startswith('RPL')]

# 2. Calculate rRNA counts per cell
# If your matrix is sparse, use .X.A to convert to dense; or use .X directly if it's already dense
ribo_counts = np.array(adata[:, ribo_genes].X.sum(axis=1)).flatten()

# 3. Store in adata.obs
adata.obs['ribo_counts'] = ribo_counts

# (Optional) Calculate percentage of total RNA
total_counts = np.array(adata.X.sum(axis=1)).flatten()
adata.obs['percent_ribo'] = (ribo_counts / total_counts) * 100


import matplotlib.pyplot as plt
import seaborn as sns
import math
import pandas as pd

# Extract data
df = adata.obs[['orig_ident', 'percent_ribo']].copy()

# Sort or order by group if desired
ordered_groups = df['orig_ident'].unique()  # or sorted(df['orig_ident'].unique())

# Plot settings
n_per_row = 30
total_groups = len(ordered_groups)
n_rows = math.ceil(total_groups / n_per_row)

# Set figure size
fig_width = 16  # adjust as needed
fig_height = 4.3 * n_rows  # increase height per row
fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, fig_height), squeeze=False)

# Flatten axes for easier indexing
axes = axes.flatten()

# Plot in chunks
for i in range(n_rows):
    start = i * n_per_row
    end = min((i + 1) * n_per_row, total_groups)
    groups_subset = ordered_groups[start:end]

    ax = axes[i]
    sns.violinplot(
        data=df[df['orig_ident'].isin(groups_subset)],
        x='orig_ident',
        y='percent_ribo',
        order=groups_subset,
        ax=ax,
        inner='quartile',
        scale='width'
    )
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel('percent_ribo')
    ax.tick_params(axis='x', rotation=90)

plt.tight_layout()
plt.savefig('percent_ribo_violinplot_by_orig_ident.png', dpi=900)
plt.show()



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
adata

adata_subset=adata[adata.obs['tech']=='CosMx'].copy()

adata_subset.obs['immune_ME'].value_counts()

adata_subset=adata_subset[adata_subset.obs['immune_ME']!='Unknown'].copy()
adata_subset

adata_subset=adata_subset[adata_subset.obs['immune_cell_annotation_combined']=='B'].copy()
adata_subset

adata_subset.X=adata_subset.layers['counts'].copy()
#sc.pp.filter_genes(adata_subset, min_cells=3)
#sc.pp.normalize_total(adata_subset, target_sum=1e4)
#sc.pp.log1p(adata_subset)
#adata_subset

from sklearn.preprocessing import MinMaxScaler

# Step 1: Normalize and log-transform
sc.pp.normalize_total(adata_subset, target_sum=1e4)
sc.pp.log1p(adata_subset)

# Step 2: Extract the required genes (make sure they exist in the var_names)
required_genes = ['IGHM', 'IGHD', 'IGHA1', 'IGHG1', 'IGHG2']

# Step 3: Extract the expression values for these genes
gene_expr_df = pd.DataFrame(
    adata_subset[:, required_genes].X.toarray(),
    columns=required_genes,
    index=adata_subset.obs_names
)
print(gene_expr_df)
# Step 4: Scale the gene expression values (z-score)
scaler = MinMaxScaler()
gene_expr_scaled_df = pd.DataFrame(
    scaler.fit_transform(gene_expr_df),
    columns=gene_expr_df.columns,
    index=gene_expr_df.index
)

print(gene_expr_scaled_df)

gene_expr_scaled_df['sum_igg_iga']=gene_expr_scaled_df['IGHG2']+gene_expr_scaled_df['IGHG1']+gene_expr_scaled_df['IGHA1']
gene_expr_scaled_df['sum_igd']=gene_expr_scaled_df['IGHD']
gene_expr_scaled_df['sum_igm']=gene_expr_scaled_df['IGHM']

gene_expr_scaled_df.head(5)

def classify_igg_iga_igd_igm(row):
    if row['sum_igg_iga'] == 0 and row['sum_igd'] == 0 and row['sum_igm'] == 0:
        return 'negative'
    elif row['sum_igg_iga'] > row['sum_igd'] or row['sum_igg_iga'] > row['sum_igm']:
        return 'IgA/G_positive'
    elif row['sum_igd'] > row['sum_igg_iga'] or row['sum_igd'] > row['sum_igm']:
        return 'IgD_positive'
    elif row['sum_igm'] > row['sum_igg_iga'] or row['sum_igm'] > row['sum_igd']:
        return 'IgM_positive'
    else:
        return 'unclear'

# Apply the classification
gene_expr_scaled_df['class'] = gene_expr_scaled_df.apply(classify_igg_iga_igd_igm, axis=1)

# Optional: see distribution
print(gene_expr_scaled_df['class'].value_counts())


adata_subset.obs['class']=gene_expr_scaled_df['class'].copy()
print(adata_subset.obs['class'].value_counts())

print(adata_subset.obs['class'].value_counts())

adata_subset.obs['immune_ME'].value_counts()

adata_subset.obs['immune_ME']=adata_subset.obs['immune_ME'].replace('Immune Immune 1', 'T predom. Immune ME')

palette={
    'T predom. Immune ME': '#8d584e',
    'Fibro Immune ME': '#EE7402',
    'B predom. Immune ME': '#007ABA',
    'Inj. Tubular Immune ME': '#d279af',
    'Residential Immune ME': '#BBBE09',
    'Vascular Immune ME': '#00AFC5',
    'Glomerular Immune ME': '#E4051E',
}##'#92D2DF',

# Create the same summary dataframe
summary_df = (
    adata_subset.obs
    .groupby(['immune_ME', 'class'])
    .size()
    .unstack(fill_value=0)
)
summary_df['sum'] = summary_df[['IgA/G_positive', 'IgD_positive', 'IgM_positive', 'negative']].sum(axis=1)
summary_df['IgD_positive_fraction'] = summary_df['IgD_positive'] / summary_df['sum']
summary_df = summary_df[['IgD_positive_fraction']]  # keep only the relevant column
print(summary_df)
summary_df_reset = summary_df.reset_index()

# Create bar colors based on immune_subcluster_ME
bar_colors = summary_df_reset['immune_ME'].map(palette)
print(summary_df_reset)
fig, ax = plt.subplots(figsize=(4, 6))
ax.bar(
    summary_df_reset['immune_ME'],
    summary_df_reset['IgD_positive_fraction'],
    color=bar_colors
)

ax.set_ylabel("Fraction of IgD-positive B Cells")
ax.set_xlabel("")
ax.set_title("IgD+ fraction by Immune ME")
plt.xticks(rotation=90)
plt.tick_params(axis='x', which='both', length=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
# Save figure
plt.savefig("IgD_positive_fraction_by_immune_ME.png", dpi=900, bbox_inches='tight')
plt.show()


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sample_col = "unique_sample_identifier"

# 1) Per-sample, per-immune_ME IgD fraction
counts = (
    adata_subset.obs
    .groupby([sample_col, "immune_ME", "class"])
    .size()
    .unstack(fill_value=0)
)

# ensure columns exist
for c in ["IgA/G_positive", "IgD_positive", "IgM_positive", "negative"]:
    if c not in counts.columns:
        counts[c] = 0

counts["sum"] = counts[["IgA/G_positive", "IgD_positive", "IgM_positive", "negative"]].sum(axis=1)
counts["IgD_positive_fraction"] = counts["IgD_positive"] / counts["sum"].replace(0, np.nan)

points_df = (
    counts[["IgD_positive_fraction"]]
    .reset_index()
    .dropna(subset=["IgD_positive_fraction"])
)

# 2) Bar heights = mean across samples (or median if you prefer)
bar_df = (
    points_df.groupby("immune_ME", as_index=False)["IgD_positive_fraction"]
    .mean()
)

# consistent order
order = bar_df.sort_values("immune_ME")["immune_ME"].tolist()
bar_df["immune_ME"] = pd.Categorical(bar_df["immune_ME"], categories=order, ordered=True)
points_df["immune_ME"] = pd.Categorical(points_df["immune_ME"], categories=order, ordered=True)
bar_df = bar_df.sort_values("immune_ME")
points_df = points_df.sort_values("immune_ME")

# 3) Plot bars + overlaid sample points (with jitter)
fig, ax = plt.subplots(figsize=(4, 6))

bar_colors = bar_df["immune_ME"].map(palette)
ax.bar(bar_df["immune_ME"].astype(str), bar_df["IgD_positive_fraction"], color=bar_colors)

# jittered x positions
xpos = np.arange(len(order))
xmap = {me: i for i, me in enumerate(order)}
x_points = points_df["immune_ME"].map(xmap).astype(float).values
jitter = np.random.uniform(-0.15, 0.15, size=len(points_df))
ax.scatter(x_points + jitter, points_df["IgD_positive_fraction"].values,
           s=15, color="black", alpha=0.6, linewidths=0)

ax.set_xticks(xpos)
ax.set_xticklabels(order, rotation=90)
ax.set_ylabel("Fraction of IgD-positive B Cells")
ax.set_xlabel("")
ax.set_title("IgD+ fraction by Immune ME")

ax.tick_params(axis="x", which="both", length=0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.ylim(0,0.2)
plt.tight_layout()
plt.savefig("IgD_positive_fraction_by_immune_ME_individual_dots_sample.png", dpi=900, bbox_inches="tight")
plt.show()


import scipy.stats as stats
adata_subset.obs['is_IgD_positive'] = adata_subset.obs['class'] == 'IgD_positive'

b_cell_me = adata_subset[adata_subset.obs['immune_ME'] == 'B predom. Immune ME']
other_me = adata_subset[adata_subset.obs['immune_ME'] != 'B predom. Immune ME']

stat, pval = stats.mannwhitneyu(
    b_cell_me.obs['is_IgD_positive'].astype(int),
    other_me.obs['is_IgD_positive'].astype(int),
    alternative='two-sided'
)
print(f"Mann–Whitney U test: U={stat:.4f}, p-value={pval:.4e}")


adata_subset.obs['class'].value_counts()

import scipy.stats as stats
adata_subset.obs['is_IgM_positive'] = adata_subset.obs['class'] == 'IgM_positive'

b_cell_me = adata_subset[adata_subset.obs['immune_ME'] == 'B predom. Immune ME']
other_me = adata_subset[adata_subset.obs['immune_ME'] != 'B predom. Immune ME']

stat, pval = stats.mannwhitneyu(
    b_cell_me.obs['is_IgM_positive'].astype(int),
    other_me.obs['is_IgM_positive'].astype(int),
    alternative='two-sided'
)
print(f"Mann–Whitney U test: U={stat:.4f}, p-value={pval:.4e}")


import scipy.stats as stats
adata_subset.obs['is_IgA_G_positive'] = adata_subset.obs['class'] == 'IgA/G_positive	'

b_cell_me = adata_subset[adata_subset.obs['immune_ME'] == 'B predom. Immune ME']
other_me = adata_subset[adata_subset.obs['immune_ME'] != 'B predom. Immune ME']

stat, pval = stats.mannwhitneyu(
    b_cell_me.obs['is_IgA_G_positive'].astype(int),
    other_me.obs['is_IgA_G_positive'].astype(int),
    alternative='two-sided'
)
print(f"Mann–Whitney U test: U={stat:.4f}, p-value={pval:.4e}")


# Create the same summary dataframe
summary_df = (
    adata_subset.obs
    .groupby(['immune_ME', 'class'])
    .size()
    .unstack(fill_value=0)
)
summary_df['sum'] = summary_df[['IgA/G_positive', 'IgD_positive', 'IgM_positive', 'negative']].sum(axis=1)
summary_df['IgM_positive_fraction'] = summary_df['IgM_positive'] / summary_df['sum']
summary_df = summary_df[['IgM_positive_fraction']]  # keep only the relevant column
print(summary_df)
summary_df_reset = summary_df.reset_index()

# Create bar colors based on immune_subcluster_ME
bar_colors = summary_df_reset['immune_ME'].map(palette)
print(summary_df_reset)
fig, ax = plt.subplots(figsize=(4, 6))
ax.bar(
    summary_df_reset['immune_ME'],
    summary_df_reset['IgM_positive_fraction'],
    color=bar_colors
)

ax.set_ylabel("Fraction of IgM-positive B Cells")
ax.set_xlabel("")
ax.set_title("IgM+ fraction by Immune ME")
plt.xticks(rotation=90)
plt.tick_params(axis='x', which='both', length=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
# Save figure
plt.savefig("IgM_positive_fraction_by_immune_ME.png", dpi=900, bbox_inches='tight')
plt.show()


# Create the same summary dataframe
summary_df = (
    adata_subset.obs
    .groupby(['immune_ME', 'class'])
    .size()
    .unstack(fill_value=0)
)
summary_df['sum'] = summary_df[['IgA/G_positive', 'IgD_positive', 'IgM_positive', 'negative']].sum(axis=1)
summary_df['IgA/G_positive_fraction'] = summary_df['IgA/G_positive'] / summary_df['sum']
summary_df = summary_df[['IgA/G_positive_fraction']]  # keep only the relevant column
print(summary_df)
summary_df_reset = summary_df.reset_index()

# Create bar colors based on immune_subcluster_ME
bar_colors = summary_df_reset['immune_ME'].map(palette)
print(summary_df_reset)
fig, ax = plt.subplots(figsize=(4, 6))
ax.bar(
    summary_df_reset['immune_ME'],
    summary_df_reset['IgA/G_positive_fraction'],
    color=bar_colors
)

ax.set_ylabel("Fraction of IgG/IgA-positive B Cells")
ax.set_xlabel("")
ax.set_title("IgG/IgA fraction by Immune ME")
plt.xticks(rotation=90)
plt.tick_params(axis='x', which='both', length=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
# Save figure
plt.savefig("IgA_G_positive_fraction_by_immune_ME.png", dpi=900, bbox_inches='tight')
plt.show()

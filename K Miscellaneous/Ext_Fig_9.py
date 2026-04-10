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

adata_subset=adata[adata.obs['immune_ME']=='Residential Immune ME'].copy()
adata_subset

samples = adata_subset.obs['unique_sample_identifier']
cluster_labels = adata_subset.obs['immune_cell_annotation_combined']

order= [
    '1072_CosMx', '1009_Xenium', '1071_CosMx', '1068_CosMx',  '1005_Xenium',
    '1013_Xenium', '1011_Xenium', '1012_Xenium','1003_Xenium', '1006_Xenium', '1062_CosMx',
    '1008_Xenium', 'HK2844_2_CosMx', 'HK2844_CosMx','1004_Xenium', '1010_Xenium',
    '1064_CosMx', 'HK3069_CosMx', '1001_Xenium', 'HK3421_CosMx', '1007_Xenium',
    'HK3035_2_CosMx', 'HK3035_CosMx', 'HK3606_CosMx', 'HK2873_CosMx', 'HK2841_CosMx',
    'HK3474_CosMx', 'HK2695_Xenium', 'HK3631_CosMx', 'HK2695_CosMx', 'HK3066_CosMx',
    'HK3535_CosMx', 'HK3594_CosMx', 'HK3614_CosMx', 'HK3070_CosMx', 'HK2874_CosMx',
    'HK3068_CosMx', 'HK3591_CosMx', 'HK3588_CosMx', 'HK3469_CosMx', 'HK3542_CosMx',
    'HK3039_CosMx', 'HK2924_CosMx', 'HK2753_CosMx', 'HK2753_Xenium', '1061_CosMx',
    'HK3106_Xenium', 'HK3106_CosMx', 'HK3647_CosMx', 'HK3624_CosMx', 'HK3531_CosMx',
    'HK3531_2_CosMx', 'HK3063_CosMx', 'HK2989_CosMx', 'HK3623_CosMx', '1063_CosMx',
    'HK3616_CosMx', 'HK3626_Xenium', 'HK3626_CosMx', 'HK2990_CosMx', 'HK3612_CosMx',
     'HK3565_CosMx', 'Pediatric1_CosMx'
]

# Create a DataFrame from the relevant columns
data = pd.DataFrame({'unique_sample_identifier': samples, 'immune_cell_annotation_combined': cluster_labels})

# Group by 'sample' and 'niches_annotation_based' and count the occurrences
grouped_data = data.groupby(['unique_sample_identifier', 'immune_cell_annotation_combined']).size().unstack(fill_value=0)

# Get sample totals in the same order as grouped_data's index
sample_totals = data['unique_sample_identifier'].value_counts()
sample_totals = sample_totals.loc[grouped_data.index]  # Match the index order

# Normalize the counts
normalized_count_df = grouped_data.div(sample_totals, axis='index')
ordered_df = normalized_count_df.loc[order]
ordered_df.index = ordered_df.index.str.replace('Xenium', 'X').str.replace('CosMx', 'C')
ax = ordered_df.plot(kind='bar', stacked=True, figsize=(16, 3.8))
ax.grid(False)
# Set plot labels and title
plt.xlabel('')
plt.ylabel('Relative Amount')
plt.title('')
plt.legend(title='ME Labels', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
plt.tick_params(axis='x', which='both', length=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('Kindey_residential_immune_ME.png', dpi=900, bbox_inches='tight')
plt.show()



df = adata.obs[['unique_sample_identifier', 'immune_ME']].copy()

df['immune_ME']=df['immune_ME'].replace('Unknown', 'other')
df['immune_ME']=df['immune_ME'].replace('Immune Immune 1', 'other')
df['immune_ME']=df['immune_ME'].replace('Inj. Tubular Immune ME', 'other')
df['immune_ME']=df['immune_ME'].replace('Fibro Immune ME', 'other')
df['immune_ME']=df['immune_ME'].replace('Vascular Immune ME', 'other')
df['immune_ME']=df['immune_ME'].replace('Glomerular Immune ME', 'other')
df['immune_ME']=df['immune_ME'].replace('B predom. Immune ME', 'other')
df['immune_ME']=df['immune_ME'].replace('Residential Immune ME', 'residential')
df

ordered_df

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Group and count per sample and immune_ME type
counts = df.groupby(['unique_sample_identifier', 'immune_ME']).size().unstack(fill_value=0)

# Calculate fraction of residential cells
counts['res_fraction'] = counts['residential'] / counts.sum(axis=1)

# Reset for plotting
plot_df = counts[['res_fraction']].reset_index()
plot_df=plot_df[plot_df['unique_sample_identifier']!='1070_CosMx']
plot_df=plot_df.set_index('unique_sample_identifier')
ordered_df = plot_df.loc[order]
ordered_df.index = ordered_df.index.str.replace('Xenium', 'X').str.replace('CosMx', 'C')
# Plot
plt.figure(figsize=(10, 3))
sns.barplot(data=ordered_df, x='unique_sample_identifier', y='res_fraction', color='skyblue')
plt.xticks(rotation=90)
plt.ylabel('Fraction')
plt.xlabel('')
plt.title('Residential Immune Cells / Total Cell Fraction')
plt.tight_layout()
plt.savefig('Kindey_residential_immune_ME_fractions.png', dpi=900, bbox_inches='tight')
plt.show()


import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

# --- 1) rebuild the sample-level table ---
df = adata.obs[['unique_sample_identifier','immune_ME','Condition']].copy()
df['immune_ME'] = (df['immune_ME']
                   .replace(['Unknown','Immune Immune 1','Inj. Tubular Immune ME',
                             'Fibro Immune ME','Vascular Immune ME',
                             'Glomerular Immune ME','B predom. Immune ME'], 'other')
                   .replace('Residential Immune ME','residential'))

# --- 2) per-sample counts ---
counts = df.groupby(['unique_sample_identifier','immune_ME']).size().unstack(fill_value=0)

# --- 3) total cells per sample ---
counts['total_cells'] = counts.sum(axis=1)

# --- 4) fraction of residential cells ---
counts['res_fraction'] = counts.get('residential',0) / counts['total_cells']

# --- 5) keep only samples with ≥300 residential cells ---
counts = counts[counts.get('residential',0) >= 300].copy()

# --- 6) add Condition ---
cond_map = df.drop_duplicates('unique_sample_identifier') \
             .set_index('unique_sample_identifier')['Condition']
counts['Condition'] = counts.index.map(cond_map)

# --- 7) keep only Control & DKD ---
counts = counts[counts['Condition'].isin(['Control','DKD'])].copy()

# --- 8) compare fractions ---
ctrl = counts.loc[counts['Condition']=='Control','res_fraction']
dkd  = counts.loc[counts['Condition']=='DKD','res_fraction']

stat, pval = mannwhitneyu(ctrl, dkd, alternative='two-sided')
print(f"Mann–Whitney U test: p = {pval:.3g}, n_Control = {len(ctrl)}, n_DKD = {len(dkd)}")

# --- 9) optional boxplot ---
import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(3.2,3.2))
sns.boxplot(data=counts, x='Condition', y='res_fraction', palette='Set2', order=['Control','DKD'])
sns.stripplot(data=counts, x='Condition', y='res_fraction', color='k', size=4,
              alpha=0.6, jitter=True, order=['Control','DKD'])
plt.ylabel('Fraction of residential immune ME')
plt.title(f'p = {pval:.3g}')
plt.tight_layout()
plt.savefig('Residential_Immune_ME_fraction_Control_vs_DKD.png', dpi=900, bbox_inches='tight')
plt.show()

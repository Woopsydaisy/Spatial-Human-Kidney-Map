import pandas as pd

##glomerular niche
correlations_glom_dkd=pd.read_csv('/content/hotspot_local_correlations_glom_niche_dkd.csv')
correlations_glom_dkd.set_index('Unnamed: 0',inplace=True)
correlations_glom_dkd.index.name=None
#correlations_glom_dkd

correlations_glom_control=pd.read_csv('/content/hotspot_local_correlations_glomerular_niche_control.csv')
correlations_glom_control.set_index('Unnamed: 0',inplace=True)
correlations_glom_control.index.name=None
#correlations_glom_control

# Step 2: Find the common genes (row and column labels must match)
common_genes = correlations_glom_control.index.intersection(correlations_glom_dkd.index)

# Step 3: Subset both matrices to the common genes and ensure the same order
control_subset = correlations_glom_control.loc[common_genes, common_genes]
disease_subset = correlations_glom_dkd.loc[common_genes, common_genes]

correlation_difference = disease_subset - control_subset

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.clustermap(
    scaled_correlation_control,
    cmap="coolwarm",
    vmin=-0.5,
    vmax=0.5,
    center=0,
    figsize=(20, 20),
    linewidths=0,
    cbar_kws={"label": "Δ correlation (disease - control)"},
    dendrogram_ratio=(0.001, 0.001),
    tree_kws={"linewidths": 0},
    cbar_pos=None
)

#plt.suptitle("Clustered Gene–Gene Correlation Differences (Disease vs. Control)", y=1.02)
#plt.savefig("clustered_correlation_difference_heatmap.png", dpi=300, bbox_inches='tight')
plt.close()
#plt.show()


from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(control_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_control = pd.DataFrame(scaled_array, index=control_subset.index, columns=control_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(disease_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_dkd = pd.DataFrame(scaled_array, index=disease_subset.index, columns=disease_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(correlation_difference.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_difference = pd.DataFrame(scaled_array, index=correlation_difference.index, columns=correlation_difference.columns)

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_dkd.index[row_order]
col_labels = scaled_correlation_dkd.columns[col_order]

# Reorder DKD matrix
scaled_correlation_dkd_ordered = scaled_correlation_dkd.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_control.index[row_order]
col_labels = scaled_correlation_control.columns[col_order]

# Reorder DKD matrix
scaled_correlation_control_ordered = scaled_correlation_control.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_difference.index[row_order]
col_labels = scaled_correlation_difference.columns[col_order]

# Reorder DKD matrix
scaled_correlation_difference_ordered = scaled_correlation_difference.loc[row_labels, col_labels]

##correlated in control and dkd
import numpy as np
threshold = 0.99
high_corr_disease = scaled_correlation_dkd.abs() > threshold
high_corr_control = scaled_correlation_control.abs() > threshold
high_in_both = high_corr_disease & high_corr_control
high_pairs = np.argwhere(high_in_both.values)
gene_pairs = [(high_in_both.index[i], high_in_both.columns[j]) for i, j in high_pairs if i < j]  # avoid duplicates
len(gene_pairs)

gene_pairs

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
      ('NPHS2', 'ITGA3'),
      ('NPHS2', 'AIF1'),
    ('NPHS2', 'ANXA2'),
    ('IGFBP7', 'TEK'),
    ('IGFBP7', 'KDR'),
    ('IGFBP7', 'TIE1'),
    ('IGFBP7', 'CDH5'),
    ('ENG', 'ITGA8'),
    ('PECAM1', 'ADGRF5'),
      ('ITGB6', 'NOS1'),
      ('KLF2', 'MEIS2'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_control_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_control_ordered,
    cmap="coolwarm",
    vmin=-0.5, vmax=0.5, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_control_glomerular_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
     ('NPHS2', 'ITGA3'),
      ('NPHS2', 'AIF1'),
    ('NPHS2', 'ANXA2'),
    ('IGFBP7', 'TEK'),
    ('IGFBP7', 'KDR'),
    ('IGFBP7', 'TIE1'),
    ('IGFBP7', 'CDH5'),
    ('ENG', 'ITGA8'),
    ('PECAM1', 'ADGRF5'),
      ('ITGB6', 'NOS1'),
      ('KLF2', 'MEIS2'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_dkd_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_dkd_ordered,
    cmap="coolwarm",
    vmin=-0.5, vmax=0.5, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_dkd_glomerular_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


# Flatten the matrix
flat_diff = correlation_difference.stack().reset_index()
flat_diff.columns = ['Gene1', 'Gene2', 'CorrelationDifference']

# Remove self-correlations (diagonal)
flat_diff = flat_diff[flat_diff['Gene1'] != flat_diff['Gene2']]

# Drop duplicates for symmetrical pairs (e.g., A–B and B–A)
flat_diff['SortedPair'] = flat_diff.apply(lambda row: tuple(sorted([row['Gene1'], row['Gene2']])), axis=1)
flat_diff = flat_diff.drop_duplicates('SortedPair').drop(columns='SortedPair')

# Get top N positive and negative differences
top_n = 30
top_positive = flat_diff.sort_values(by='CorrelationDifference', ascending=False).head(top_n)
top_negative = flat_diff.sort_values(by='CorrelationDifference', ascending=True).head(top_n)

print("Top positive correlation differences:")
print(top_positive)

print("\nTop negative correlation differences:")
print(top_negative)


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('NPHS2', 'FOS'),
    ('PTGDS', 'FOS'),
    ('COL4A1', 'MT2A'),
    ('COL3A1', 'MT2A'),
    ('COL3A1', 'ITGB8'),
    ('SPOCK2', 'CD24'),
    ('ANXA2', 'RBM47'),
    ('CHI3L1', 'JUN'),
    ('VIM', 'RBM47'),
    ('COL4A2', 'MT2A'),
    ('IGFBP5', 'CD24')
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_difference_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_difference_ordered,
    cmap="coolwarm",
    vmin=-0.5, vmax=0.5, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_differences_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()








correlations_iPT_dkd=pd.read_csv('/content/hotspot_local_correlations_iPT_niche_dkd.csv')
correlations_iPT_dkd.set_index('Unnamed: 0',inplace=True)
correlations_iPT_dkd.index.name=None
#correlations_iPT_dkd

correlations_iPT_control=pd.read_csv('/content/hotspot_local_correlations_iPT_niche_control.csv')
correlations_iPT_control.set_index('Unnamed: 0',inplace=True)
correlations_iPT_control.index.name=None
#correlations_iPT_control

# Step 2: Find the common genes (row and column labels must match)
common_genes = correlations_iPT_control.index.intersection(correlations_iPT_dkd.index)

# Step 3: Subset both matrices to the common genes and ensure the same order
control_subset = correlations_iPT_control.loc[common_genes, common_genes]
disease_subset = correlations_iPT_dkd.loc[common_genes, common_genes]

correlation_difference = disease_subset - control_subset

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(control_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_control = pd.DataFrame(scaled_array, index=control_subset.index, columns=control_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(disease_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_dkd = pd.DataFrame(scaled_array, index=disease_subset.index, columns=disease_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(correlation_difference.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_difference = pd.DataFrame(scaled_array, index=correlation_difference.index, columns=correlation_difference.columns)

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.clustermap(
    scaled_correlation_control,
    cmap="coolwarm",
    vmin=-0.5,
    vmax=0.5,
    center=0,
    figsize=(20, 20),
    linewidths=0,
    cbar_kws={"label": "Δ correlation (disease - control)"},
    dendrogram_ratio=(0.001, 0.001),
    tree_kws={"linewidths": 0},
    cbar_pos=None
)

#plt.suptitle("Clustered Gene–Gene Correlation Differences (Disease vs. Control)", y=1.02)
#plt.savefig("clustered_correlation_difference_heatmap.png", dpi=300, bbox_inches='tight')
plt.close()
#plt.show()


# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_dkd.index[row_order]
col_labels = scaled_correlation_dkd.columns[col_order]

# Reorder DKD matrix
scaled_correlation_dkd_ordered = scaled_correlation_dkd.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_control.index[row_order]
col_labels = scaled_correlation_control.columns[col_order]

# Reorder DKD matrix
scaled_correlation_control_ordered = scaled_correlation_control.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_difference.index[row_order]
col_labels = scaled_correlation_difference.columns[col_order]

# Reorder DKD matrix
scaled_correlation_difference_ordered = scaled_correlation_difference.loc[row_labels, col_labels]

##correlated in control and dkd
import numpy as np
threshold = 0.99
high_corr_disease = scaled_correlation_dkd.abs() > threshold
high_corr_control = scaled_correlation_control.abs() > threshold
high_in_both = high_corr_disease & high_corr_control
high_pairs = np.argwhere(high_in_both.values)
gene_pairs = [(high_in_both.index[i], high_in_both.columns[j]) for i, j in high_pairs if i < j]  # avoid duplicates
len(gene_pairs)

gene_pairs

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
      ('IFITM3', 'B2M'),
      ('AZGP1', 'ACE'),
    ('CUBN', 'AZGP1'),
('IL32', 'ACSM2B'),
 ('IL32', 'TACSTD2'),
 ('IL32', 'COTL1'),
 ('ACSM2B', 'APOE'),
 ('ACSM2B', 'MAF'),
 ('ACSM2B', 'LRP2'),
 ('ACSM2B', 'GPX3'),
      ('MAF', 'TNFRSF12A'),
      ('CRYAB', 'SERPINA1'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_control_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_control_ordered,
    cmap="coolwarm",
    vmin=-0.5, vmax=0.5, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_control_iPT_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
     ('IFITM3', 'B2M'),
      ('AZGP1', 'ACE'),
    ('CUBN', 'AZGP1'),
('IL32', 'ACSM2B'),
 ('IL32', 'TACSTD2'),
 ('IL32', 'COTL1'),
 ('ACSM2B', 'APOE'),
 ('ACSM2B', 'MAF'),
 ('ACSM2B', 'LRP2'),
 ('ACSM2B', 'GPX3'),
      ('MAF', 'TNFRSF12A'),
      ('CRYAB', 'SERPINA1'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_dkd_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_dkd_ordered,
    cmap="coolwarm",
    vmin=-0.5, vmax=0.5, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_dkd_iPT_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


# Flatten the matrix
flat_diff = correlation_difference.stack().reset_index()
flat_diff.columns = ['Gene1', 'Gene2', 'CorrelationDifference']

# Remove self-correlations (diagonal)
flat_diff = flat_diff[flat_diff['Gene1'] != flat_diff['Gene2']]

# Drop duplicates for symmetrical pairs (e.g., A–B and B–A)
flat_diff['SortedPair'] = flat_diff.apply(lambda row: tuple(sorted([row['Gene1'], row['Gene2']])), axis=1)
flat_diff = flat_diff.drop_duplicates('SortedPair').drop(columns='SortedPair')

# Get top N positive and negative differences
top_n = 30
top_positive = flat_diff.sort_values(by='CorrelationDifference', ascending=False).head(top_n)
top_negative = flat_diff.sort_values(by='CorrelationDifference', ascending=True).head(top_n)

print("Top positive correlation differences:")
print(top_positive)

print("\nTop negative correlation differences:")
print(top_negative)


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('IL32', 'CRYAB'),
    ('ACSM2B', 'CRYAB'),
    ('APOE', 'CRYAB'),
    ('LRP2', 'CRYAB'),
    ('IGHG1', 'IGHG2'),
    ('IGHG1', 'IGKC'),
    ('IGKC', 'IGHG2'),
    ('IL32', 'HAVCR1'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_difference_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_difference_ordered,
    cmap="coolwarm",
    vmin=-0.5, vmax=0.5, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_differences_iPT_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()










correlations_PT_dkd=pd.read_csv('/content/hotspot_local_correlations_pt_niche_dkd.csv')
correlations_PT_dkd.set_index('Unnamed: 0',inplace=True)
correlations_PT_dkd.index.name=None
#correlations_PT_dkd

correlations_PT_control=pd.read_csv('/content/hotspot_local_correlations_pt_niche_control.csv')
correlations_PT_control.set_index('Unnamed: 0',inplace=True)
correlations_PT_control.index.name=None
#correlations_PT_control

# Step 2: Find the common genes (row and column labels must match)
common_genes = correlations_PT_control.index.intersection(correlations_PT_dkd.index)

# Step 3: Subset both matrices to the common genes and ensure the same order
control_subset = correlations_PT_control.loc[common_genes, common_genes]
disease_subset = correlations_PT_dkd.loc[common_genes, common_genes]

correlation_difference = disease_subset - control_subset

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(control_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_control = pd.DataFrame(scaled_array, index=control_subset.index, columns=control_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(disease_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_dkd = pd.DataFrame(scaled_array, index=disease_subset.index, columns=disease_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(correlation_difference.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_difference = pd.DataFrame(scaled_array, index=correlation_difference.index, columns=correlation_difference.columns)

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.clustermap(
    scaled_correlation_control,
    cmap="coolwarm",
    vmin=-0.5,
    vmax=0.5,
    center=0,
    figsize=(20, 20),
    linewidths=0,
    cbar_kws={"label": "Δ correlation (disease - control)"},
    dendrogram_ratio=(0.001, 0.001),
    tree_kws={"linewidths": 0},
    cbar_pos=None
)

#plt.suptitle("Clustered Gene–Gene Correlation Differences (Disease vs. Control)", y=1.02)
#plt.savefig("clustered_correlation_difference_heatmap.png", dpi=300, bbox_inches='tight')
plt.close()
#plt.show()


# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_dkd.index[row_order]
col_labels = scaled_correlation_dkd.columns[col_order]

# Reorder DKD matrix
scaled_correlation_dkd_ordered = scaled_correlation_dkd.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_control.index[row_order]
col_labels = scaled_correlation_control.columns[col_order]

# Reorder DKD matrix
scaled_correlation_control_ordered = scaled_correlation_control.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_difference.index[row_order]
col_labels = scaled_correlation_difference.columns[col_order]

# Reorder DKD matrix
scaled_correlation_difference_ordered = scaled_correlation_difference.loc[row_labels, col_labels]

##correlated in control and dkd
import numpy as np
threshold = 0.99
high_corr_disease = scaled_correlation_dkd.abs() > threshold
high_corr_control = scaled_correlation_control.abs() > threshold
high_in_both = high_corr_disease & high_corr_control
high_pairs = np.argwhere(high_in_both.values)
gene_pairs = [(high_in_both.index[i], high_in_both.columns[j]) for i, j in high_pairs if i < j]  # avoid duplicates
len(gene_pairs)

gene_pairs

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('APOE', 'SLC13A3'),
    ('S100A6', 'CUBN'),
    ('ACSM2B', 'CUBN'),
    ('ACSM2B', 'CD24'),
    ('GPX3', 'SLC13A3'),
    ('APOE', 'ITM2A'),
    ('COTL1', 'FGF9'),
    ('S100A6', 'IGFBP7'),
    ('GPX3', 'IGFBP3'),
    ('GPX3', 'NEAT1'),
    ('GPX3', 'COTL1'),
    ('ACSM2B', 'SIGIRR'),
     ('GPX3', 'SLC4A4'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_control_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_control_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_control_PT_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('APOE', 'SLC13A3'),
    ('S100A6', 'CUBN'),
    ('ACSM2B', 'CUBN'),
    ('ACSM2B', 'CD24'),
    ('GPX3', 'SLC13A3'),
    ('APOE', 'ITM2A'),
    ('COTL1', 'FGF9'),
    ('S100A6', 'IGFBP7'),
    ('GPX3', 'IGFBP3'),
    ('GPX3', 'NEAT1'),
    ('GPX3', 'COTL1'),
    ('ACSM2B', 'SIGIRR'),
     ('GPX3', 'SLC4A4'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_dkd_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_dkd_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_dkd_PT_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


# Flatten the matrix
flat_diff = correlation_difference.stack().reset_index()
flat_diff.columns = ['Gene1', 'Gene2', 'CorrelationDifference']

# Remove self-correlations (diagonal)
flat_diff = flat_diff[flat_diff['Gene1'] != flat_diff['Gene2']]

# Drop duplicates for symmetrical pairs (e.g., A–B and B–A)
flat_diff['SortedPair'] = flat_diff.apply(lambda row: tuple(sorted([row['Gene1'], row['Gene2']])), axis=1)
flat_diff = flat_diff.drop_duplicates('SortedPair').drop(columns='SortedPair')

# Get top N positive and negative differences
top_n = 30
top_positive = flat_diff.sort_values(by='CorrelationDifference', ascending=False).head(top_n)
top_negative = flat_diff.sort_values(by='CorrelationDifference', ascending=True).head(top_n)

print("Top positive correlation differences:")
print(top_positive)

print("\nTop negative correlation differences:")
print(top_negative)


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('DUSP1', 'MT2A'),
    ('COL3A1', 'MT2A'),
    ('MT2A', 'IGFBP5'),
    ('FOS', 'SPP1'),
    ('COL3A1', 'SPP1'),
    ('MT1X', 'IGFBP7'),
    ('MT1X', 'JUN'),
    ('SPP1', 'IGFBP5'),
    ('SPP1', 'B2M')
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_difference_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_difference_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_differences_PT_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()












correlations_PC_dkd=pd.read_csv('/content/hotspot_local_correlations_pc_niche_dkd.csv')
correlations_PC_dkd.set_index('Unnamed: 0',inplace=True)
correlations_PC_dkd.index.name=None
#correlations_PC_dkd

correlations_PC_control=pd.read_csv('/content/hotspot_local_correlations_pc_niche_control.csv')
correlations_PC_control.set_index('Unnamed: 0',inplace=True)
correlations_PC_control.index.name=None
#correlations_PC_control

# Step 2: Find the common genes (row and column labels must match)
common_genes = correlations_PC_control.index.intersection(correlations_PC_dkd.index)

# Step 3: Subset both matrices to the common genes and ensure the same order
control_subset = correlations_PC_control.loc[common_genes, common_genes]
disease_subset = correlations_PC_dkd.loc[common_genes, common_genes]

correlation_difference = disease_subset - control_subset

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(control_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_control = pd.DataFrame(scaled_array, index=control_subset.index, columns=control_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(disease_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_dkd = pd.DataFrame(scaled_array, index=disease_subset.index, columns=disease_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(correlation_difference.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_difference = pd.DataFrame(scaled_array, index=correlation_difference.index, columns=correlation_difference.columns)

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.clustermap(
    scaled_correlation_control,
    cmap="coolwarm",
    vmin=-0.5,
    vmax=0.5,
    center=0,
    figsize=(20, 20),
    linewidths=0,
    cbar_kws={"label": "Δ correlation (disease - control)"},
    dendrogram_ratio=(0.001, 0.001),
    tree_kws={"linewidths": 0},
    cbar_pos=None
)

#plt.suptitle("Clustered Gene–Gene Correlation Differences (Disease vs. Control)", y=1.02)
#plt.savefig("clustered_correlation_difference_heatmap.png", dpi=300, bbox_inches='tight')
plt.close()
#plt.show()


# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_dkd.index[row_order]
col_labels = scaled_correlation_dkd.columns[col_order]

# Reorder DKD matrix
scaled_correlation_dkd_ordered = scaled_correlation_dkd.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_control.index[row_order]
col_labels = scaled_correlation_control.columns[col_order]

# Reorder DKD matrix
scaled_correlation_control_ordered = scaled_correlation_control.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_difference.index[row_order]
col_labels = scaled_correlation_difference.columns[col_order]

# Reorder DKD matrix
scaled_correlation_difference_ordered = scaled_correlation_difference.loc[row_labels, col_labels]

##correlated in control and dkd
import numpy as np
threshold = 0.99
high_corr_disease = scaled_correlation_dkd.abs() > threshold
high_corr_control = scaled_correlation_control.abs() > threshold
high_in_both = high_corr_disease & high_corr_control
high_pairs = np.argwhere(high_in_both.values)
gene_pairs = [(high_in_both.index[i], high_in_both.columns[j]) for i, j in high_pairs if i < j]  # avoid duplicates
len(gene_pairs)

gene_pairs

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('AQP2', 'AQP3'),
    ('CALB1', 'VEGFA'),
    ('AQP2', 'CXCL14'),
    ('SPP1', 'GATA3'),
    ('CXCL14', 'ADIRF'),
    ('SLC12A3', 'TRPM6')
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_control_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_control_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_control_PC_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
      ('AQP2', 'AQP3'),
    ('CALB1', 'VEGFA'),
    ('AQP2', 'CXCL14'),
    ('SPP1', 'GATA3'),
    ('CXCL14', 'ADIRF'),
    ('SLC12A3', 'TRPM6')
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_dkd_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_dkd_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_dkd_PC_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


# Flatten the matrix
flat_diff = correlation_difference.stack().reset_index()
flat_diff.columns = ['Gene1', 'Gene2', 'CorrelationDifference']

# Remove self-correlations (diagonal)
flat_diff = flat_diff[flat_diff['Gene1'] != flat_diff['Gene2']]

# Drop duplicates for symmetrical pairs (e.g., A–B and B–A)
flat_diff['SortedPair'] = flat_diff.apply(lambda row: tuple(sorted([row['Gene1'], row['Gene2']])), axis=1)
flat_diff = flat_diff.drop_duplicates('SortedPair').drop(columns='SortedPair')

# Get top N positive and negative differences
top_n = 30
top_positive = flat_diff.sort_values(by='CorrelationDifference', ascending=False).head(top_n)
top_negative = flat_diff.sort_values(by='CorrelationDifference', ascending=True).head(top_n)

print("Top positive correlation differences:")
print(top_positive)

print("\nTop negative correlation differences:")
print(top_negative)


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('CALB1', 'HSPA1A/B'),
    ('VEGFA', 'HSPA1A/B'),
    ('VEGFA', 'HSPB1'),
    ('DUSP1', 'KRT19'),
    ('KRT7', 'CRYAB'),
    ('FOS', 'S100A2')
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_difference_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_difference_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_differences_PC_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()










correlations_LOH_dkd=pd.read_csv('/content/hotspot_local_correlations_loh_niche_dkd.csv')
correlations_LOH_dkd.set_index('Unnamed: 0',inplace=True)
correlations_LOH_dkd.index.name=None
#correlations_LOH_dkd

correlations_LOH_control=pd.read_csv('/content/hotspot_local_correlations_loh_niche_control.csv')
correlations_LOH_control.set_index('Unnamed: 0',inplace=True)
correlations_LOH_control.index.name=None
#correlations_LOH_control

# Step 2: Find the common genes (row and column labels must match)
common_genes = correlations_LOH_control.index.intersection(correlations_LOH_dkd.index)

# Step 3: Subset both matrices to the common genes and ensure the same order
control_subset = correlations_LOH_control.loc[common_genes, common_genes]
disease_subset = correlations_LOH_dkd.loc[common_genes, common_genes]

correlation_difference = disease_subset - control_subset

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(control_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_control = pd.DataFrame(scaled_array, index=control_subset.index, columns=control_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(disease_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_dkd = pd.DataFrame(scaled_array, index=disease_subset.index, columns=disease_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(correlation_difference.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_difference = pd.DataFrame(scaled_array, index=correlation_difference.index, columns=correlation_difference.columns)

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.clustermap(
    scaled_correlation_control,
    cmap="coolwarm",
    vmin=-0.5,
    vmax=0.5,
    center=0,
    figsize=(20, 20),
    linewidths=0,
    cbar_kws={"label": "Δ correlation (disease - control)"},
    dendrogram_ratio=(0.001, 0.001),
    tree_kws={"linewidths": 0},
    cbar_pos=None
)

#plt.suptitle("Clustered Gene–Gene Correlation Differences (Disease vs. Control)", y=1.02)
#plt.savefig("clustered_correlation_difference_heatmap.png", dpi=300, bbox_inches='tight')
plt.close()
#plt.show()


# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_dkd.index[row_order]
col_labels = scaled_correlation_dkd.columns[col_order]

# Reorder DKD matrix
scaled_correlation_dkd_ordered = scaled_correlation_dkd.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_control.index[row_order]
col_labels = scaled_correlation_control.columns[col_order]

# Reorder DKD matrix
scaled_correlation_control_ordered = scaled_correlation_control.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_difference.index[row_order]
col_labels = scaled_correlation_difference.columns[col_order]

# Reorder DKD matrix
scaled_correlation_difference_ordered = scaled_correlation_difference.loc[row_labels, col_labels]

##correlated in control and dkd
import numpy as np
threshold = 0.99
high_corr_disease = scaled_correlation_dkd.abs() > threshold
high_corr_control = scaled_correlation_control.abs() > threshold
high_in_both = high_corr_disease & high_corr_control
high_pairs = np.argwhere(high_in_both.values)
gene_pairs = [(high_in_both.index[i], high_in_both.columns[j]) for i, j in high_pairs if i < j]  # avoid duplicates
len(gene_pairs)

gene_pairs

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('UMOD', 'SLC12A1'),
 ('UMOD', 'AQP3'),
 ('UMOD', 'TPM1'),
 ('UMOD', 'CALD1'),
 ('UMOD', 'JAK1'),
    ('APOE', 'EZR'),
 ('SLC12A1', 'EGF'),
    ('FOS', 'JUNB'),
 ('FOS', 'ZFP36'),

]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_control_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_control_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_control_LOH_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
     ('UMOD', 'SLC12A1'),
 ('UMOD', 'AQP3'),
 ('UMOD', 'TPM1'),
 ('UMOD', 'CALD1'),
 ('UMOD', 'JAK1'),
    ('APOE', 'EZR'),
 ('SLC12A1', 'EGF'),
    ('FOS', 'JUNB'),
 ('FOS', 'ZFP36'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_dkd_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_dkd_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_dkd_loh_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


# Flatten the matrix
flat_diff = correlation_difference.stack().reset_index()
flat_diff.columns = ['Gene1', 'Gene2', 'CorrelationDifference']

# Remove self-correlations (diagonal)
flat_diff = flat_diff[flat_diff['Gene1'] != flat_diff['Gene2']]

# Drop duplicates for symmetrical pairs (e.g., A–B and B–A)
flat_diff['SortedPair'] = flat_diff.apply(lambda row: tuple(sorted([row['Gene1'], row['Gene2']])), axis=1)
flat_diff = flat_diff.drop_duplicates('SortedPair').drop(columns='SortedPair')

# Get top N positive and negative differences
top_n = 30
top_positive = flat_diff.sort_values(by='CorrelationDifference', ascending=False).head(top_n)
top_negative = flat_diff.sort_values(by='CorrelationDifference', ascending=True).head(top_n)

print("Top positive correlation differences:")
print(top_positive)

print("\nTop negative correlation differences:")
print(top_negative)


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap


# Define gene pairs
highlight_pairs = [
    ('COL3A1', 'MT2A'),
    ('COL1A1', 'MT2A'),
    ('TPT1', 'TNFRSF12A'),
    ('VEGFA', 'SOX4'),
    ('TNFRSF12A', 'RPL34'),
    ('COL3A1', 'SPP1'),
    ('TPT1', 'SOX4'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_difference_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_difference_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_differences_loh_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()
s









correlations_immune_dkd=pd.read_csv('/content/hotspot_local_correlations_immune_niche_dkd.csv')
correlations_immune_dkd.set_index('Unnamed: 0',inplace=True)
correlations_immune_dkd.index.name=None
#correlations_immune_dkd

correlations_immune_control=pd.read_csv('/content/hotspot_local_correlations_immune_niche_control.csv')
correlations_immune_control.set_index('Unnamed: 0',inplace=True)
correlations_immune_control.index.name=None
#correlations_immune_control

# Step 2: Find the common genes (row and column labels must match)
common_genes = correlations_immune_control.index.intersection(correlations_immune_dkd.index)

# Step 3: Subset both matrices to the common genes and ensure the same order
control_subset = correlations_immune_control.loc[common_genes, common_genes]
disease_subset = correlations_immune_dkd.loc[common_genes, common_genes]

correlation_difference = disease_subset - control_subset

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(control_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_control = pd.DataFrame(scaled_array, index=control_subset.index, columns=control_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(disease_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_dkd = pd.DataFrame(scaled_array, index=disease_subset.index, columns=disease_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(correlation_difference.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_difference = pd.DataFrame(scaled_array, index=correlation_difference.index, columns=correlation_difference.columns)

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.clustermap(
    scaled_correlation_control,
    cmap="coolwarm",
    vmin=-0.5,
    vmax=0.5,
    center=0,
    figsize=(20, 20),
    linewidths=0,
    cbar_kws={"label": "Δ correlation (disease - control)"},
    dendrogram_ratio=(0.001, 0.001),
    tree_kws={"linewidths": 0},
    cbar_pos=None
)

#plt.suptitle("Clustered Gene–Gene Correlation Differences (Disease vs. Control)", y=1.02)
#plt.savefig("clustered_correlation_difference_heatmap.png", dpi=300, bbox_inches='tight')
plt.close()
#plt.show()


# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_dkd.index[row_order]
col_labels = scaled_correlation_dkd.columns[col_order]

# Reorder DKD matrix
scaled_correlation_dkd_ordered = scaled_correlation_dkd.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_control.index[row_order]
col_labels = scaled_correlation_control.columns[col_order]

# Reorder DKD matrix
scaled_correlation_control_ordered = scaled_correlation_control.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_difference.index[row_order]
col_labels = scaled_correlation_difference.columns[col_order]

# Reorder DKD matrix
scaled_correlation_difference_ordered = scaled_correlation_difference.loc[row_labels, col_labels]

##correlated in control and dkd
import numpy as np
threshold = 0.99
high_corr_disease = scaled_correlation_dkd.abs() > threshold
high_corr_control = scaled_correlation_control.abs() > threshold
high_in_both = high_corr_disease & high_corr_control
high_pairs = np.argwhere(high_in_both.values)
gene_pairs = [(high_in_both.index[i], high_in_both.columns[j]) for i, j in high_pairs if i < j]  # avoid duplicates
len(gene_pairs)

gene_pairs

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('IGKC', 'SLC8A1'),
    ('GPX3', 'IRF4'),
    ('IGFBP7', 'P2RX5'),
    ('IGKC', 'SLC8A1'),
    ('IGKC', 'VEGFA'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_control_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_control_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_control_immune_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
       ('IGKC', 'SLC8A1'),
    ('GPX3', 'IRF4'),
    ('IGFBP7', 'P2RX5'),
    ('IGKC', 'SLC8A1'),
    ('IGKC', 'VEGFA'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_dkd_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_dkd_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_dkd_immune_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


# Flatten the matrix
flat_diff = correlation_difference.stack().reset_index()
flat_diff.columns = ['Gene1', 'Gene2', 'CorrelationDifference']

# Remove self-correlations (diagonal)
flat_diff = flat_diff[flat_diff['Gene1'] != flat_diff['Gene2']]

# Drop duplicates for symmetrical pairs (e.g., A–B and B–A)
flat_diff['SortedPair'] = flat_diff.apply(lambda row: tuple(sorted([row['Gene1'], row['Gene2']])), axis=1)
flat_diff = flat_diff.drop_duplicates('SortedPair').drop(columns='SortedPair')

# Get top N positive and negative differences
top_n = 30
top_positive = flat_diff.sort_values(by='CorrelationDifference', ascending=False).head(top_n)
top_negative = flat_diff.sort_values(by='CorrelationDifference', ascending=True).head(top_n)

print("Top positive correlation differences:")
print(top_positive)

print("\nTop negative correlation differences:")
print(top_negative)


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('IGHG1', 'IGHG2'),
    ('IGKC', 'IGHM'),
    ('IGHG2', 'IGKC'),
    ('IGHG1', 'IGKC'),
    ('IGHG2', 'IGHA1'),
    ('CD37', 'IGHD'),
    ('MS4A1', 'IGHD'),
    ('CD37', 'IGHM'),
    ('CD37', 'P2RX5'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_difference_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_difference_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_differences_immune_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()








correlations_iloh_dkd=pd.read_csv('/content/hotspot_local_correlations_iLOH_niche_dkd.csv')
correlations_iloh_dkd.set_index('Unnamed: 0',inplace=True)
correlations_iloh_dkd.index.name=None
#correlations_iloh_dkd

correlations_iloh_control=pd.read_csv('/content/hotspot_local_correlations_iLOH_niche_control.csv')
correlations_iloh_control.set_index('Unnamed: 0',inplace=True)
correlations_iloh_control.index.name=None
#correlations_iloh_control

# Step 2: Find the common genes (row and column labels must match)
common_genes = correlations_iloh_control.index.intersection(correlations_iloh_dkd.index)

# Step 3: Subset both matrices to the common genes and ensure the same order
control_subset = correlations_iloh_control.loc[common_genes, common_genes]
disease_subset = correlations_iloh_dkd.loc[common_genes, common_genes]

correlation_difference = disease_subset - control_subset

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(control_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_control = pd.DataFrame(scaled_array, index=control_subset.index, columns=control_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(disease_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_dkd = pd.DataFrame(scaled_array, index=disease_subset.index, columns=disease_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(correlation_difference.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_difference = pd.DataFrame(scaled_array, index=correlation_difference.index, columns=correlation_difference.columns)

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.clustermap(
    scaled_correlation_control,
    cmap="coolwarm",
    vmin=-0.5,
    vmax=0.5,
    center=0,
    figsize=(20, 20),
    linewidths=0,
    cbar_kws={"label": "Δ correlation (disease - control)"},
    dendrogram_ratio=(0.001, 0.001),
    tree_kws={"linewidths": 0},
    cbar_pos=None
)

#plt.suptitle("Clustered Gene–Gene Correlation Differences (Disease vs. Control)", y=1.02)
#plt.savefig("clustered_correlation_difference_heatmap.png", dpi=300, bbox_inches='tight')
plt.close()
#plt.show()


# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_dkd.index[row_order]
col_labels = scaled_correlation_dkd.columns[col_order]

# Reorder DKD matrix
scaled_correlation_dkd_ordered = scaled_correlation_dkd.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_control.index[row_order]
col_labels = scaled_correlation_control.columns[col_order]

# Reorder DKD matrix
scaled_correlation_control_ordered = scaled_correlation_control.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_difference.index[row_order]
col_labels = scaled_correlation_difference.columns[col_order]

# Reorder DKD matrix
scaled_correlation_difference_ordered = scaled_correlation_difference.loc[row_labels, col_labels]

##correlated in control and dkd
import numpy as np
threshold = 0.99
high_corr_disease = scaled_correlation_dkd.abs() > threshold
high_corr_control = scaled_correlation_control.abs() > threshold
high_in_both = high_corr_disease & high_corr_control
high_pairs = np.argwhere(high_in_both.values)
gene_pairs = [(high_in_both.index[i], high_in_both.columns[j]) for i, j in high_pairs if i < j]  # avoid duplicates
len(gene_pairs)

gene_pairs



import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('S100A2', 'SPP1'),
 ('S100A2', 'ITGB6'),
 ('UMOD', 'CLU'),
 ('UMOD', 'SLPI'),
 ('SLC12A1', 'S100A6'),
 ('S100A6', 'VCAM1'),
 ('SPP1', 'COL9A2'),
 ('IGKC', 'IGHA1'),
 ('IGKC', 'PPIA'),
 ('IGKC', 'JCHAIN'),
 ('TM4SF1', 'ADGRL2'),
 ('EGF', 'CXCL12'),
 ('ITGA3', 'MET'),
 ('IL32', 'VCAM1'),
 ('SERPINA1', 'CRYAB'),
 ('CYSTM1', 'MIF'),
 ('CDKN1A', 'EPHA2'),
 ('ENO1', 'PPIA'),
 ('HLA-DRB', 'CD74'),
 ('HLA-DRB', 'HLA-DRA'),
 ('CRIP1', 'EZR'),
 ('ATP5F1E', 'GLUL')
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_control_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_control_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_control_iLOH_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
           ('S100A2', 'SPP1'),
 ('S100A2', 'ITGB6'),
 ('UMOD', 'CLU'),
 ('UMOD', 'SLPI'),
 ('SLC12A1', 'S100A6'),
 ('S100A6', 'VCAM1'),
 ('SPP1', 'COL9A2'),
 ('IGKC', 'IGHA1'),
 ('IGKC', 'PPIA'),
 ('IGKC', 'JCHAIN'),
 ('TM4SF1', 'ADGRL2'),
 ('EGF', 'CXCL12'),
 ('ITGA3', 'MET'),
 ('IL32', 'VCAM1'),
 ('SERPINA1', 'CRYAB'),
 ('CYSTM1', 'MIF'),
 ('CDKN1A', 'EPHA2'),
 ('ENO1', 'PPIA'),
 ('HLA-DRB', 'CD74'),
 ('HLA-DRB', 'HLA-DRA'),
 ('CRIP1', 'EZR'),
 ('ATP5F1E', 'GLUL')
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_dkd_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_dkd_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_dkd_iLOH_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


# Flatten the matrix
flat_diff = correlation_difference.stack().reset_index()
flat_diff.columns = ['Gene1', 'Gene2', 'CorrelationDifference']

# Remove self-correlations (diagonal)
flat_diff = flat_diff[flat_diff['Gene1'] != flat_diff['Gene2']]

# Drop duplicates for symmetrical pairs (e.g., A–B and B–A)
flat_diff['SortedPair'] = flat_diff.apply(lambda row: tuple(sorted([row['Gene1'], row['Gene2']])), axis=1)
flat_diff = flat_diff.drop_duplicates('SortedPair').drop(columns='SortedPair')

# Get top N positive and negative differences
top_n = 30
top_positive = flat_diff.sort_values(by='CorrelationDifference', ascending=False).head(top_n)
top_negative = flat_diff.sort_values(by='CorrelationDifference', ascending=True).head(top_n)

print("Top positive correlation differences:")
print(top_positive)

print("\nTop negative correlation differences:")
print(top_negative)




import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
     ('MMP7', 'SLPI'),
    ('IGKC', 'IGHG1'),
    ('TACSTD2', 'SLPI'),
    ('IGKC', 'IGHG2'),
    ('IGKC', 'IGHA1'),
    ('CLU', 'SLPI'),
    ('CLU', 'MMP7'),
    ('TACSTD2', 'MMP7'),
    ('TACSTD2', 'CLU'),
    ('THBS1', 'CCL2'),
    ('TACSTD2', 'CRIP1'),
    ('CRIP1', 'EZR'),
    ('KRT7', 'TACSTD2'),
    ('KRT7', 'SLPI'),
    ('KRT7', 'CRIP1'),
    ('TPM1', 'SOD2'),
    ('TM4SF1', 'ATF3')
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_difference_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_difference_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_differences_iLOH_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()








correlations_fibro_dkd=pd.read_csv('/content/hotspot_local_correlations_fibroblast_niche_dkd.csv')
correlations_fibro_dkd.set_index('Unnamed: 0',inplace=True)
correlations_fibro_dkd.index.name=None
#correlations_fibro_dkd

correlations_fibro_control=pd.read_csv('/content/hotspot_local_correlations_fibroblast_niche_control.csv')
correlations_fibro_control.set_index('Unnamed: 0',inplace=True)
correlations_fibro_control.index.name=None
#correlations_fibro_control

# Step 2: Find the common genes (row and column labels must match)
common_genes = correlations_fibro_control.index.intersection(correlations_fibro_dkd.index)

# Step 3: Subset both matrices to the common genes and ensure the same order
control_subset = correlations_fibro_control.loc[common_genes, common_genes]
disease_subset = correlations_fibro_dkd.loc[common_genes, common_genes]

correlation_difference = disease_subset - control_subset

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(control_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_control = pd.DataFrame(scaled_array, index=control_subset.index, columns=control_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(disease_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_dkd = pd.DataFrame(scaled_array, index=disease_subset.index, columns=disease_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(correlation_difference.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_difference = pd.DataFrame(scaled_array, index=correlation_difference.index, columns=correlation_difference.columns)

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.clustermap(
    scaled_correlation_control,
    cmap="coolwarm",
    vmin=-0.5,
    vmax=0.5,
    center=0,
    figsize=(20, 20),
    linewidths=0,
    cbar_kws={"label": "Δ correlation (disease - control)"},
    dendrogram_ratio=(0.001, 0.001),
    tree_kws={"linewidths": 0},
    cbar_pos=None
)

#plt.suptitle("Clustered Gene–Gene Correlation Differences (Disease vs. Control)", y=1.02)
#plt.savefig("clustered_correlation_difference_heatmap.png", dpi=300, bbox_inches='tight')
plt.close()
#plt.show()


# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_dkd.index[row_order]
col_labels = scaled_correlation_dkd.columns[col_order]

# Reorder DKD matrix
scaled_correlation_dkd_ordered = scaled_correlation_dkd.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_control.index[row_order]
col_labels = scaled_correlation_control.columns[col_order]

# Reorder DKD matrix
scaled_correlation_control_ordered = scaled_correlation_control.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_difference.index[row_order]
col_labels = scaled_correlation_difference.columns[col_order]

# Reorder DKD matrix
scaled_correlation_difference_ordered = scaled_correlation_difference.loc[row_labels, col_labels]

##correlated in control and dkd
import numpy as np
threshold = 0.99
high_corr_disease = scaled_correlation_dkd.abs() > threshold
high_corr_control = scaled_correlation_control.abs() > threshold
high_in_both = high_corr_disease & high_corr_control
high_pairs = np.argwhere(high_in_both.values)
gene_pairs = [(high_in_both.index[i], high_in_both.columns[j]) for i, j in high_pairs if i < j]  # avoid duplicates
len(gene_pairs)

gene_pairs

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
     ('COL1A1', 'COL1A2'),
    ('COL1A1', 'COL5A1'),
    ('COL3A1', 'COL1A1'),
    ('TAGLN', 'ACTA2'),
    ('TAGLN', 'MGP'),
    ('TAGLN', 'DDIT3'),
    ('RGS5', 'NOTCH3'),
    ('RGS5', 'MYH11'),
    ('IGFBP7', 'TNFRSF14'),
    ('IGFBP7', 'CRYAB'),
    ('IGFBP7', 'IL6R'),
    ('IGFBP7', 'RORA'),
    ('DCN', 'PDGFRA'),
    ('DCN', 'EIF5A/L1'),
    ('MGP', 'PTGES2'),
    ('ENO1', 'PDGFRB'),
    ('TM4SF1', 'ADGRA2'),
    ('VIM', 'ITGA2')
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_control_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_control_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_control_fibro_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
           ('COL1A1', 'COL1A2'),
    ('COL1A1', 'COL5A1'),
    ('COL3A1', 'COL1A1'),
    ('TAGLN', 'ACTA2'),
    ('TAGLN', 'MGP'),
    ('TAGLN', 'DDIT3'),
    ('RGS5', 'NOTCH3'),
    ('RGS5', 'MYH11'),
    ('IGFBP7', 'TNFRSF14'),
    ('IGFBP7', 'CRYAB'),
    ('IGFBP7', 'IL6R'),
    ('IGFBP7', 'RORA'),
    ('DCN', 'PDGFRA'),
    ('DCN', 'EIF5A/L1'),
    ('MGP', 'PTGES2'),
    ('ENO1', 'PDGFRB'),
    ('TM4SF1', 'ADGRA2'),
    ('VIM', 'ITGA2')
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_dkd_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_dkd_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_dkd_fibro_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


# Flatten the matrix
flat_diff = correlation_difference.stack().reset_index()
flat_diff.columns = ['Gene1', 'Gene2', 'CorrelationDifference']

# Remove self-correlations (diagonal)
flat_diff = flat_diff[flat_diff['Gene1'] != flat_diff['Gene2']]

# Drop duplicates for symmetrical pairs (e.g., A–B and B–A)
flat_diff['SortedPair'] = flat_diff.apply(lambda row: tuple(sorted([row['Gene1'], row['Gene2']])), axis=1)
flat_diff = flat_diff.drop_duplicates('SortedPair').drop(columns='SortedPair')

# Get top N positive and negative differences
top_n = 30
top_positive = flat_diff.sort_values(by='CorrelationDifference', ascending=False).head(top_n)
top_negative = flat_diff.sort_values(by='CorrelationDifference', ascending=True).head(top_n)

print("Top positive correlation differences:")
print(top_positive)

print("\nTop negative correlation differences:")
print(top_negative)


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('JUNB', 'CCL2'),
    ('MT2A', 'JUNB'),
    ('IGHG1', 'IGHG2'),
    ('THBS1', 'JUNB'),
    ('THBS1', 'CCL2'),
    ('MT2A', 'THBS1'),
    ('FOS', 'MT2A'),
    ('VIM', 'MT2A'),
    ('COL1A1', 'MT2A'),
    ('JUNB', 'SERPINE1'),
    ('FOS', 'SERPINE1'),
    ('MT2A', 'CCL2'),
    ('COL3A1', 'MT2A'),
    ('MT2A', 'COL4A1'),
    ('CCL2', 'SERPINE1'),
    ('VIM', 'JUNB'),
    ('VIM', 'CCL2'),
    ('COL1A1', 'SERPINE1'),
    ('GDF15', 'SLPI'),
    ('FOS', 'SRGN'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_difference_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_difference_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_differences_fibro_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()








correlations_dct_dkd=pd.read_csv('/content/hotspot_local_correlations_dct_niche_dkd.csv')
correlations_dct_dkd.set_index('Unnamed: 0',inplace=True)
correlations_dct_dkd.index.name=None
#correlations_dct_dkd

correlations_dct_control=pd.read_csv('/content/hotspot_local_correlations_dct_niche_control.csv')
correlations_dct_control.set_index('Unnamed: 0',inplace=True)
correlations_dct_control.index.name=None
#correlations_dct_control

# Step 2: Find the common genes (row and column labels must match)
common_genes = correlations_dct_control.index.intersection(correlations_dct_dkd.index)

# Step 3: Subset both matrices to the common genes and ensure the same order
control_subset = correlations_dct_control.loc[common_genes, common_genes]
disease_subset = correlations_dct_dkd.loc[common_genes, common_genes]

correlation_difference = disease_subset - control_subset

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(control_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_control = pd.DataFrame(scaled_array, index=control_subset.index, columns=control_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(disease_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_dkd = pd.DataFrame(scaled_array, index=disease_subset.index, columns=disease_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(correlation_difference.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_difference = pd.DataFrame(scaled_array, index=correlation_difference.index, columns=correlation_difference.columns)

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.clustermap(
    scaled_correlation_control,
    cmap="coolwarm",
    vmin=-0.5,
    vmax=0.5,
    center=0,
    figsize=(20, 20),
    linewidths=0,
    cbar_kws={"label": "Δ correlation (disease - control)"},
    dendrogram_ratio=(0.001, 0.001),
    tree_kws={"linewidths": 0},
    cbar_pos=None
)

#plt.suptitle("Clustered Gene–Gene Correlation Differences (Disease vs. Control)", y=1.02)
#plt.savefig("clustered_correlation_difference_heatmap.png", dpi=300, bbox_inches='tight')
plt.close()
#plt.show()


# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_dkd.index[row_order]
col_labels = scaled_correlation_dkd.columns[col_order]

# Reorder DKD matrix
scaled_correlation_dkd_ordered = scaled_correlation_dkd.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_control.index[row_order]
col_labels = scaled_correlation_control.columns[col_order]

# Reorder DKD matrix
scaled_correlation_control_ordered = scaled_correlation_control.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_difference.index[row_order]
col_labels = scaled_correlation_difference.columns[col_order]

# Reorder DKD matrix
scaled_correlation_difference_ordered = scaled_correlation_difference.loc[row_labels, col_labels]

##correlated in control and dkd
import numpy as np
threshold = 0.99
high_corr_disease = scaled_correlation_dkd.abs() > threshold
high_corr_control = scaled_correlation_control.abs() > threshold
high_in_both = high_corr_disease & high_corr_control
high_pairs = np.argwhere(high_in_both.values)
gene_pairs = [(high_in_both.index[i], high_in_both.columns[j]) for i, j in high_pairs if i < j]  # avoid duplicates
len(gene_pairs)

gene_pairs

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
      ('CALB1', 'UMOD'),
    ('CALB1', 'SLC8A1'),
    ('CALB1', 'ALDH1A2'),
    ('SLC12A3', 'EGF'),
    ('SLC12A3', 'SPP1'),
    ('SLC12A3', 'SLC12A1'),
    ('SLC12A3', 'SAT1'),
    ('SLC12A3', 'LHX1'),
    ('SLC12A3', 'EGFR'),
    ('EGF', 'EGFR'),
    ('EGF', 'RARRES2'),
    ('SPP1', 'SCNN1G'),
    ('AQP2', 'AQP3'),
    ('AQP2', 'SCNN1G'),
    ('ACSM2B', 'FGFR3'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_control_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_control_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_control_dct_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
      ('CALB1', 'UMOD'),
    ('CALB1', 'SLC8A1'),
    ('CALB1', 'ALDH1A2'),
    ('SLC12A3', 'EGF'),
    ('SLC12A3', 'SPP1'),
    ('SLC12A3', 'SLC12A1'),
    ('SLC12A3', 'SAT1'),
    ('SLC12A3', 'LHX1'),
    ('SLC12A3', 'EGFR'),
    ('EGF', 'EGFR'),
    ('EGF', 'RARRES2'),
    ('SPP1', 'SCNN1G'),
    ('AQP2', 'AQP3'),
    ('AQP2', 'SCNN1G'),
    ('ACSM2B', 'FGFR3'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_dkd_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_dkd_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_dkd_dct_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


# Flatten the matrix
flat_diff = correlation_difference.stack().reset_index()
flat_diff.columns = ['Gene1', 'Gene2', 'CorrelationDifference']

# Remove self-correlations (diagonal)
flat_diff = flat_diff[flat_diff['Gene1'] != flat_diff['Gene2']]

# Drop duplicates for symmetrical pairs (e.g., A–B and B–A)
flat_diff['SortedPair'] = flat_diff.apply(lambda row: tuple(sorted([row['Gene1'], row['Gene2']])), axis=1)
flat_diff = flat_diff.drop_duplicates('SortedPair').drop(columns='SortedPair')

# Get top N positive and negative differences
top_n = 30
top_positive = flat_diff.sort_values(by='CorrelationDifference', ascending=False).head(top_n)
top_negative = flat_diff.sort_values(by='CorrelationDifference', ascending=True).head(top_n)

print("Top positive correlation differences:")
print(top_positive)

print("\nTop negative correlation differences:")
print(top_negative)


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('COL6A1', 'SLC13A3'),
    ('AZGP1', 'ADGRG1'),
    ('RPL34', 'IGHG1'),
    ('SLC13A3', 'ADGRG1'),
    ('SLC13A3', 'CD24'),
    ('TPT1', 'KRT17'),
    ('LRP2', 'CD24'),
    ('COL6A1', 'LRP2'),
    ('TPT1', 'POU5F1'),
    ('IGFBP7', 'SLC13A3'),
    ('TPT1', 'IGHG1'),
    ('SPP1', 'AQP2'),
    ('CD24', 'AZGP1'),
    ('IGHG1', 'RPL22'),
    ('IGFBP7', 'LRP2'),
    ('LRP2', 'ADGRG1'),
    ('ADGRG1', 'CRYAB'),
    ('S100A6', 'SLC13A3'),
    ('CUBN', 'ADGRG1'),
    ('RPL37', 'IGHG1'),
    ('RPL34', 'KRT16')

]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_difference_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_difference_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_differences_dct_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()












correlations_cnt_dkd=pd.read_csv('/content/hotspot_local_correlations_cnt_niche_dkd.csv')
correlations_cnt_dkd.set_index('Unnamed: 0',inplace=True)
correlations_cnt_dkd.index.name=None
#correlations_cnt_dkd

correlations_cnt_control=pd.read_csv('/content/hotspot_local_correlations_cnt_niche_control.csv')
correlations_cnt_control.set_index('Unnamed: 0',inplace=True)
correlations_cnt_control.index.name=None
#correlations_cnt_control

# Step 2: Find the common genes (row and column labels must match)
common_genes = correlations_cnt_control.index.intersection(correlations_cnt_dkd.index)

# Step 3: Subset both matrices to the common genes and ensure the same order
control_subset = correlations_cnt_control.loc[common_genes, common_genes]
disease_subset = correlations_cnt_dkd.loc[common_genes, common_genes]

correlation_difference = disease_subset - control_subset

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(control_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_control = pd.DataFrame(scaled_array, index=control_subset.index, columns=control_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(disease_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_dkd = pd.DataFrame(scaled_array, index=disease_subset.index, columns=disease_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(correlation_difference.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_difference = pd.DataFrame(scaled_array, index=correlation_difference.index, columns=correlation_difference.columns)

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.clustermap(
    scaled_correlation_control,
    cmap="coolwarm",
    vmin=-0.5,
    vmax=0.5,
    center=0,
    figsize=(20, 20),
    linewidths=0,
    cbar_kws={"label": "Δ correlation (disease - control)"},
    dendrogram_ratio=(0.001, 0.001),
    tree_kws={"linewidths": 0},
    cbar_pos=None
)

#plt.suptitle("Clustered Gene–Gene Correlation Differences (Disease vs. Control)", y=1.02)
#plt.savefig("clustered_correlation_difference_heatmap.png", dpi=300, bbox_inches='tight')
plt.close()
#plt.show()


# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_dkd.index[row_order]
col_labels = scaled_correlation_dkd.columns[col_order]

# Reorder DKD matrix
scaled_correlation_dkd_ordered = scaled_correlation_dkd.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_control.index[row_order]
col_labels = scaled_correlation_control.columns[col_order]

# Reorder DKD matrix
scaled_correlation_control_ordered = scaled_correlation_control.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_difference.index[row_order]
col_labels = scaled_correlation_difference.columns[col_order]

# Reorder DKD matrix
scaled_correlation_difference_ordered = scaled_correlation_difference.loc[row_labels, col_labels]

##correlated in control and dkd
import numpy as np
threshold = 0.99
high_corr_disease = scaled_correlation_dkd.abs() > threshold
high_corr_control = scaled_correlation_control.abs() > threshold
high_in_both = high_corr_disease & high_corr_control
high_pairs = np.argwhere(high_in_both.values)
gene_pairs = [(high_in_both.index[i], high_in_both.columns[j]) for i, j in high_pairs if i < j]  # avoid duplicates
len(gene_pairs)

gene_pairs

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('AQP2', 'AQP3'),
    ('CALB1', 'SLC8A1'),
    ('CALB1', 'UMOD'),
    ('AQP3', 'SLC12A3'),
    ('SLC12A3', 'EGF'),
    ('SPP1', 'AQP2'),
    ('SPP1', 'AQP3'),
    ('SPP1', 'SCNN1G'),
    ('AQP2', 'SCNN1G'),
    ('AQP2', 'KRT19'),
    ('AQP2', 'ST6GAL1')
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_control_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_control_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_control_cnt_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
       ('AQP2', 'AQP3'),
    ('CALB1', 'SLC8A1'),
    ('CALB1', 'UMOD'),
    ('AQP3', 'SLC12A3'),
    ('SLC12A3', 'EGF'),
    ('SPP1', 'AQP2'),
    ('SPP1', 'AQP3'),
    ('SPP1', 'SCNN1G'),
    ('AQP2', 'SCNN1G'),
    ('AQP2', 'KRT19'),
    ('AQP2', 'ST6GAL1'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_dkd_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_dkd_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_dkd_cnt_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


# Flatten the matrix
flat_diff = correlation_difference.stack().reset_index()
flat_diff.columns = ['Gene1', 'Gene2', 'CorrelationDifference']

# Remove self-correlations (diagonal)
flat_diff = flat_diff[flat_diff['Gene1'] != flat_diff['Gene2']]

# Drop duplicates for symmetrical pairs (e.g., A–B and B–A)
flat_diff['SortedPair'] = flat_diff.apply(lambda row: tuple(sorted([row['Gene1'], row['Gene2']])), axis=1)
flat_diff = flat_diff.drop_duplicates('SortedPair').drop(columns='SortedPair')

# Get top N positive and negative differences
top_n = 30
top_positive = flat_diff.sort_values(by='CorrelationDifference', ascending=False).head(top_n)
top_negative = flat_diff.sort_values(by='CorrelationDifference', ascending=True).head(top_n)

print("Top positive correlation differences:")
print(top_positive)

print("\nTop negative correlation differences:")
print(top_negative)


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
   ('TNFRSF12A', 'VEGFA'),
    ('COL6A1', 'TNFRSF12A'),
    ('SLC8A1', 'TNFRSF12A'),
    ('CALB1', 'APOE'),
    ('CALB1', 'ACSM2B'),
    ('TNFRSF12A', 'GAS6'),
    ('CALB1', 'GPX3'),
    ('SOD2', 'SLPI'),
    ('TNFRSF12A', 'ITGA6'),
    ('CALB1', 'SLC13A3'),
    ('B2M', 'TNFRSF12A'),
    ('B2M', 'MT2A'),
    ('MIF', 'TNFRSF12A'),
    ('CALB1', 'MAF'),
    ('TNFRSF12A', 'COL18A1'),
    ('TNFRSF12A', 'SLC26A7'),
    ('CALB1', 'AZGP1'),
    ('TNFRSF12A', 'SPINK1'),
    ('TNFRSF12A', 'ADGRF5'),
    ('CALB1', 'CXCL14'),
    ('TNFRSF12A', 'ADGRG1'),
    ('MMP7', 'SLPI'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_difference_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_difference_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_differences_cnt_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()






import pandas as pd

correlations_vascular_dkd=pd.read_csv('/content/hotspot_local_correlations_vascular_niche_dkd.csv')
correlations_vascular_dkd.set_index('Unnamed: 0',inplace=True)
correlations_vascular_dkd.index.name=None
#correlations_vascular_dkd

correlations_vascular_control=pd.read_csv('/content/hotspot_local_correlations_vascular_niche_control.csv')
correlations_vascular_control.set_index('Unnamed: 0',inplace=True)
correlations_vascular_control.index.name=None
#correlations_vascular_control

# Step 2: Find the common genes (row and column labels must match)
common_genes = correlations_vascular_control.index.intersection(correlations_vascular_dkd.index)

# Step 3: Subset both matrices to the common genes and ensure the same order
control_subset = correlations_vascular_control.loc[common_genes, common_genes]
disease_subset = correlations_vascular_dkd.loc[common_genes, common_genes]

#control_subset

#disease_subset

correlation_difference = disease_subset - control_subset
correlation_difference

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(control_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_control = pd.DataFrame(scaled_array, index=control_subset.index, columns=control_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(disease_subset.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_dkd = pd.DataFrame(scaled_array, index=disease_subset.index, columns=disease_subset.columns)

from sklearn.preprocessing import MinMaxScaler
import pandas as pd

# Scale the full matrix to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaled_array = scaler.fit_transform(correlation_difference.values)

# Convert back to DataFrame with original index/columns
scaled_correlation_difference = pd.DataFrame(scaled_array, index=correlation_difference.index, columns=correlation_difference.columns)

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.clustermap(
    scaled_correlation_control,
    cmap="coolwarm",
    vmin=-0.5,
    vmax=0.5,
    center=0,
    figsize=(20, 20),
    linewidths=0,
    cbar_kws={"label": "Δ correlation (disease - control)"},
    dendrogram_ratio=(0.001, 0.001),
    tree_kws={"linewidths": 0},
    cbar_pos=None
)

#plt.suptitle("Clustered Gene–Gene Correlation Differences (Disease vs. Control)", y=1.02)
#plt.savefig("clustered_correlation_difference_heatmap.png", dpi=300, bbox_inches='tight')
plt.close()
#plt.show()


# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_dkd.index[row_order]
col_labels = scaled_correlation_dkd.columns[col_order]

# Reorder DKD matrix
scaled_correlation_dkd_ordered = scaled_correlation_dkd.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_control.index[row_order]
col_labels = scaled_correlation_control.columns[col_order]

# Reorder DKD matrix
scaled_correlation_control_ordered = scaled_correlation_control.loc[row_labels, col_labels]

# Extract the row and column order from the dendrogram
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Use them to reorder your DKD matrix
row_labels = scaled_correlation_difference.index[row_order]
col_labels = scaled_correlation_difference.columns[col_order]

# Reorder DKD matrix
scaled_correlation_difference_ordered = scaled_correlation_difference.loc[row_labels, col_labels]

##correlated in control and dkd
import numpy as np
threshold = 0.99
high_corr_disease = scaled_correlation_dkd.abs() > threshold
high_corr_control = scaled_correlation_control.abs() > threshold
high_in_both = high_corr_disease & high_corr_control
high_pairs = np.argwhere(high_in_both.values)
gene_pairs = [(high_in_both.index[i], high_in_both.columns[j]) for i, j in high_pairs if i < j]  # avoid duplicates
len(gene_pairs)

gene_pairs



import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
('MYH11', 'SPARCL1'),
    ('ACTA2', 'TAGLN'),
    ('ACTA2', 'NOTCH3'),
    ('TPM2', 'NOTCH3'),
    ('TAGLN', 'CAV1'),
    ('TAGLN', 'ADGRG1'),
    ('TAGLN', 'VHL'),
    ('VIM', 'EZR'),
    ('CXCL12', 'LHX1'),
    ('PECAM1', 'ICAM2'),
    ('PECAM1', 'TEK'),
    ('PECAM1', 'ADGRL4'),
    ('GPX3', 'PDGFC'),
    ('COL4A1', 'COL4A2'),
    ('COL18A1', 'ACKR1'),
    ('RBM47', 'VEGFA'),
    ('TM4SF1', 'PECAM1'),
    ('JAG1', 'RAMP3'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_control_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_control_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_control_vascular_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
('MYH11', 'SPARCL1'),
    ('ACTA2', 'TAGLN'),
    ('ACTA2', 'NOTCH3'),
    ('TPM2', 'NOTCH3'),
    ('TAGLN', 'CAV1'),
    ('TAGLN', 'ADGRG1'),
    ('TAGLN', 'VHL'),
    ('VIM', 'EZR'),
    ('CXCL12', 'LHX1'),
    ('PECAM1', 'ICAM2'),
    ('PECAM1', 'TEK'),
    ('PECAM1', 'ADGRL4'),
    ('GPX3', 'PDGFC'),
    ('COL4A1', 'COL4A2'),
    ('COL18A1', 'ACKR1'),
    ('RBM47', 'VEGFA'),
    ('TM4SF1', 'PECAM1'),
    ('JAG1', 'RAMP3'),
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_dkd_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_dkd_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_dkd_vascular_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()


# Flatten the matrix
flat_diff = correlation_difference.stack().reset_index()
flat_diff.columns = ['Gene1', 'Gene2', 'CorrelationDifference']

# Remove self-correlations (diagonal)
flat_diff = flat_diff[flat_diff['Gene1'] != flat_diff['Gene2']]

# Drop duplicates for symmetrical pairs (e.g., A–B and B–A)
flat_diff['SortedPair'] = flat_diff.apply(lambda row: tuple(sorted([row['Gene1'], row['Gene2']])), axis=1)
flat_diff = flat_diff.drop_duplicates('SortedPair').drop(columns='SortedPair')

# Get top N positive and negative differences
top_n = 30
top_positive = flat_diff.sort_values(by='CorrelationDifference', ascending=False).head(top_n)
top_negative = flat_diff.sort_values(by='CorrelationDifference', ascending=True).head(top_n)

print("Top positive correlation differences:")
print(top_positive)

print("\nTop negative correlation differences:")
print(top_negative)


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_hex
from matplotlib.cm import get_cmap

# Define gene pairs
highlight_pairs = [
    ('THBS1', 'JUNB'),
    ('FOS', 'MT2A'),
    ('COL1A1', 'MT2A'),
    ('MT2A', 'JUNB'),
    ('JUNB', 'SERPINE1'),
    ('FOS', 'THBS1'),
    ('JUNB', 'CCL2'),
    ('COL1A1', 'SERPINE1'),
    ('THBS1', 'KLF2'),
    ('DUSP1', 'THBS1'),
    ('MT2A', 'COL4A1'),
    ('COL4A1', 'SERPINE1'),
    ('CXCL12', 'MHC I'),
    ('JUNB', 'HLA-DRB'),
    ('JUNB', 'CD74'),
    ('COL1A1', 'THBS1'),
    ('RGS5', 'FOS')
]

# Assign one color per pair; map each gene only once based on first appearance
tab20 = get_cmap('tab20', len(highlight_pairs))
gene_to_color = {}
for i, (g1, g2) in enumerate(highlight_pairs):
    color = to_hex(tab20(i))
    for gene in (g1, g2):
        if gene not in gene_to_color:
            gene_to_color[gene] = color

# Extract gene ordering from DataFrame
all_genes = scaled_correlation_difference_ordered.columns.tolist()

# Get indices of x and y genes in the DataFrame
x_genes = [pair[0] for pair in highlight_pairs]
y_genes = [pair[1] for pair in highlight_pairs]
x_indices = [i for i, g in enumerate(all_genes) if g in x_genes]
y_indices = [i for i, g in enumerate(all_genes) if g in y_genes]

# Plot heatmap with no default tick labels
plt.figure(figsize=(12, 12))
ax = sns.heatmap(
    scaled_correlation_difference_ordered,
    cmap="coolwarm",
    vmin=-0.8, vmax=0.8, center=0,
    xticklabels=False,
    yticklabels=False,
    cbar=False,
    linewidths=0
)

# Set custom ticks and labels
ax.set_xticks(x_indices)
ax.set_xticklabels([all_genes[i] for i in x_indices], rotation=90, fontsize=12)
ax.set_yticks(y_indices)
ax.set_yticklabels([all_genes[i] for i in y_indices], rotation=0, fontsize=12)

# Define 5 stagger positions for offsetting labels
x_offsets = [0, -0.14, -0.28, -0.42, -0.56]
y_offsets = [0, -0.14, -0.28, -0.42, -0.56]

# Color and stagger x-axis tick labels
for i, label in enumerate(ax.get_xticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_y(x_offsets[i % 5])

# Color and stagger y-axis tick labels
for i, label in enumerate(ax.get_yticklabels()):
    gene = label.get_text()
    if gene in gene_to_color:
        label.set_color(gene_to_color[gene])
        label.set_x(-0.6 - y_offsets[i % 5])  # Further left

plt.title("")
plt.tight_layout()
plt.savefig("correlation_differences_vascular_highlighted_axes.png", dpi=900, bbox_inches='tight')
plt.show()

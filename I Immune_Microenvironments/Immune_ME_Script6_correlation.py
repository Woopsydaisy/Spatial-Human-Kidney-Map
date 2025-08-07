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

# Load and preprocess sample metadata
sample_metadata = pd.read_csv('/content/diagnosis.csv', sep=';')
sample_metadata = sample_metadata.set_index('Sample ID')
sample_metadata = sample_metadata.drop(['Disease', 'Race', 'Sex', 'DM', 'HTN'], axis=1)

# Define conditions to exclude
excluded_conditions = [
    'DM', 'HTN', 'CKD', 'C3GN', 'AA amyloid', 'IgA', 'FSGS', 
    'MN', 'DKD+FSGS', 'TMA', 'DM/HTN'
]

# Filter out excluded conditions
sample_metadata_updated = sample_metadata[~sample_metadata['Condition'].isin(excluded_conditions)]

# Drop unnecessary columns
sample_metadata_updated = sample_metadata_updated.drop(['Condition', 'Age'], axis=1)

sample_metadata_updated


df=pd.read_csv('/content/statistical_analysis_interaction_scores_04_21_2025_152700.txt',sep='\t', low_memory=False)

df=df.set_index('interacting_pair')

df_subset=df.filter(regex='^Glomerular Immune ME')

niches_to_filter=['^Glomerular Immune ME']
celltypes_to_send=['CD8+']
celltypes_to_receive=['MC1']
for niche in niches_to_filter:
  df_filter=df.filter(regex=niche)
  for receiver in celltypes_to_receive:
    df_filter_receiver=df_filter.loc[:, df_filter.columns.str.endswith(f'_{receiver}')]
    for sender in celltypes_to_send:
      print(f'{niche}:{sender} to {receiver}:')
      if sender == receiver:
        # Require two instances of the cell type in the column name
        df_filter_receiver_sender = df_filter_receiver.loc[:, df_filter_receiver.columns.str.count(f'_{sender}') == 2]
      else:
        df_filter_receiver_sender = df_filter_receiver.loc[:, df_filter_receiver.columns.str.contains(f'_{sender}')]
      if df_filter_receiver_sender.empty:
        print('empty')
        continue
      extracted = df_filter_receiver_sender.columns.str.extract(r'ME_(.*?)_(Xenium|CosMx)')
      # Combine them to form 'Sample_Platform' format
      df_filter_receiver_sender.columns = extracted[0] + '_' + extracted[1]
      df_pivoted = df_filter_receiver_sender.T
      df_pivoted.index.name = 'Sample ID'  # Just to be explicit, optional
      # Merge on index (Sample ID)
      merged_df = df_pivoted.merge(sample_metadata_updated, left_index=True, right_index=True)

      filtered_df_merged = merged_df.loc[:, (merged_df != 0).sum(axis=0) > 8]
      filtered_df_merged = filtered_df_merged.loc[:, filtered_df_merged.mean(axis=0) > 0.2]
      gfr_correlations = filtered_df_merged.corrwith(filtered_df_merged['GFR'])

      # Optionally: sort correlations
      gfr_correlations_sorted = gfr_correlations.sort_values(ascending=False)
      #gfr_correlations_sorted = gfr_correlations_sorted[gfr_correlations_sorted.index.str.contains('FGF1')]
      if len(gfr_correlations_sorted) == 0:
        continue

      # Get top 20 features most correlated with GFR (excluding GFR itself)
      top_20_features = gfr_correlations_sorted.drop(labels='GFR').head(20).index #
      # Set up the plot grid
      n_cols = 5
      n_rows = (len(top_20_features) + n_cols - 1) // n_cols
      fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 2.8), constrained_layout=True)

      # Flatten axes for easier indexing
      axes = axes.flatten()

      # Plot each feature
      for i, feature in enumerate(top_20_features):
          ax = axes[i]
          sns.regplot(data=merged_df, x=feature, y='GFR', ax=ax,
            scatter_kws={'s': 10, 'color': 'black'},
            line_kws={"color": "black"}
          )
          ax.set_xlabel('Interaction score', color='black', fontsize=12)
          ax.set_ylabel('GFR', color='black', fontsize=12)
          ax.set_ylim(0, 120)
          ax.set_title(f"{feature}\nr = {gfr_correlations_sorted[feature]:.2f}", color='black', fontsize=10)

      # Style the plot borders and ticks
      ax.tick_params(colors='black')
      for spine in ax.spines.values():
          spine.set_color('black')

      # Turn off any unused subplots
      for j in range(i + 1, len(axes)):
          axes[j].axis('off')

      # Show or save the plot
      plt.suptitle(f"Top 20 in {niche} for {sender} signaling to {receiver} corr with GFR", fontsize=16, y=1.02)
      #plt.savefig("top20_gfr_correlations.png", dpi=900, bbox_inches='tight')
      plt.show()
      # Get top 20 features most correlated with GFR (excluding GFR itself)
      top_20_features = gfr_correlations_sorted.drop(labels='GFR').tail(60).index #


      # Set up the plot grid
      n_cols = 5
      n_rows = (len(top_20_features) + n_cols - 1) // n_cols
      fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 2.8), constrained_layout=True)

      # Flatten axes for easier indexing
      axes = axes.flatten()

      for i, feature in enumerate(top_20_features):
          ax = axes[i]
          sns.regplot(data=merged_df, x=feature, y='GFR', ax=ax,
            scatter_kws={'s': 10, 'color': 'black'},
            line_kws={"color": "black"}
          )
          ax.set_xlabel('Interaction score', color='black', fontsize=12)
          ax.set_ylabel('GFR', color='black', fontsize=12)
          ax.set_ylim(0, 120)
          ax.set_title(f"{feature}\nr = {gfr_correlations_sorted[feature]:.2f}", color='black', fontsize=10)

      # Style the plot borders and ticks
      ax.tick_params(colors='black')
      for spine in ax.spines.values():
          spine.set_color('black')

      # Turn off any unused subplots
      for j in range(i + 1, len(axes)):
          axes[j].axis('off')

      # Show or save the plot
      plt.suptitle(f"Top 20 in {niche} for {sender} signaling to {receiver} corr with GFR", fontsize=16, y=1.02)
      #plt.savefig("top20_gfr_negative_correlations_mc1_to_macro_glom_immune_me.png", dpi=900, bbox_inches='tight')
      plt.show()

import csv
import gzip
import os
import scipy.io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

df=pd.read_csv('/content/statistical_analysis_interaction_scores_04_11_2025_160422.txt',sep='\t', low_memory=False)
df=df.set_index('interacting_pair')

##this is the expression of HAVCR1 and COL1A1 for the correlation with interaction scores:
expression_data=pd.read_csv('/content/geneexpression_microenvironment_for_correlation.csv')
expression_data=expression_data.set_index('microenvironment_cell_types')
expression_data_iPT=expression_data.drop(['COL1A1'], axis=1)
expression_data_iPT=expression_data_iPT.T
expression_data_iPT
expression_data_fibro=expression_data.drop(['HAVCR1'], axis=1)
expression_data_fibro=expression_data_fibro.T
expression_data_iPT=expression_data_iPT.filter(regex='^Profibrotic ME')
expression_data_iPT=expression_data_iPT.filter(regex='(_iPT)$')
extracted = expression_data_iPT.columns.str.extract(r'ME_(.*?)_(Xenium|CosMx)')
expression_data_iPT.columns = extracted[0] + '_' + extracted[1]
expression_data_iPT=expression_data_iPT.T
expression_data_iPT
expression_data_fibro=expression_data_fibro.filter(regex='^Profibrotic ME')
expression_data_fibro=expression_data_fibro.filter(regex='(_Fibroblast)$')
extracted = expression_data_fibro.columns.str.extract(r'ME_(.*?)_(Xenium|CosMx)')
expression_data_fibro.columns = extracted[0] + '_' + extracted[1]
expression_data_fibro=expression_data_fibro.T
expression_data_fibro

niches_to_filter=['^Profibrotic ME']
celltypes_to_send=['iPT']
celltypes_to_receive=['Fibroblast']
for niche in niches_to_filter:
  df_filter=df.filter(regex=niche)
  for receiver in celltypes_to_receive:
    df_filter_receiver=df_filter.filter(regex=f'(_{receiver})$')
    for sender in celltypes_to_send:
      print(f'{niche}:{sender} to {receiver}:')
      if sender == receiver:
        # Require two instances of the cell type in the column name
        df_filter_receiver_sender = df_filter_receiver.loc[:, df_filter_receiver.columns.str.count(f'_{sender}') == 2]
      else:
        df_filter_receiver_sender = df_filter_receiver.loc[:, df_filter_receiver.columns.str.contains(f'_{sender}')]
      if df_filter_receiver_sender.empty:
        continue
      extracted = df_filter_receiver_sender.columns.str.extract(r'ME_(.*?)_(Xenium|CosMx)')
      # Combine them to form 'Sample_Platform' format
      df_filter_receiver_sender.columns = extracted[0] + '_' + extracted[1]
      df_pivoted = df_filter_receiver_sender.T
      df_pivoted.index.name = 'Sample ID'  # Just to be explicit, optional
      # Merge on index (Sample ID)
      merged_df = df_pivoted.merge(expression_data_fibro, left_index=True, right_index=True)

      filtered_df_merged = merged_df.loc[:, (merged_df != 0).sum(axis=0) > 8]
      gfr_correlations = filtered_df_merged.corrwith(filtered_df_merged['COL1A1'])

      # Optionally: sort correlations
      gfr_correlations_sorted = gfr_correlations.sort_values(ascending=False)

      # Get top 20 features most correlated with GFR (excluding GFR itself)
      top_20_features = gfr_correlations_sorted.drop(labels=expression_data_fibro.columns).head(20).index #'COL1A1'
      # Set up the plot grid
      n_cols = 5
      n_rows = (len(top_20_features) + n_cols - 1) // n_cols
      fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 2, n_rows * 2.9), constrained_layout=True)

      # Flatten axes for easier indexing
      axes = axes.flatten()

      # Plot each feature
      for i, feature in enumerate(top_20_features):
          ax = axes[i]
          sns.regplot(data=merged_df, x='COL1A1', y=feature, ax=ax, scatter_kws={'s': 15, 'color': 'black'}, line_kws={"color": "black"})
          ax.set_title(f"{feature}\nr={gfr_correlations_sorted[feature]:.2f}")
          ax.set_xlabel('COL1A1', fontsize=16)
          ax.set_ylabel('Interaction score', fontsize=10)

      # Turn off any unused subplots
      for j in range(i + 1, len(axes)):
          axes[j].axis('off')

      # Show or save the plot
      plt.suptitle(f"Top 20 in {niche} for {sender} signaling to {receiver} corr with COL1A1",y=1.02)
      plt.savefig("top20_col1a1_correlations.png", dpi=900, bbox_inches='tight')
      plt.show()
      # Get top 20 features most correlated with GFR (excluding GFR itself)
      top_20_features = gfr_correlations_sorted.drop(labels='COL1A1').tail(20).index

      # Set up the plot grid
      n_cols = 5
      n_rows = (len(top_20_features) + n_cols - 1) // n_cols
      fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 5, n_rows * 5), constrained_layout=True)

      # Flatten axes for easier indexing
      axes = axes.flatten()

      # Plot each feature
      for i, feature in enumerate(top_20_features):
          ax = axes[i]
          sns.regplot(data=merged_df, x='COL1A1', y=feature, ax=ax, scatter_kws={'s': 15}, line_kws={"color": "black"})
          ax.set_title(f"{feature}\nr={gfr_correlations_sorted[feature]:.2f}")
          ax.set_xlabel('COL1A1')
          ax.set_ylabel('Interaction score')

      # Turn off any unused subplots
      for j in range(i + 1, len(axes)):
          axes[j].axis('off')

      # Show or save the plot
      plt.suptitle(f"Top 20 in {niche} for {sender} signaling to {receiver} corr with COL1A1", y=1.02)
      plt.show()


niches_to_filter=['^Profibrotic ME']
celltypes_to_send=['iPT']
celltypes_to_receive=['Fibroblast']
for niche in niches_to_filter:
  df_filter=df.filter(regex=niche)
  for receiver in celltypes_to_receive:
    df_filter_receiver=df_filter.filter(regex=f'(_{receiver})$')
    for sender in celltypes_to_send:
      print(f'{niche}:{sender} to {receiver}:')
      if sender == receiver:
        # Require two instances of the cell type in the column name
        df_filter_receiver_sender = df_filter_receiver.loc[:, df_filter_receiver.columns.str.count(f'_{sender}') == 2]
      else:
        df_filter_receiver_sender = df_filter_receiver.loc[:, df_filter_receiver.columns.str.contains(f'_{sender}')]
      if df_filter_receiver_sender.empty:
        continue
      extracted = df_filter_receiver_sender.columns.str.extract(r'ME_(.*?)_(Xenium|CosMx)')
      # Combine them to form 'Sample_Platform' format
      df_filter_receiver_sender.columns = extracted[0] + '_' + extracted[1]
      df_pivoted = df_filter_receiver_sender.T
      df_pivoted.index.name = 'Sample ID'  # Just to be explicit, optional
      # Merge on index (Sample ID)
      merged_df = df_pivoted.merge(expression_data_iPT, left_index=True, right_index=True)

      filtered_df_merged = merged_df.loc[:, (merged_df != 0).sum(axis=0) > 8]
      gfr_correlations = filtered_df_merged.corrwith(filtered_df_merged['HAVCR1'])

      # Optionally: sort correlations
      gfr_correlations_sorted = gfr_correlations.sort_values(ascending=False)

      # Get top 20 features most correlated with GFR (excluding GFR itself)
      top_20_features = gfr_correlations_sorted.drop(labels=expression_data_iPT.columns).head(20).index
      # Set up the plot grid
      n_cols = 5
      n_rows = (len(top_20_features) + n_cols - 1) // n_cols
      fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 2, n_rows * 2.9), constrained_layout=True)

      # Flatten axes for easier indexing
      axes = axes.flatten()

      # Plot each feature
      for i, feature in enumerate(top_20_features):
          ax = axes[i]
          sns.regplot(data=merged_df, x='HAVCR1', y=feature, ax=ax, scatter_kws={'s': 15, 'color': 'black'}, line_kws={"color": "black"})
          ax.set_title(f"{feature}\nr={gfr_correlations_sorted[feature]:.2f}")
          ax.set_xlabel('HAVCR1', fontsize=16)
          ax.set_ylabel('Interaction score', fontsize=16)

      # Turn off any unused subplots
      for j in range(i + 1, len(axes)):
          axes[j].axis('off')

      # Show or save the plot
      plt.suptitle(f"Top 20 in {niche} for {sender} signaling to {receiver} corr with HAVCR1", fontsize=12, y=1.02)
      plt.savefig("top20_correlations_with_havcr1.png", dpi=900, bbox_inches='tight')
      plt.show()
      # Get top 20 features most correlated with GFR (excluding GFR itself)
      top_20_features = gfr_correlations_sorted.drop(labels='HAVCR1').tail(20).index

      # Set up the plot grid
      n_cols = 5
      n_rows = (len(top_20_features) + n_cols - 1) // n_cols
      fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 5), constrained_layout=True)

      # Flatten axes for easier indexing
      axes = axes.flatten()

      # Plot each feature
      for i, feature in enumerate(top_20_features):
          ax = axes[i]
          sns.regplot(data=merged_df, x='HAVCR1', y=feature, ax=ax, scatter_kws={'s': 10, 'color': 'black'}, line_kws={"color": "black"})
          ax.set_title(f"{feature}\nr={gfr_correlations_sorted[feature]:.2f}")
          ax.set_xlabel('HAVCR1')
          ax.set_ylabel('Interaction score')

      # Turn off any unused subplots
      for j in range(i + 1, len(axes)):
          axes[j].axis('off')

      # Show or save the plot
      plt.suptitle(f"Top 20 in {niche} for {sender} signaling to {receiver} corr with HAVCR1", fontsize=16, y=1.02)
      plt.show()

import csv
import gzip
import os
import scipy.io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

##load cellphoneDB output
df=pd.read_csv('/content/statistical_analysis_interaction_scores_04_07_2025_181923.txt', sep='\t', low_memory=False)
df=df.set_index('interacting_pair')

# Load and clean sample metadata
sample_metadata = pd.read_csv('/content/diagnosis.csv', sep=';')
sample_metadata = sample_metadata.set_index('Sample ID')
# Drop unnecessary columns
columns_to_drop = ['Disease', 'Race', 'Sex', 'DM', 'HTN']
sample_metadata = sample_metadata.drop(columns=columns_to_drop)
# Exclude unwanted conditions
excluded_conditions = [
    'DM', 'HTN', 'CKD', 'C3GN', 'AA amyloid', 'IgA', 'FSGS',
    'MN', 'DKD+FSGS', 'TMA', 'DM/HTN'
]
sample_metadata = sample_metadata[~sample_metadata['Condition'].isin(excluded_conditions)]

# Remove rows with missing or zero GFR
sample_metadata = sample_metadata.fillna(0)
sample_metadata = sample_metadata[sample_metadata['GFR'] != 0]

# Optional: inspect condition distribution
print(sample_metadata['Condition'].value_counts())

# Drop final columns
sample_metadata = sample_metadata.drop(columns=['Condition', 'Age'])

# Display cleaned DataFrame
sample_metadata


##Glomerular Niche
niches_to_filter=['^Glomerular niche']
celltypes_to_send=['EC_glom']
celltypes_to_receive=['Podo']         
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
      extracted = df_filter_receiver_sender.columns.str.extract(r'niche_(.*?)_(Xenium|CosMx)')
      # Combine them to form 'Sample_Platform' format
      df_filter_receiver_sender.columns = extracted[0] + '_' + extracted[1]
      df_pivoted = df_filter_receiver_sender.T
      df_pivoted.index.name = 'Sample ID'  # Just to be explicit, optional
      # Merge on index (Sample ID)
      merged_df = df_pivoted.merge(sample_metadata, left_index=True, right_index=True)

      filtered_df_merged = merged_df.loc[:, (merged_df != 0).sum(axis=0) > 8]
      gfr_correlations = filtered_df_merged.corrwith(filtered_df_merged['GFR'])

      # Optionally: sort correlations
      gfr_correlations_sorted = gfr_correlations.sort_values(ascending=False)

      # Get top 20 features most correlated with GFR (excluding GFR itself)
      top_20_features = gfr_correlations_sorted.drop(labels='GFR').head(5).index
      # Set up the plot grid
      n_cols = 5
      n_rows = (len(top_20_features) + n_cols - 1) // n_cols
      fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 3), constrained_layout=True)

      # Flatten axes for easier indexing
      axes = axes.flatten()

      # Plot each feature
      for i, feature in enumerate(top_20_features):
          ax = axes[i]
          sns.regplot(data=merged_df, x='GFR', y=feature, ax=ax, scatter_kws={'s': 10}, line_kws={"color": "red"})
          ax.set_title(f"{feature}\nr={gfr_correlations_sorted[feature]:.2f}")
          ax.set_xlabel('GFR')
          ax.set_ylabel('Interaction score')

      # Turn off any unused subplots
      for j in range(i + 1, len(axes)):
          axes[j].axis('off')

      # Show or save the plot
      plt.suptitle(f"Top 20 in {niche} for {sender} signaling to {receiver} corr with GFR", fontsize=16, y=1.02)
      #plt.savefig("top20_gfr_correlations.png", dpi=300, bbox_inches='tight')
      plt.show()
      # Get top 20 features most correlated with GFR (excluding GFR itself)
      top_20_features = gfr_correlations_sorted.drop(labels='GFR').tail(5).index

      # Set up the plot grid
      n_cols = 5
      n_rows = (len(top_20_features) + n_cols - 1) // n_cols
      fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 3), constrained_layout=True)

      # Flatten axes for easier indexing
      axes = axes.flatten()

      # Plot each feature
      for i, feature in enumerate(top_20_features):
          ax = axes[i]
          sns.regplot(data=merged_df, x='GFR', y=feature, ax=ax, scatter_kws={'s': 10}, line_kws={"color": "red"})
          ax.set_title(f"{feature}\nr={gfr_correlations_sorted[feature]:.2f}")
          ax.set_xlabel('GFR')
          ax.set_ylabel('Interaction score')

      # Turn off any unused subplots
      for j in range(i + 1, len(axes)):
          axes[j].axis('off')

      # Show or save the plot
      plt.suptitle(f"Top 20 in {niche} for {sender} signaling to {receiver} corr with GFR", fontsize=16, y=1.02)
      #plt.savefig("top20_gfr_negative_correlations.png", dpi=300, bbox_inches='tight')
      plt.show()

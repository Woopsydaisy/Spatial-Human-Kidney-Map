import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from nichecompass.models import NicheCompass

# Set paths and parameters
fig_dir = "/home/bcd/revision_nature/nichecompass/revision_atlas/plots_script4"
model_dir = "/home/bcd/revision_nature/nichecompass/revision_atlas/model"
adata_file = "postnichecompass_adata.h5ad"
sc.settings.figdir = fig_dir
sc.set_figure_params(dpi_save=900)

# Load AnnData and model
adata = sc.read_h5ad(f"{model_dir}/{adata_file}")
model = NicheCompass.load(dir_path=model_dir, adata=None, adata_file_name=adata_file)

# Compute neighbors and UMAP using latent space
latent_key = "nichecompass_latent"
sc.pp.neighbors(model.adata, use_rep=latent_key, key_added=latent_key)
sc.tl.umap(model.adata, neighbors_key=latent_key)
adata.obsm['X_umap_nichecompass'] = model.adata.obsm['X_umap'].copy()

# Annotate niches and plot UMAPs
annotation_map = {
    0: 'PT Niche NC', 1: 'LOH Niche NC', 2: 'CNT_PC Niche NC', 3: 'Fibro_Immune_Vascular Niche NC',
    4: 'iPT Niche NC', 5: 'PT Niche NC', 6: 'DCT Niche NC', 7: 'iLOH Niche NC', 8: 'Glomerular Niche NC',
    9: 'CNT_PC Niche NC', 10: 'CNT_PC Niche NC', 11: 'LOH Niche NC', 12: 'Fibro_Immune_Vascular Niche NC',
    13: 'PT Niche NC'
}
adata.obs['nichecompass_niches'] = adata.obs[latent_cluster_key].astype(int).map(annotation_map).astype('category')
model.adata.obs['nichecompass_niches'] = adata.obs['nichecompass_niches']
model.adata.obs['niches_annotation_based'] = adata.obs['niches_annotation_based']

# UMAP plots
for color_key in ['nichecompass_niches', 'niches_annotation_based', 'annotation_updated', 'unique_sample_identifier']:
    sc.pl.umap(model.adata, color=color_key, title=f"{color_key} in Latent Space", save=f"_{color_key}.png")
    sc.pl.embedding(adata, basis='umap_nichecompass', color=color_key, title=f"{color_key} in Latent Space", save=f"_{color_key}_adata_only.png")

# Re-run GP testing using updated niche annotation
enriched_gps = model.run_differential_gp_tests(
    cat_key='nichecompass_niches',
    selected_cats=None,
    comparison_cats="rest",
    log_bayes_factor_thresh=log_bayes_thresh
)

# Heatmap
df = model.adata.obs[['nichecompass_niches'] + enriched_gps].groupby('nichecompass_niches').mean()
norm_df = pd.DataFrame(MinMaxScaler().fit_transform(df), columns=df.columns, index=df.index)
norm_df.to_csv(f"{fig_dir}/normalized_df_gp_across_niches_updated.csv")
sns.heatmap(norm_df, cmap='viridis')
plt.xticks(rotation=45, fontsize=8, ha="right")
plt.xlabel("Gene Programs", fontsize=16)
plt.savefig(f"{fig_dir}/enriched_gps_heatmap_updated_niche_key.png", bbox_inches="tight", dpi=900)

# Updated gene program summary
summary_df = gp_summary[gp_summary["gp_name"].isin(enriched_gps)].copy()
summary_df["gp_name"] = pd.Categorical(summary_df["gp_name"], categories=enriched_gps, ordered=True)
summary_df = summary_df.sort_values("gp_name")
summary_df.to_csv(f"{fig_dir}/gp_summary_updated_niche_key.csv", index=False)

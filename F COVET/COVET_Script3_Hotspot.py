import hotspot
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import mplscience
import os

# Load your preprocessed data (ensure this path is correct)
st_data_sst = sc.read_h5ad("/home/liranmao/katalin/envi/0_revision/envi_new_object/revision_object_for_colab_ENVI_envi_all.h5ad")  # Replace with actual file


def load_umap_numpy(st_data, filename="umap_covet_coordinates.npy"):
    """Load UMAP coordinates from NumPy array"""
    umap_coords = np.load(filename)
    st_data.obsm['UMAP_COVET'] = umap_coords
    print(f"UMAP coordinates loaded from {filename}")
    return st_data

# Load your preprocessed data (ensure this path is correct)
st_data_sst = load_umap_numpy(st_data_sst, filename = '/home/liranmao/katalin/envi/0_revision/envi_new_object/0_envi_16reso_results/umap_covet_coordinates_min_dist_0.1.npy')

st_data_sst = st_data_sst[st_data_sst.obs['Condition'].isin(['DKD', 'Control'])]

st_data_sst.X=st_data_sst.layers['counts'].copy() 

sc.pp.filter_genes(st_data_sst, min_cells=3)

# st_data_sst.X=st_data_sst.layers['counts'].copy() 
conditions = st_data_sst.obs['Condition'].unique()
nich = st_data_sst.obs['niches_annotation_based'].unique()

for sample_condition in conditions:
    for sample_niche in nich:
        dir_name = f"/home/liranmao/katalin/envi/0_revision/hotspot_new_object/v2_no_gene_filter/v2_hotspot_results_nich_latent_all_gene_{sample_condition}_{sample_niche}"
        final_output = f"{dir_name}/hotspot_local_correlations_with_gene.pdf"
        # Check if final output exists
        if os.path.exists(final_output):
            print(f"Skipping {sample_condition} {sample_niche} - already processed")
            continue
   
        try:
            # Create directory for this condition and niche
            
            os.makedirs(dir_name, exist_ok=True)
            
            # Subset data
            adata = st_data_sst[
                (st_data_sst.obs['Condition'] == sample_condition) & 
                (st_data_sst.obs['niches_annotation_based'] == sample_niche)
            ]
            
            if len(adata) > 0:
                    # Initialize Hotspot
                    hs = hotspot.Hotspot(
                        adata,
                        layer_key="counts",
                        model='bernoulli',
                        latent_obsm_key="UMAP_COVET"
                    )
                    
                    # Run analysis and save results
                    hs.create_knn_graph(weighted_graph=False, n_neighbors=300)
                    
                    hs_results = hs.compute_autocorrelations(jobs=4)
                    hs_results.to_csv(f"{dir_name}/hotspot_autocorrelations.csv")
                    
                    hs_genes = hs_results.index[hs_results.FDR < 0.05]
                    hs_genes.to_series().to_csv(f"{dir_name}/hotspot_significant_genes.csv")
                    
                    lcz = hs.compute_local_correlations(hs_genes, jobs=20)
                    lcz.to_csv(f"{dir_name}/hotspot_local_correlations.csv")
                    
                    
                    modules = hs.create_modules(min_gene_threshold=50, core_only=False, fdr_threshold=0.05)
                    modules.to_csv(f"{dir_name}/hotspot_modules.csv")

                    module_scores = hs.calculate_module_scores()
                    module_scores.to_csv(f"{dir_name}/hotspot_module_scores.csv")

                    # Save visualizations
                    plt.figure(figsize=(10, 12))
                    hs.plot_local_correlations()
                    plt.savefig(f"{dir_name}/hotspot_local_correlations.png", dpi=300, bbox_inches='tight')
                    plt.savefig(f"{dir_name}/hotspot_local_correlations.pdf", format='pdf', bbox_inches='tight')
                    plt.close()

                    plt.figure(figsize=(10, 12))
                    hs.plot_local_correlations(yticklabels=True)
                    plt.savefig(f"{dir_name}/hotspot_local_correlations_with_gene.png", dpi=300, bbox_inches='tight')
                    plt.savefig(f"{dir_name}/hotspot_local_correlations_with_gene.pdf", format='pdf', bbox_inches='tight')
                    plt.close()

                    print(f"Finished {sample_condition} {sample_niche}")
        except Exception as e:
            print(f"Error processing {sample_condition} {sample_niche}: {str(e)}")
            # Optionally log the full traceback for debugging
            import traceback
            print(traceback.format_exc())
            continue
       

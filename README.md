  # Spatial-Human-Kidney-Map

## A. Exporting Spatial Data

- **Exporting_CosMx_Script1.py**  
  Using exported data from AtoMx for the spatial datasets, this script generates an h5ad object of the entire slide for further usage.  
  *Container used:* `docker://10jll/spatial:version2`

- **Exporting_CosMx_Script2.py**  
  Partitions the slide data into individual samples and generates QC information.  
  *Container used:* `docker://10jll/spatial:version2`

- **Exporting_Xenium.py**  
  Exports Xenium spatial datasets from the raw output and generates QC information.  
  *Container used:* `docker://10jll/xenium:version4`

---

## B. Integration

### Integration_CosMx

- **Integration_CosMx_Script1.py**  
  Prepares all CosMx files for integration and concatenates with the reference SN atlas.  
  *Container used:* `docker://10jll/spatial:version2`

- **Integration_CosMx_Script2.py**  
  Runs SCVI and SCANVI to filter the CosMx dataset.  
  *Container used:* `docker://10jll/scvi_cuda12:version4`

### Integration_Xenium

- **Integration_Xenium_Script1.py**  
  Prepares all Xenium files for integration and concatenates with the reference SN atlas.  
  *Container used:* `docker://10jll/spatial:version2`

- **Integration_Xenium_Script2.py**  
  Runs SCVI and SCANVI to filter the Xenium dataset.  
  *Container used:* `docker://10jll/scvi_cuda12:version4`

### Final_Combined_Integration

- **Final_Combined_Integration_Script1_Concat.py**  
  Concatenates filtered CosMx, Xenium, and SN reference atlas.  
  *Container used:* `docker://10jll/spatial:version2`

- **Final_Combined_Integration_Script2_SCVI_SCANVI.py**  
  Runs final integration using SCVI and SCANVI.  
  *Container used:* `docker://10jll/scvi_cuda12:version4`

- **Final_Combined_Integration_Script3_Annotation.py**  
  Performs Leiden clustering and annotation.  
  *Container used:* `docker://10jll/spatial:version2`

- **Final_Combined_Integration_Script4_Imputation.py**  
  Imputes spatial dataset using SCANVI latent space.  
  *Container used:* `docker://10jll/spatial:version2`

---

## C. SCIB (after Luecken et al. 2021)

- **SCIB_Script1_Harmony_Scanorama.py**  
  Integrates spatial dataset with SN-Seq using Harmony and Scanorama.

- **SCIB_Script2_Harmony_Pyliger.py**  
  Integrates using Pyliger.

- **SCIB_Script3_scvi_scanvi.py**  
  Integrates using SCVI and SCANVI.

- **SCIB_Script4_benchmarking.py**  
  Calculates benchmarking parameters.

> These scripts were run on Google Colab.

---

## D. Niche Analysis

- **Niche_Analysis_Script1_call_neighbors.py**  
  Calls neighbors in CosMx and Xenium datasets.  
  *Container used:* `docker://10jll/spatial:version2`

- **Niche_Analysis_Script2_Neighbor_Dataframe.py**  
  Generates the neighbor dataframe.  
  *Container used:* `docker://10jll/spatial:version2`

- **Niche_Analysis_Script3_Kmeans.py**  
  K-means clustering of the neighbor dataframe.  
  *Container used:* `docker://10jll/spatial:version2`

- **Niche_Analysis_Script4_niche_annotation.py**  
  Niche annotation.  
  *Container used:* `docker://10jll/spatial:version2`

- **Niche_Analysis_Script5_cellphoneDB.py**  
  Prepares and runs CellPhoneDB using imputed expression and microenvironment key.  
  *Container used:* `docker://10jll/cellphonedb:version1`

- **Niche_Analysis_Script6_GFR_correlation.py**  
  Correlates CellPhoneDB interaction scores with disease severity in the glomerular niche.  
  *Container used:* `docker://10jll/spatial:version2`

- **Niche_Analysis_Script7_DEGs.py**  
  Calculates DEGs between DKD and control patients at niche level.  
  *Container used:* `docker://10jll/spatial:version2`

- **Niche_Analysis_Script8_spatial_genesignatures.py**  
  Generates pseudobulk counts and spatial gene signatures.  
  *Container used:* `docker://10jll/spatial:version2`

- **Niche_Analysis_Script9_DESeq2.R**  
  Runs DESeq2 on pseudobulk counts.  
  *Run on:* Google Colab

---

## E. NicheCompass

- **Nichecompass_Script1_setup.py**  
  Sets up model and computes adjacency matrices.  
  *Container used:* `docker://10jll/nichecompass:version1`

- **Nichecompass_Script2_model.py**  
  Runs the NicheCompass model.  
  *Container used:* `docker://10jll/nichecompass:version1`

- **Nichecompass_Script3_clustering_annotation.py**  
  Performs clustering, annotation, and gene program comparison.  
  *Container used:* `docker://10jll/nichecompass:version1`

---

## F. COVET

- **COVET_Script1.py**  
  Runs COVET on the CosMx dataset.

- **COVET_Script2_Visualization.py**  
  Downstream visualizations.

- **COVET_Script3_Hotspot.py**  
  Hotspot analysis based on latent niche representation.  
  *Run in local environment*

- **COVET_Script4_Hotspot_Visualization.py**  
  Visualizes gene-gene correlation matrices.  
  *Run on:* Google Colab

---

## G. Injured Tubular Microenvironments

- **Tubular_ME_Script1_Kmeans.py**  
  K-means clustering of neighbor dataframe.  
  *Container used:* `docker://10jll/spatial:version2`

- **Tubular_ME_Script2_annotation.py**  
  Annotation of microenvironments.  
  *Container used:* `docker://10jll/spatial:version2`

- **Tubular_ME_Script3_20um_Env.py**  
  Infers 20um neighborhood for each ME.  
  *Container used:* `docker://10jll/spatial:version2`

- **Tubular_ME_Script4_cellphoneDB.py**  
  Calculates ligand-receptor interactions using imputed expression.  
  *Container used:* `docker://10jll/cellphonedb:version1`

- **Tubular_ME_Script5_correlation.py**  
  Correlates iPT–Fibroblast interaction scores with disease severity.  
  *Container used:* `docker://10jll/spatial:version2`

- **Tubular_ME_Script6_cellfractions.py**  
  Statistical comparison of cell types across injured tubular MEs.  
  *Container used:* `docker://10jll/spatial:version2`

---

## H. Immune Cell Atlas

- **Script1_Prepare_SCVI_Xenium.py**  
  Prepares Xenium and SN-Seq input for SCVI/SCANVI.  
  *Container used:* `docker://10jll/spatial:version2`

- **Script2_Xenium_scvi_scanvi.py**  
  Runs SCVI/SCANVI on Xenium and SN-Seq.  
  *Container used:* `docker://10jll/scvi_cuda12:version4`

- **Script3_Add_CosMx.py**  
  Adds CosMx to integrated Xenium and SN for SCVI/SCANVI.  
  *Container used:* `docker://10jll/spatial:version2`

- **Script4_scvi_scanvi.py**  
  Runs SCVI/SCANVI on CosMx, Xenium, and SN-Seq.  
  *Container used:* `docker://10jll/scvi_cuda12:version4`

- **Script5_Annotation_Imputation.py**  
  Performs clustering, annotation, and gene expression imputation from SN-Seq.  
  *Container used:* `docker://10jll/spatial:version2`

---

## I. Immune Microenvironments

- **Immune_ME_Script1_call_neighbors.py**  
  Generates neighbor dataframe using immune cell annotations.  
  *Container used:* `docker://10jll/spatial:version2`

- **Immune_ME_Script2_Kmeans.py**  
  Clusters immune cells into microenvironments.  
  *Container used:* `docker://10jll/spatial:version2`

- **Immune_ME_Script3_20um_Env.py**  
  Infers 20um environments.  
  *Container used:* `docker://10jll/spatial:version2`

- **Immune_ME_Script4_prepare_cellphoneDB.py**  
  Prepares CellPhoneDB with imputed expression.  
  *Container used:* `docker://10jll/spatial:version2`

- **Immune_ME_Script5_cellphoneDB.py**  
  Runs CellPhoneDB with a microenvironment key.  
  *Container used:* `docker://10jll/cellphonedb:version1`

- **Immune_ME_Script6_correlation.py**  
  Correlates interaction scores in immune MEs with disease severity.  
  *Container used:* `docker://10jll/spatial:version2`

- **Immune_ME_Script7_call_neighbors.py**  
  Visualizes cell neighborhoods and compares immune cell fractions between glomerular and tubular immune MEs.  
  *Container used:* `docker://10jll/spatial:version2`


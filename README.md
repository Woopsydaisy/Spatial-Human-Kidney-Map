# Spatial-Human-Kidney-Map
A Exporting Spatial Data:
  Exporting_CosMx_Script1.py: Using exported data from AtoMx for the spatial datasets, this script generates an h5ad object of the entire slide for 
  further usage. Container used: docker://10jll/spatial:version2
  Exporting_CosMx_Script2.py: This script partitions the slide data from the above script into individual samples and generates QC information. 
  Container used: docker://10jll/spatial:version2
  Exporting_Xenium.py: This script exports Xenium spatial datasets from the raw Xenium output and generates QC information.
  Container used: docker://10jll/xenium:version4
B Integration
  Integration_CosMx
    Integration_CosMx_Script1.py: This prepares all CosMx files for integration and concats with reference SN atlas. 
    Container used: docker://10jll/spatial:version2
    Integration_CosMx_Script2.py: This shows how to run through scvi and scanvi to filter the CosMx dataset. 
    Container used: docker://10jll/scvi_cuda12:version4
  Integration_Xenium
    Integration_Xenium_Script1.py: This prepares all Xenium files for integration and concats with reference SN atlas. 
    Container used: docker://10jll/spatial:version2
    Integration_Xenium_Script2.py: This shows how to run through scvi and scanvi to filter the Xenium dataset. 
    Container used: docker://10jll/scvi_cuda12:version4
  Final_Combined_Integration
    Final_Combined_Integration_Script1_Concat.py: This concats filtered CosMx, Xenium with the SN reference atlas.
    Container used: docker://10jll/spatial:version2
    Final_Combined_Integration_Script2_SCVI_SCANVI.py: This shows the integration using SCVI and SCANVI of the datasets. 
    Container used: docker://10jll/scvi_cuda12:version4
    Final_Combined_Integration_Script3_Annotation.py: This shows leidenclustering and subsequent annotation. 
    Container used: docker://10jll/spatial:version2
    Final_Combined_Integration_Script4_Imputation.py: This shows imputation of spatial dataset using the scanvi latent space. 
    Container used: docker://10jll/spatial:version2
C SCIB (after Luecken et al. 2021)
  SCIB_Script1_Harmony_Scanorama.py: This script shows the integration of spatial dataset with SN-Seq using Harmony and Scanorama. 
  SCIB_Script2_Harmony_Pyliger.py: This script shows the integration of spatial dataset with SN-Seq using Pyliger.
  SCIB_Script3_scvi_scanvi.py: This script shows the integration of spatial dataset with SN-Seq using scvi and scanvi.
  SCIB_Script4_benchmarking.py: This script shows the calculation of benchmarking parameters.
  These scripts were run on Google Colab.
D Niche_Analysis
  Niche_Analysis_Script1_call_neighbors.py: This script shows how to call neighbors in the CosMx and Xenium datasets. 
  Container used: docker://10jll/spatial:version2
  Niche_Analysis_Script2_Neighbor_Dataframe.py: This script shows how to generate the neighbordataframe for subsequent clustering. 
  Container used: docker://10jll/spatial:version2
  Niche_Analysis_Script3_Kmeans.py: Kmeans clustering of the neighbordataframe. 
  Container used: docker://10jll/spatial:version2
  Niche_Analysis_Script4_niche_annotation.py: Niche Annotation. 
  Container used: docker://10jll/spatial:version2
  Niche_Analysis_Script5_cellphoneDB.py: Using imputed geneexpression this script prepares and runs cellphoneDB with a Microenvironment key. 
  Container used: docker://10jll/cellphonedb:version1
  Niche_Analysis_Script6_GFR_correlation.py: This script shows how to correlate the cellphone DB interaction scores with disease severity in the glomerular niche between glomerular endothelial cells and podocytes. 
  Container used: docker://10jll/spatial:version2
  Niche_Analysis_Script7_DEGs.py: This script shows how DEGs between DKD and control patients on a niche level with imputed gene expression were calculated. 
  Container used: docker://10jll/spatial:version2
  Niche_Analysis_Script8_spatial_genesignatures.py: This script shows how to generate pseudobulk counts and extract the spatial genesignature. 
  Container used: docker://10jll/spatial:version2
  Niche_Analysis_Script9_DESeq2.R: This script shows how to analyze pseudobulk counts using DESeq2. 
  This was run on Google Colab. 
E Nichecompass
  Nichecompass_Script1_setup.py: This script shows how the Nichecompassmodel is setup and adjacency matrices are calculated.
  Container used: docker://10jll/nichecompass:version1
  Nichecompass_Script2_model.py: This script shows how the Nichecompass is run.
  Container used: docker://10jll/nichecompass:version1
  Nichecompass_Script1_setup.py: This script shows downstream clustering, annotation and comparison of geneprograms.
  Container used: docker://10jll/nichecompass:version1
F COVET
  COVET_Script1.py: This shows how COVET was run on the CosMx dataset. 
  COVET_Script2_Visualization.py: This shows the downstream visualization.
  COVET_Script3_Hotspot.py: This shows the Hotspot analysis based on the latent niche representation.
  These Scripts were run in a local environment. 
  COVET_Script4_Hotspot_Visualization.py: This shows how gene-gene correlation matrices were visualized. 
  This script was run on Google Colab. 
G Injured_Tubular_Microenvironments:
  Tubular_ME_Script1_Kmeans.py: This Script shows how the neighbordataframe generated in the nicheanalysis under E was clustered.
  Container used: docker://10jll/spatial:version2
  Tubular_ME_Script2_annotation.py: This Script shows the annotation of Microenvironments.
  Container used: docker://10jll/spatial:version2
  Tubular_ME_Script3_20um_Env.py: This Script shows how the 20um environment for each tubular ME was inferred.
  Container used: docker://10jll/spatial:version2
  Tubular_ME_Script4_Kmeans.py: This Script shows how using imputed geneexpression and a microenvironment key ligand receptor interactions were calculated.
  Container used: docker://10jll/cellphonedb:version1
  Tubular_ME_Script5_correlation.py: This Script shows how the resulting interaction scores of interactions between iPT and Fibroblasts were correlated with disease severity in the profibrotic microenvironment.
  Container used: docker://10jll/spatial:version2
  Tubular_ME_Script6_cellfractions.py: This Script shows the statistical comparison of key celltypes between injured tubular microenvironments. 
  Container used: docker://10jll/spatial:version2
H Immune_Cell_Atlas:
  Script1_Prepare_SCVI_Xenium.py: Preparing input for subsequent scvi, scanvi of xenium and SN-Seq. 
  Container used: docker://10jll/spatial:version2
  Script2_Xenium_scvi_scanvi.py: Scvi, scanvi of xenium and SN-Seq. 
  Container used: docker://10jll/scvi_cuda12:version4
  Script3_Add_CosMx.py: Adding CosMx to integrated Xenium and SN for subsequent scvi, scanvi of CosMx, Xenium and SN-Seq. 
  Container used: docker://10jll/spatial:version2
  Script4_scvi_scanvi.py: Scvi, scanvi of CosMx, Xenium and SN-Seq. 
  Container used: docker://10jll/scvi_cuda12:version4
  Script5_Annotation_Imputation.py: This script shows the clustering and subsequent annotation of the immune cell atlas and imputation fo gene expression from SN-Seq to the spatial dataset. 
  Container used: docker://10jll/spatial:version2
I Immune_Microenvironments
  Immune_ME_Script1_call_neighbors.py: This Script shows the generation of the neighbordataframe using the immune cell annotations from the immune cell atlas. 
  Container used: docker://10jll/spatial:version2
  Immune_ME_Script2_Kmeans.py: This Script shows the clustering into microenvironments. 
  Container used: docker://10jll/spatial:version2
  Immune_ME_Script3_20um_Env.py: This Script shows how the neighboring 20um cells were inferred. 
  Container used: docker://10jll/spatial:version2
  Immune_ME_Script4_prepare_cellphoneDB.py: This Script shows the preparation of cellphone DB using imputed geneexpression. 
  Container used: docker://10jll/spatial:version2
  Immune_ME_Script5_cellphoneDB.py: This Script shows running of cellphoneDB with a microenvironment key. 
  Container used: docker://10jll/cellphonedb:version1
  Immune_ME_Script6_correlation.py: This Script shows how the resulting interaction scores of interactions in the injured tubular immune ME and glomerular immune ME were correlated with disease severity. 
  Container used: docker://10jll/spatial:version2
  Immune_ME_Script7_call_neighbors.py: This Script shows the generation of pie charts to visualize the center cell and 20um environment as well as the comparison of cell fractions of immune cell types in glomerulus and injured tubular immune ME. 
  Container used: docker://10jll/spatial:version2

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

# Read in counts and metadata
counts <- read.csv("/content/pseudobulk_counts_niches.csv", row.names = 1)
metadata <- read.csv("/content/pseudobulk_metadata_niches.csv", row.names = 1)

# Make sure metadata rownames match count column names
all(rownames(metadata) == rownames(counts))  # should return TRUE


counts_t <- t(round(as.matrix(counts)))
storage.mode(counts_t) <- "integer"

any(counts < 0)

# Ensure the sample names match
stopifnot(all(colnames(counts_t) == rownames(metadata)))

dim(counts_t)     # should return (genes x samples)
dim(metadata)

dds_list <- list()
res_list <- list()

# Get all niche levels
niche_levels <- unique(metadata$group)

niche_levels

# Loop through each niche
for (niche in niche_levels) {
  # Create a binary group column
  metadata$binary_group <- ifelse(metadata$group == niche, niche, paste0("non_", niche))
  metadata$binary_group <- factor(metadata$binary_group)

  # Re-create DESeq object
  dds_tmp <- DESeqDataSetFromMatrix(countData = t(counts),  # genes as rows
                                    colData = metadata,
                                    design = ~ binary_group)

  dds_tmp <- estimateSizeFactors(dds_tmp, type = "poscounts")
  dds_tmp <- DESeq(dds_tmp)
  # Run results for this niche vs rest
  res <- results(dds_tmp, contrast = c("binary_group", niche, paste0("non_", niche)))
  res_ordered <- res[order(res$padj), ]

  # Store
  dds_list[[niche]] <- dds_tmp
  res_list[[niche]] <- res_ordered
}


for (niche in names(res_list)) {
  write.csv(as.data.frame(res_list[[niche]]), file = paste0("DE_niche_vs_rest_", niche, ".csv"))
}

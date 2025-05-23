# This script performs differential gene expression analysis (DGEA) on a single-cell RNA-seq dataset of ageing Drosophila, 
# comparing gene expression between 5-day-old flies and older age groups within each sex. 
# It uses the MAST method via Seurat, with sex_age as the grouping variable and batch as a latent variable. 
# For each cell type in the dataset, it tests all valid comparisons where a 5-day group (male or female) and another group are both present. 
# Results for each comparison are appended to combined_markers.csv, including Benjamini-Hochberg adjusted p-values.

# Load libraries
library(Seurat)
library(SeuratDisk)
library(reticulate)
library(zellkonverter)
library(dplyr)
library(MAST)

# Set the working directory
setwd("~/projects/ageing_flies")

# Load the anndata object
ad <- readH5AD("Data/adata_body_filtered_2.h5ad")

# Convert the anndata object to a Seurat object
adata_Seurat <- as.Seurat(ad, counts = "X", data = NULL)

# Reformat the metadata 
adata_Seurat$age <- as.factor(adata_Seurat$age)
adata_Seurat$batch <- as.factor(adata_Seurat$batch)
adata_Seurat$sex_age <- as.factor(adata_Seurat$sex_age)
adata_Seurat$afca_annotation <- as.factor(adata_Seurat$afca_annotation)

# Define groups
young_groups <- c("male_5", "female_5")
all_groups <- unique(adata_Seurat$sex_age)
other_groups <- setdiff(all_groups, young_groups)

# Get unique cell types (as character vector)
cell_types <- unique(as.character(adata_Seurat$afca_annotation))

# Set output file path
results_path <- "Results/combined_markers.csv"


# Begin looping over cell types and comparisons
for (ct in cell_types) {
  cell_subset <-  subset(adata_Seurat, subset = afca_annotation == as.character(ct))
  cat(paste("\n Covering the cell type:", ct, "\n"))
  
  groups_present <- unique(cell_subset$sex_age)
  
  for (young in young_groups) {
    if (young %in% groups_present) {
      for (grp in other_groups) {
        if (grp %in% groups_present) {
          cat(paste("  Running DGEA for:", young, "vs", grp, "\n"))
          try({
            markers <- FindMarkers(cell_subset,
                                   ident.1 = young,
                                   ident.2 = grp,
                                   group.by = "sex_age",
                                   test.use = "MAST",
                                   latent.vars = "batch")
            
            markers$cell_type <- ct
            markers$comparison <- paste(young, "vs", grp)
            markers$gene <- rownames(markers)
            markers$adjusted_p_val <- p.adjust(markers$p_val, method = "BH")
            
            # Write to CSV file (with header if first write)
            write.table(markers, results_path,
                        sep = ",",
                        row.names = FALSE,
                        col.names = !file.exists(results_path),
                        append = file.exists(results_path))
          }, silent = TRUE)
        }
      }
    }
  }
}

cat("All results saved incrementally to:\n", results_path, "\n")

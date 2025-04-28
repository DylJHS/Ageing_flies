#!/usr/bin/env Rscript

library(Seurat)
library(zellkonverter)
library(dplyr)
library(fs)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: run_dgea_single.R <input_file> <output_dir>")
}
input_file <- args[1]
output_dir <- args[2]

# Check input file
if (!file.exists(input_file)) {
  stop(paste("Input file does not exist:", input_file))
}

# Setup output file
file_name <- path_file(input_file)
ct_name <- path_ext_remove(file_name)
results_file <- file.path(output_dir, paste0("markers_", ct_name, ".csv"))

# Read and convert
suppressWarnings({
  ad <- readH5AD(input_file)
  adata_Seurat <- as.Seurat(ad, counts = "X", data = NULL)
})
rm(ad)

# Convert variables
adata_Seurat$age <- as.factor(adata_Seurat$age)
adata_Seurat$sex_age <- as.factor(adata_Seurat$sex_age)
adata_Seurat$afca_annotation <- as.factor(adata_Seurat$afca_annotation)

# Define groups
young_group <- "5"
groups_present <- unique(adata_Seurat$age)
other_groups <- setdiff(groups_present, young_group)

# Run DGEA
res_list <- list()
if (young_group %in% groups_present) {
  for (grp in other_groups) {
    if (grp %in% groups_present) {
      cat(paste("  Running DGEA for:", young_group, "vs", grp, "\n"))
      markers <- FindMarkers(adata_Seurat,
                             ident.1 = young_group,
                             ident.2 = grp,
                             group.by = "age",
                             test.use = "wilcox")

      markers$cell_type <- ct_name
      markers$comparison <- paste(young_group, "vs", grp)
      markers$gene <- rownames(markers)
      markers <- markers %>% rename(p_val_BF = p_val_adj)
      markers$p_val_adj <- p.adjust(markers$p_val, method = "BH")

      res_list[[grp]] <- markers
    }
  }
}

# Combine and write
if (length(res_list) > 0) {
  all_results <- do.call(rbind, res_list)
  write.table(all_results, results_file, sep = ",", row.names = FALSE, col.names = TRUE)
}

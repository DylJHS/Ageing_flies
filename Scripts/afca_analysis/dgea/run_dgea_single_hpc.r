#!/usr/bin/env Rscript

# This script performs differential gene expression analysis (DGEA) for a single cell type using the pre-filtered afca .h5ad file as input. 
# It is designed to be run on an HPC system as part of an array job, where each task processes one cell type. 
# The script converts the input AnnData file to a Seurat object, sets relevant metadata factors, and uses Seurat’s FindMarkers with the 
# MAST test to compare 5-day-old flies against older age groups, adjusting for individual variation using indiv as a latent variable.
# The results for each comparison are saved as a .csv file named after the input file. 
# This script is intended to be called by run_dgea_array.sh, which distributes the workload across multiple SLURM array jobs.

# Load required libraries
library(Seurat)
library(zellkonverter)
library(dplyr)
library(fs)
library(MAST)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: run_dgea_single_hpc.r <input_file> <output_dir>")
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
cat(paste("\t\tRunning DGEA on:", ct_name, "\n"))

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
adata_Seurat$indiv <- as.factor(adata_Seurat$indiv)

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
      markers <- FindMarkers(
        adata_Seurat,
        ident.1 = young_group,
        ident.2 = grp,
        group.by = "age",
        test.use = "MAST", # Using MAST to account for batch effects
        latent.vars = "indiv" # Adjust for individual variability
      ) 
      
      # Check markers
      cat("Shape of markers:", paste(dim(markers), collapse = " x "), "\n")
      # print the first few rows of markers
      print(head(markers))
      
      
      markers$cell_type <- ct_name
      markers$comparison <- paste(young_group, "vs", grp)
      markers$gene <- rownames(markers)
      # markers <- markers %>% rename(p_val_BF = p_val_adj)
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
cat("\n\nFinished Running DGEA on Cell Type !\n")
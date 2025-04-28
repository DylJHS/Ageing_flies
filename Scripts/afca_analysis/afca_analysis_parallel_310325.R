# Load libraries
library(anndata)
library(SeuratData)
library(Seurat)
library(SeuratDisk)
library(reticulate)
library(zellkonverter)
library(dplyr)
library(MAST)
library(future)
library(future.apply)

# Load the anndata object
ad <- readH5AD("../Data/adata_body_filtered_2.h5ad")

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

# Pre-split the main Seurat object by cell type to avoid exporting the full object each time.
cell_type_list <- lapply(cell_types, function(ct) {
  subset(adata_Seurat, subset = as.character(afca_annotation) == ct)
})
names(cell_type_list) <- cell_types

# Set up parallel backend
plan(multisession, workers = 4)  # Adjust 'workers' as needed

# Define a function to run DGEA for a given cell type subset object
run_dgea_for_subset <- function(cell_subset) {
  # Determine cell type label from the subset (should be unique)
  ct <- unique(as.character(cell_subset$afca_annotation))[1]
  cat(paste("\nCovering the cell type:", ct, "\n"))
  
  # Determine which sex_age groups are present in this subset
  groups_present <- unique(cell_subset$sex_age)
  
  # Container for marker results for this cell type
  ct_markers <- list()
  
  # Loop over each young group
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
            
            # Append additional information to the results
            markers$cell_type <- ct
            markers$comparison <- paste(young, "vs", grp)
            markers$gene <- rownames(markers)
            
            # Use a unique key for each comparison
            key <- paste(ct, young, grp, sep = "_")
            ct_markers[[key]] <- markers
          }, silent = TRUE)
        }
      }
    }
  }
  
  return(ct_markers)
}

# Run the DGEA in parallel over the pre-split list of cell type objects
results_list <- future_lapply(cell_type_list, run_dgea_for_subset)

# Flatten the nested list into one list of data frames
all_markers <- do.call(c, results_list)

# Combine all results into a single data frame
combined_markers <- bind_rows(all_markers)

# Apply global p-value adjustment across all comparisons using BH correction (on raw p-values)
combined_markers$adjusted_p_val <- p.adjust(combined_markers$p_val, method = "BH")

# Print first few rows of the results
print(head(combined_markers))

# Save the results to a CSV file
write.csv(combined_markers,
          file = "/hpc/shared/onco_janssen/dhaynessimmons/results/ageing_flies/combined_markers.csv",
          row.names = FALSE)

cat("Results saved to combined_markers.csv\n")

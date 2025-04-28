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
adata_Seurat$afca_annotation <- as.factor(adata_Seurat$afca_annotation)

# Define groups
young_group <- "5"
all_groups <- unique(adata_Seurat$age)
other_groups <- setdiff(all_groups, young_group)

# Get unique cell types (as character vector)
cell_types <- unique(as.character(adata_Seurat$afca_annotation))

# Set output file path
results_path <- "Results/combined_markers.csv"


# Begin looping over cell types and comparisons
for (ct in cell_types) {
  cell_subset <-  subset(adata_Seurat, subset = afca_annotation == as.character(ct))
  cat(paste("\n Covering the cell type:", ct, "\n"))
  
  groups_present <- unique(cell_subset$age)
  
  if (young_group %in% groups_present) {
    for (grp in other_groups) {
      if (grp %in% groups_present) {
        cat(paste("  Running DGEA for:", young_group, "vs", grp, "\n"))
        try({
          markers <- FindMarkers(cell_subset,
                                 ident.1 = young_group,
                                 ident.2 = grp,
                                 group.by = "age",
                                 test.use = "MAST",
                                 latent.vars = c("batch", "sex"))
          
          markers$cell_type <- ct
          markers$comparison <- paste(young_group, "vs", grp)
          markers$gene <- rownames(markers)
          
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

# update the df by adjusting the pvalues with BH correction
combined_markers <- read.csv( "Results/combined_markers.csv") %>%
  group_by(cell_type) %>%
  mutate(adj_p_val_within_cell_type = p.adjust(p_val, method = "BH")) %>%
  ungroup()

# Save the final version of the df
write.csv(combined_markers, "Results/full_combined_markers.csv", row.names = FALSE)


cat("All results saved incrementally to:\n", results_path, "\n")

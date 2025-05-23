---
title: "afca_analysis_visualisation"
description: This script creates volcano plots from differential gene expression results
  across ageing comparisons in the enterocyte grouping in the afca data which was
  preprocessed by afca_analysis_entero_geneset_data_prep_160525.ipynb. The plot will
  show the difference of the young (5day old) group versus the other
date: "2025-05-17"
output:
  pdf_document: default
  html_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "D:/projects/ageing_flies")
```

Load in the libraries
```{r}
library(dplyr)
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)
library(patchwork)
library(cowplot)
```

Load in the data and set paths
```{r}
# date <- format(Sys.Date(), "%Y%m%d")
date <- "20250517"
input_path <- file.path("Results/age_scdgea", date, "entero_markers.csv")
fig_path <- paste0("Figures/afca_analysis/", date)
entero_markers <- read.csv(input_path, header = TRUE, stringsAsFactors = FALSE)
```

Define the variables
```{r}
genes_of_interest <- c("Su(var)205","Su(var)3-9","G9a", "HP1b", "HP1c", "HP4", "HP5", "HP6", "ADD1", "Su(var)2-HP2", "Su(var)3-7", "Lam", "LamC", "LBR", "Kdm4A", "Kdm4B", "His2Av", "His3.3A", "His3.3B")
```

Create the volcano plot for the group
```{r}
age_comps <- unique(entero_markers$comparison)
plots <- list()

for (comp in age_comps) {
  cat(paste0("\n\n age comparison: ", comp, "\n"))
  subset_comp <- subset(entero_markers, entero_markers$comparison == comp)
  
  # Order data so genes of interest come last (plotted on top)
  goi_flag <- subset_comp$gene %in% genes_of_interest
  subset_comp <- rbind(
    subset_comp[!goi_flag, ],
    subset_comp[goi_flag, ]
  )
  
  # Map colours: red for genes of interest, grey otherwise
  keyvals.colour <- ifelse(subset_comp$gene %in% genes_of_interest, '#920000',
    'grey')
  
  names(keyvals.colour) <- ifelse(subset_comp$gene %in% genes_of_interest, 'Gene of Interest',
    'Other')
  
  # Map sizes: 4 for genes of interest, 1 for others
  keyvals.size <- ifelse(
    subset_comp$gene %in% genes_of_interest, 4, 1
  )
  names(keyvals.size) <- ifelse(
    subset_comp$gene %in% genes_of_interest, 'Gene of Interest', 'Other'
  )
  
  
  goi_sub <- subset_comp %>%
    subset(., subset_comp$gene %in% genes_of_interest) %>%
    mutate(abs_FC = abs(avg_log2FC))
  
  # Find the max pval from the set
  min_pval <- min(goi_sub$p_val_adj)
  plt_min <- ifelse(min_pval < 0.005, 80/100*(min_pval), 0.005)
  
   # Find the max FC from the set
  max_FC <- max(goi_sub$abs_FC)
  plt_FC <- ifelse(max_FC > 2.5, 150/100*(max_FC), 2.5)
  
  # Significance of the comparisons
  both_insign <- all(goi_sub$abs_FC<1 | goi_sub$p_val_adj>0.1)
  signif <- ifelse(both_insign, "Not significant", "Some significant")
  
  cat(paste0(
    "\n All points insignif: ", both_insign,
    "\n Both conditions significant: ", signif
  ))
  
  vol_plt <- EnhancedVolcano(
    subset_comp,
    lab = ifelse(subset_comp$gene %in% genes_of_interest, subset_comp$gene, ""),
    selectLab = genes_of_interest,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    subtitle = signif,
    pCutoff = 0.1,
    xlim = (c(-plt_FC, plt_FC)),
    ylim = c(0, -log10(plt_min)),
    title = comp,
    colCustom = keyvals.colour,
    pointSize = keyvals.size,
    legendPosition = "none"
  )
  
  plots[[comp]] <- vol_plt
}

combined_plot <- wrap_plots(plots) + 
  plot_annotation(title = paste("Enterocyte volcano plot"))

if (!dir.exists(fig_path)) {
  dir.create(fig_path, recursive = TRUE)
}

file_name <- paste0(fig_path, "/volcano_enterocyte.png")
ggsave(file_name, plot = combined_plot, width = 20, height = 9, dpi = 300)
```


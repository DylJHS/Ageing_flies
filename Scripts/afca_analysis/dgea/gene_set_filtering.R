# Script to filter the gene set of interest from the larger afca entero marker
# dataset and save it to a new CSV file

library(dplyr)
library(tidyr)

print(getwd())
data <- read.csv("D:/projects/ageing_flies/Results/age_scdgea/20250517/entero_markers.csv")

genes_of_interest <- c(
  "Su(var)205","Su(var)3-9","G9a", "HP1b", 
  "HP1c", "HP4", "HP5", "HP6", "ADD1", "Su(var)2-HP2", 
  "Su(var)3-7", "Lam", "LamC", "LBR", "Kdm4A", 
  "Kdm4B", "His2Av", "His3.3A", "His3.3B"
)

reduced_gene_set_markers <- data %>% 
  filter(gene %in% genes_of_interest) 

#save 
write.csv(
  reduced_gene_set_markers, 
  "D:/projects/ageing_flies/Results/age_scdgea/20250517/reduced_gene_set_markers.csv", 
  row.names = FALSE
)


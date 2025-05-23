library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(writexl)

setwd("P:/Janssen/chet_repair/dhaynessimmons/data_analysis/Ageing flies") 

raw_tracks <- import("Data/ISC_chromatin_states.updt.bed") %>% as.data.frame()
trks_gr <- GRanges(raw_tracks)

# Split by chromosome and track name
gr_list <- split(trks_gr, interaction(seqnames(trks_gr), mcols(trks_gr)$name))

# Reduce overlapping or adjacent regions per group
# reduced <- x[, as.data.table(reduce(IRanges(start, end))), by = .(chrom, hgnc)]
reduced_gr <- lapply(gr_list, GenomicRanges::reduce)

# Calculate total span per track on each chromosome
total_span <- sapply(reduced_gr, function(gr) sum(width(gr)))

# Convert to a readable data frame
result_df <- data.frame(Chromosome_Track = names(total_span), Total_Span = total_span)


# Write to an excel file
write_xlsx(result_df, "total_span.xlsx")
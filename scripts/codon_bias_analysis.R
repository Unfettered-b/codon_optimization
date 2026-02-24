#!/usr/bin/env Rscript

#############################################
# Codon Bias Analysis using coRdon
# Genome-wide ENC + GC3 + Wright plot
#############################################

suppressPackageStartupMessages({
  library(coRdon)
  library(ggplot2)
})

# Snakemake inputs/outputs
input_fasta <- snakemake@input[[1]]
output_metrics <- snakemake@output[["metrics"]]
output_wright <- snakemake@output[["wright"]]

#############################################
# 1. Load CDS
#############################################

cds <- readSet(file = input_fasta)

#############################################
# 2. Compute ENC (Wright 1990)
#############################################
head(cds)
codon_files <- codonTable(cds)
milc <- MILC(codon_files)
head(milc)

enc_values <- ENC(codon_files)

#############################################
# 3. Compute GC3
#############################################

#############################################
# Custom GC3 function
#############################################

calculate_GC3 <- function(cds_object) {
  
  # Extract sequences as character vector
  seqs <- cds_object
  
  gc3_values <- sapply(seqs, function(seq) {
    
    seq <- toupper(seq)
    
    # Trim to full codons
    len <- nchar(seq)
    seq <- substr(seq, 1, len - (len %% 3))
    
    # Split into codons
    codons <- substring(seq,
                         seq(1, nchar(seq), 3),
                         seq(3, nchar(seq), 3))
    
    if (length(codons) == 0) return(NA)
    
    third_positions <- substr(codons, 3, 3)
    
    gc_count <- sum(third_positions %in% c("G", "C"))
    
    gc_count / length(third_positions)
  })
  
  return(gc3_values)
}

gc3_values <- calculate_GC3(cds)

#############################################
# 4. Combine into dataframe
#############################################

gene_ids <- names(cds)

metrics_df <- data.frame(
  gene_id = gene_ids,
  ENC = enc_values,
  GC3 = gc3_values
)

#############################################
# 5. Write gene metrics table
#############################################

write.table(
  metrics_df,
  file = output_metrics,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

#############################################
# 6. Wright Plot (ENC vs GC3)
#############################################

gc3_curve <- seq(0.01, 0.99, by = 0.01)

enc_expected <- 2 + gc3_curve +
  29 / (gc3_curve^2 + (1 - gc3_curve)^2)

wright_df <- data.frame(
  GC3 = metrics_df$GC3,
  ENC = metrics_df$ENC
)

curve_df <- data.frame(
  GC3 = gc3_curve,
  ENC = enc_expected
)

p <- ggplot(wright_df, aes(x = GC3, y = ENC)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_line(data = curve_df, aes(x = GC3, y = ENC),
            color = "red", linewidth = 1) +
  labs(
    title = "Wright Plot (ENC vs GC3)",
    x = "GC3",
    y = "ENC"
  ) +
  theme_minimal()

ggsave(output_wright, plot = p, width = 6, height = 5)

#############################################
# End
#############################################
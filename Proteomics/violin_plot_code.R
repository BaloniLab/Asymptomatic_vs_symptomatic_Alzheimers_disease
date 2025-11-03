# Required libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(tibble)

# Reading the files
setwd("/depot/pbaloni/data/Lab_members/Purba_Mandal/Proteomics/Brain_regions/BA9(Banner+ROSMAP)")
expr <- read.csv("Aggregated_expression_matrix.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("Final_metadata.csv")  

# Proteins of interest
proteins <- c("CPT1A", "CPT2", "ACADSB", "SDHC", "SUCLA2", "IDH2", "ECHS1", "BCKDHB","BCAT1", 
           "MRPL37", "TOMM40", "SLC25A4", "SLC25A6", "MCU", "VDAC2", "PDHB", "FIS1", "IMMT")
# Melt expression matrix
expr_long <- as.data.frame(expr[proteins, ]) %>%
  rownames_to_column("proteins") %>%
  pivot_longer(-proteins, names_to = "sample_id", values_to = "expr")

# Joining with metadata
expr_merged <- expr_long %>%
  left_join(metadata, by = "sample_id") 
expr_merged$diagnosis <- factor(expr_merged$diagnosis, levels = c("Control", "AD", "AsymAD"))

# Plotting function
plot_violin_with_stats <- function(proteins_name) {
  df <- expr_merged %>% filter(proteins == proteins_name)
  
  ggviolin(df, x = "diagnosis", y = "expr", fill = "diagnosis", 
           palette = "aaas", add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(comparisons = list(
      c("Control", "AD"),
      c("Control", "AsymAD"),
      c("AsymAD", "AD")),
      method = "t.test",
      p.adjust.method = "BH",                    # BH correction across the 3 pairwise tests
      label = "p.signif",                       
      size = 4) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold", colour = "black"),
      axis.text = element_text(size = 12, face = "bold", colour = "black"),
      axis.text.x =element_text(size = 12, face = "bold", colour = "black", angle = 45),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10, face = "bold")
    ) +
    labs(title = proteins_name, y = "Residual Expression", x = "Diagnosis") +
    theme(legend.position = "none")
}

library(ggpubr)

# Plotting all proteins
plots <- lapply(proteins, plot_violin_with_stats)

# Pplots in a grid and assign to a variable
combined_plot <- ggarrange(plotlist = plots, ncol = 3, nrow = 2)
combined_plot
# Saving the combined plot
ggsave("violin_plot_fission.png", plot = combined_plot, dpi = 600, bg = "transparent", limitsize = FALSE)


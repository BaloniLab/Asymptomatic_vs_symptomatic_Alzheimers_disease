# Cell type composition comparision using marker scores across diagnosis

# Load packages
library(ggplot2)
library(tidyr)
library(ggpubr)

# Taking marker scores df file from De_limma_cell_type_correction.R code
# First, convert to long format for ggplot
marker_scores_long <- marker_scores_df %>%
  pivot_longer(cols = c("Astrocytes", "Microglia", "Neuron", "Oligodendrocytes", "Endothelia"),
               names_to = "cell_type", values_to = "score")

# Custom color palette
custom_palette <- c("Control" = "lightblue", "AD" = "tomato", "AsymAD" = "violet")

# Define pairwise comparisons
comparisons <- list(
  c("AsymAD", "AD"),
  c("AsymAD", "Control"),
  c("Control", "AD"))
marker_scores_long$diagnosis <- factor(marker_scores_long$diagnosis, levels = c("Control", "AD", "AsymAD"))
# Single plot with all 5 cell types in one row
plot <- ggboxplot(marker_scores_long, x = "diagnosis", y = "score",
          color = "black", fill = "diagnosis", 
          palette = custom_palette, 
          facet.by = "cell_type", nrow = 1, short.panel.labs = TRUE,
          ) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") +
  labs(title = "Cell Type Marker Scores by Diagnosis",
       y = "Average Marker Expression", x = "Diagnosis") +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_blank(),       # Removes box around facet labels
    strip.text = element_text(face = "bold"), # Keeps facet title bold
    axis.line = element_line(),               # Draw axis lines
    panel.border = element_blank(),           # Remove full panel borders
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black"),
    axis.title = element_text(face = "bold", colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )
plot
ggsave("Marker_scores_boxplot.png", plot = plot, dpi = 600, bg = "transparent", limitsize = FALSE)

#========================================
# Fishers test for before and after comparison for AD vs Control

# Running Fisher's test on significant proteins before cell_type_correction. 
# Set FDR threshold
#sig proteins computed from DE_limma.R
sig_proteins <- rownames(res_table[res_table$adj.P.Val < 0.05, ])
# Background = all proteins tested
background_genes <- rownames(res_table)
# Fisher test results
fisher_results <- marker_long %>%
  filter(gene %in% background_genes) %>%
  group_by(cell_type) %>%
  summarise(
    # Number of marker genes in significant set
    marker_in_sig = sum(gene %in% sig_proteins),
    # Marker genes not in significant set
    marker_not_sig = sum(!(gene %in% sig_proteins)),
    # Non-marker genes in significant set
    non_marker_in_sig = sum(!background_genes %in% gene & background_genes %in% sig_proteins),
    # Non-marker genes not in significant set
    non_marker_not_sig = sum(!background_genes %in% gene & !background_genes %in% sig_proteins),
    
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(matrix(c(marker_in_sig, marker_not_sig, non_marker_in_sig, non_marker_not_sig), 
                                 nrow = 2))$p.value
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)
write.csv(fisher_results, "fishers_results_before_AD_Control.csv")

# Running Fisher's test on significant proteins after cell_type_correction. 
# Set FDR threshold
#sig proteins computed from DE_limma_cell_type_correction.R
sig_proteins <- rownames(res_table[res_table$adj.P.Val < 0.05, ])
# Background = all proteins tested
background_genes <- rownames(res_table)
# Fisher test results
fisher_results <- marker_long %>%
  filter(gene %in% background_genes) %>%
  group_by(cell_type) %>%
  summarise(
    # Number of marker genes in significant set
    marker_in_sig = sum(gene %in% sig_proteins),
    # Marker genes not in significant set
    marker_not_sig = sum(!(gene %in% sig_proteins)),
    # Non-marker genes in significant set
    non_marker_in_sig = sum(!background_genes %in% gene & background_genes %in% sig_proteins),
    # Non-marker genes not in significant set
    non_marker_not_sig = sum(!background_genes %in% gene & !background_genes %in% sig_proteins),
    
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(matrix(c(marker_in_sig, marker_not_sig, non_marker_in_sig, non_marker_not_sig), 
                                 nrow = 2))$p.value
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)
write.csv(fisher_results, "fishers_results_after_AD_Control.csv")

# Reading the files and combining them

fishers_results_before <- read_csv("fishers_results_before_AD_Control.csv")
fishers_results_after <- read_csv("fishers_results_after_AD_Control.csv")
fishers_results_before$Condition <- "Before"
fishers_results_after$Condition <- "After"
combined_fisher <- rbind(fishers_results_before, fishers_results_after)
combined_fisher$Condition <- factor(combined_fisher$Condition, levels = c("After", "Before"))
combined_fisher$cell_type <- factor(combined_fisher$cell_type, levels = c( "Endothelia", "Oligodendrocytes", "Astrocytes", "Microglia", "Neuron" ))

# Plotting the results
plot <- ggplot(combined_fisher, aes(x = cell_type, 
                                    y = -log10(p_adj), 
                                    fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1.0)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Before" = "steelblue", "After" = "tomato")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(title = "Cell Type Enrichment (Fisher Test)",
       x = "Cell Type", y = "-log10(FDR-adjusted p-value)") +
  coord_flip()

plot
ggsave("Fishers_combined_AD_Control.png",
       plot = plot,
       dpi = 600,
       bg = "transparent",
       width = 6.16,
       height = 3.6,
       units = "in",
       limitsize = FALSE)

#========================================
# Fishers test for before and after comparison for AsymAD vs AD

# Running Fisher's test on significant proteins before cell_type_correction. 
# Set FDR threshold
#sig proteins computed from DE_limma.R
sig_proteins <- rownames(res_table[res_table$adj.P.Val < 0.05, ])
# Background = all proteins tested
background_genes <- rownames(res_table)
# Fisher test results
fisher_results <- marker_long %>%
  filter(gene %in% background_genes) %>%
  group_by(cell_type) %>%
  summarise(
    # Number of marker genes in significant set
    marker_in_sig = sum(gene %in% sig_proteins),
    # Marker genes not in significant set
    marker_not_sig = sum(!(gene %in% sig_proteins)),
    # Non-marker genes in significant set
    non_marker_in_sig = sum(!background_genes %in% gene & background_genes %in% sig_proteins),
    # Non-marker genes not in significant set
    non_marker_not_sig = sum(!background_genes %in% gene & !background_genes %in% sig_proteins),
    
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(matrix(c(marker_in_sig, marker_not_sig, non_marker_in_sig, non_marker_not_sig), 
                                 nrow = 2))$p.value
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)
write.csv(fisher_results, "fishers_results_before_AsymAD_AD.csv")

# Running Fisher's test on significant proteins after cell_type_correction. 
# Set FDR threshold
#sig proteins computed from DE_limma_cell_type_correction.R
sig_proteins <- rownames(res_table[res_table$adj.P.Val < 0.05, ])
# Background = all proteins tested
background_genes <- rownames(res_table)
# Fisher test results
fisher_results <- marker_long %>%
  filter(gene %in% background_genes) %>%
  group_by(cell_type) %>%
  summarise(
    # Number of marker genes in significant set
    marker_in_sig = sum(gene %in% sig_proteins),
    # Marker genes not in significant set
    marker_not_sig = sum(!(gene %in% sig_proteins)),
    # Non-marker genes in significant set
    non_marker_in_sig = sum(!background_genes %in% gene & background_genes %in% sig_proteins),
    # Non-marker genes not in significant set
    non_marker_not_sig = sum(!background_genes %in% gene & !background_genes %in% sig_proteins),
    
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(matrix(c(marker_in_sig, marker_not_sig, non_marker_in_sig, non_marker_not_sig), 
                                 nrow = 2))$p.value
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)
write.csv(fisher_results, "fishers_results_after_AsymAD_AD.csv")

# Reading the files and combining them

fishers_results_before <- read_csv("fishers_results_before_AsymAD_AD.csv")
fishers_results_after <- read_csv("fishers_results_after_AsymAD_AD.csv")
fishers_results_before$Condition <- "Before"
fishers_results_after$Condition <- "After"
combined_fisher <- rbind(fishers_results_before, fishers_results_after)
combined_fisher$Condition <- factor(combined_fisher$Condition, levels = c("After", "Before"))
combined_fisher$cell_type <- factor(combined_fisher$cell_type, levels = c( "Endothelia", "Oligodendrocytes", "Astrocytes", "Microglia", "Neuron" ))

# Plotting the results
plot <- ggplot(combined_fisher, aes(x = cell_type, 
                                    y = -log10(p_adj), 
                                    fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1.0)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Before" = "steelblue", "After" = "tomato")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(title = "Cell Type Enrichment (Fisher Test)",
       x = "Cell Type", y = "-log10(FDR-adjusted p-value)") +
  coord_flip()

plot
ggsave("Fishers_combined_AsymAD_AD.png",
       plot = plot,
       dpi = 600,
       bg = "transparent",
       width = 6.16,
       height = 3.6,
       units = "in",
       limitsize = FALSE)

#========================================
# Fishers test for before and after comparison for AsymAD vs Control

# Running Fisher's test on significant proteins before cell_type_correction. 
# Set FDR threshold
#sig proteins computed from DE_limma.R
sig_proteins <- rownames(res_table[res_table$adj.P.Val < 0.05, ])
# Background = all proteins tested
background_genes <- rownames(res_table)
# Fisher test results
fisher_results <- marker_long %>%
  filter(gene %in% background_genes) %>%
  group_by(cell_type) %>%
  summarise(
    # Number of marker genes in significant set
    marker_in_sig = sum(gene %in% sig_proteins),
    # Marker genes not in significant set
    marker_not_sig = sum(!(gene %in% sig_proteins)),
    # Non-marker genes in significant set
    non_marker_in_sig = sum(!background_genes %in% gene & background_genes %in% sig_proteins),
    # Non-marker genes not in significant set
    non_marker_not_sig = sum(!background_genes %in% gene & !background_genes %in% sig_proteins),
    
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(matrix(c(marker_in_sig, marker_not_sig, non_marker_in_sig, non_marker_not_sig), 
                                 nrow = 2))$p.value
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)
write.csv(fisher_results, "fishers_results_before_AsymAD_Control.csv")

# Running Fisher's test on significant proteins after cell_type_correction. 
# Set FDR threshold
#sig proteins computed from DE_limma_cell_type_correction.R
sig_proteins <- rownames(res_table[res_table$adj.P.Val < 0.05, ])
# Background = all proteins tested
background_genes <- rownames(res_table)
# Fisher test results
fisher_results <- marker_long %>%
  filter(gene %in% background_genes) %>%
  group_by(cell_type) %>%
  summarise(
    # Number of marker genes in significant set
    marker_in_sig = sum(gene %in% sig_proteins),
    # Marker genes not in significant set
    marker_not_sig = sum(!(gene %in% sig_proteins)),
    # Non-marker genes in significant set
    non_marker_in_sig = sum(!background_genes %in% gene & background_genes %in% sig_proteins),
    # Non-marker genes not in significant set
    non_marker_not_sig = sum(!background_genes %in% gene & !background_genes %in% sig_proteins),
    
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(matrix(c(marker_in_sig, marker_not_sig, non_marker_in_sig, non_marker_not_sig), 
                                 nrow = 2))$p.value
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)
write.csv(fisher_results, "fishers_results_after_AsymAD_Control.csv")

# Reading the files and combining them

fishers_results_before <- read_csv("fishers_results_before_AsymAD_Control.csv")
fishers_results_after <- read_csv("fishers_results_after_AsymAD_Control.csv")
fishers_results_before$Condition <- "Before"
fishers_results_after$Condition <- "After"
combined_fisher <- rbind(fishers_results_before, fishers_results_after)
combined_fisher$Condition <- factor(combined_fisher$Condition, levels = c("After", "Before"))
combined_fisher$cell_type <- factor(combined_fisher$cell_type, levels = c( "Endothelia", "Oligodendrocytes", "Astrocytes", "Microglia", "Neuron" ))

# Plotting the results
plot <- ggplot(combined_fisher, aes(x = cell_type, 
                                    y = -log10(p_adj), 
                                    fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1.0)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Before" = "steelblue", "After" = "tomato")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(title = "Cell Type Enrichment (Fisher Test)",
       x = "Cell Type", y = "-log10(FDR-adjusted p-value)") +
  coord_flip()

plot
ggsave("Fishers_combined_AsymAD_Control.png",
       plot = plot,
       dpi = 600,
       bg = "transparent",
       width = 6.16,
       height = 3.6,
       units = "in",
       limitsize = FALSE)

#============================================


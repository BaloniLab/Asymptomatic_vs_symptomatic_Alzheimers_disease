# Loading packages

library(ggplot2)
library(ggrepel)

# Cell type groups
exc_cells <- c("Exc superficial neuron", "Exc mid-layer neuron", "Exc deep-layer neuron")
inh_cells <- c("Inh PVALB", "Inh SST", "Inh VIP", "Inh LAMP5")
glia_cells <- c("Ast CHI3L1", "Ast DPP10", "Ast GRM3", "Oligodendrocyte", "Mic P2RY12", "OPC")
# deg_df for Excitatory only
exc_df <- deg_df %>%
  filter(cell_type %in% exc_cells) %>%
  mutate(
    p_val_adj = ifelse(p_val_adj == 0, 1e-300, p_val_adj),
    log10_p = -log10(p_val_adj),
    direction = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )


# Top 20 by abs(log2FC)
top_logfc_genes <- exc_df %>%
  filter(p_val_adj < 0.01) %>%
  slice_max(order_by = abs(avg_log2FC), n = 26, with_ties = FALSE)

# Top 20 by significance (smallest p-value)
top_signif_genes <- exc_df %>%
  slice_min(order_by = p_val_adj, n = 6, with_ties = FALSE)

# Combining the gene–cell_type pairs from both
top_gene_pairs <- bind_rows(top_logfc_genes, top_signif_genes) %>%
  distinct(gene, cell_type)

# Adding label info to the main dataframe
exc_df <- exc_df %>%
  mutate(
    label = ifelse(paste0(gene, cell_type) %in% paste0(top_gene_pairs$gene, top_gene_pairs$cell_type), gene, NA),
    label_color = ifelse(!is.na(label), cell_type, NA)
  )

# Labeling colors for Exc subtypes
label_colors_exc <- c(
    "Exc superficial neuron" = "orange3",
    "Exc mid-layer neuron"   = "#377EB8",
    "Exc deep-layer neuron"  = "#4DAF4A",
    "Inh PVALB"              = "#984EA3",
    "Inh SST"                = "#FF7F00",
    "Inh VIP"                = "#A65628",
    "Inh LAMP5"              = "#F781BF",
    "Ast CHI3L1"             = "#999999",
    "Ast DPP10"              = "#66C2A5",
    "Ast GRM3"               = "#FC8D62",
    "Mic P2RY12"             = "#8DA0CB",
    "Oligodendrocyte"        = "indianred4",
    "OPC"                    = "#A6D854"
)

# Plot
plot_final <- ggplot(exc_df, aes(x = avg_log2FC, y = log10_p)) +
  # Dot color: DEG direction (mapped in geom_point)
  geom_point(aes(color = direction), alpha = 0.6, size = 2.5, show.legend = TRUE) +
  
  # Gene label: subtype colors
  geom_text_repel(
    aes(label = label, color = label_color),
    size = 6, fontface = "bold",
    show.legend = FALSE, na.rm = TRUE, max.overlaps = Inf 
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      # DEG direction
      "Up" = "firebrick3",
      "Down" = "lightblue",
      "NS" = "gray70",
      # Subtype label colors
      label_colors_exc 
    ),
    breaks = c("Up", "Down", "NS", names(label_colors_exc)),
    guide = guide_legend(override.aes = list(size = 4, shape = 16))
  ) +
  
  labs(
    title = "Volcano Plot - Excitatory",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  scale_y_continuous(limits = c(0, 310)) 

plot_final
ggsave("volcano_Exc.png", plot_final, dpi = 600, limitsize = FALSE, bg = "transparent")
#===============================================================

# Barplot for pathway
# Loading all packages
library(dplyr)
library(ggplot2)
library(forcats)
# all_df has: cell_type, Term (pathway), direction ∈ {"Upregulated","Downregulated"}, count
all_df <- rbind(Exc,Inh,Glia)
# 1) Downregulated counts negative
plot_df <- all_df %>%
  mutate(count_signed = ifelse(direction == "Downregulated", -count, count))

# Cell types order
plot_df <- plot_df %>%
  mutate(cell_type = dplyr::recode(as.character(cell_type),
                                   
                                   `Exc deep-layer`   = "Exc deep-layer neuron",
                                   `Exc mid-layer`    = "Exc mid-layer neuron",
                                   `Exc superficial`  = "Exc superficial neuron",
                                   LAMP5 = "Inh LAMP5",
                                   PV = "Inh PVALB",
                                   SST   = "Inh SST",
                                   VIP   = "Inh VIP",
                                   `Microglia (homeostatic)` = "Mic P2RY12",
                                   
                                   Oligodendrocyte = "Oligodendrocyte",
                                   OPC = "OPC",
                                   .default = as.character(cell_type)   
  ))

plot_df$cell_type <- factor(
  plot_df$cell_type,
  levels = c(
    
    "Exc superficial neuron","Exc mid-layer neuron", "Exc deep-layer neuron",
    "Inh PVALB","Inh SST","Inh VIP", "Inh LAMP5",
    "Ast CHI3L1","Ast GRM3",
    "Mic P2RY12",
    "Oligodendrocyte",
    "OPC"
  )
)


# Pathways by overall abundance 
plot_df <- plot_df %>%
  group_by(Term) %>% mutate(total_term = sum(abs(count_signed))) %>% ungroup() %>%
  mutate(Term = fct_reorder(Term, total_term, .desc = TRUE))

# 2) Plot: single panel, up above 0, down below 0, stacked by pathway
plot <- ggplot(plot_df, aes(x = cell_type, y = count_signed, fill = Term)) +
  geom_col(width = 0.8, color = "black") +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Cell type", y = "Number of  DE genes", fill = "Pathway")+
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
        axis.text.y  = element_text(face = "bold", color = "black"),
        axis.title.x  = element_text(face = "bold", color = "black"),
        legend.text = element_text(face = "bold", color = "black"),
        legend.title = element_text(face = "bold", color = "black")
  )
plot
ggsave("Pathway_plot_all_cell_types.png", plot, dpi = 600, limitsize = FALSE, bg = "transparent")



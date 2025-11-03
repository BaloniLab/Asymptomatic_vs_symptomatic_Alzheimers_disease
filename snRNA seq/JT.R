# Loading all packages
library(dplyr)
library(ggplot2)
library(clinfun)
library(readr)
# Reading the metadata which contains all the cell types annotated.
metadata_final_cell_types_14_types <- read_csv("/depot/pbaloni/data/Lab_members/Purba_Mandal/snRNA_seq_analysis/metadata_final_cell_types_14_types.csv")
meta <- metadata_final_cell_types_14_types
meta$diagnosis <- factor(meta$diagnosis, levels = c("Control", "AsymAD", "AD"))
# Grouping by projid (donor ID) to count how many cells belong to each donor (n()) and then store that as total_cells.
# Per-donor total counts
total_counts <- meta %>%
  group_by(projid) %>%
  summarise(total_cells = n(), .groups = "drop")

# Grouping by donor × diagnosis × cell type, count cells of each type per donor (n_cells) and 
# joining to total_counts so we know how many total cells that donor had.
# Calculation of prop = n_cells / total_cells = proportion of a given cell type within each donor.
# Per-donor per-cell type counts
ct_counts <- meta %>%
  group_by(projid, diagnosis, cell_type) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  left_join(total_counts, by = "projid") %>%
  mutate(prop = n_cells / total_cells)

# Function to make one plot per cell type
# Filters ct_counts for one cell type.
# Runs a Jonckheere–Terpstra test (non-parametric trend test across ordered groups) using proportions vs. numeric-coded diagnosis.
# Wrapped in tryCatch to avoid crashing if test fails.
# Creates a boxplot + jitter of donor proportions per diagnosis.
# Annotates the plot with the JT p-value in italics under the title.
plot_ct <- function(ct_name) {
  df <- ct_counts %>% filter(cell_type == ct_name)
  
  # Trend test p-value
  set.seed(123)
  jt <- tryCatch(
    jonckheere.test(df$prop, as.numeric(df$diagnosis), nperm = 10000)$p.value,
    error = function(e) NA_real_
  )
  p_lab <- ifelse(is.na(jt), "p = NA",
                  paste0("JT p = ", signif(jt, 3)))
  
  ggplot(df, aes(x = diagnosis, y = prop, fill = diagnosis)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.6, color = "black") +
    geom_jitter(width = 0.15, size = 1.8, alpha = 0.85) +
    scale_fill_manual(values = c("#99ccff", "plum3", "brown")) +
    stat_compare_means(comparisons = list(
      c("Control", "AsymAD"),
      c("Control", "AD"),
      c("AsymAD", "AD")
    ), 
    method = "wilcox", 
    label = "p.signif",
    size = 5) +
    labs(
      title = ct_name,
      subtitle = p_lab,
      x = NULL, y = "Proportion"
    ) +
    theme_classic(base_size = 18) +
    theme(legend.position = "none",
          plot.subtitle = element_text(face = "italic"),
          axis.text.y = element_text(face = "bold",  color = "black"),
          axis.text.x = element_text(face = "bold", color = "black", angle = 45, vjust = 1, hjust = 1),
          axis.title.y  = element_blank()
          )
}

# List of cell types
cell_types <- unique(ct_counts$cell_type)

# Plots for each cell type
plots <- lapply(cell_types, plot_ct)

# Function to split a list into chunks of given size
split_into_chunks <- function(lst, chunk_size) {
  split(lst, ceiling(seq_along(lst) / chunk_size))
}

# Splitting into groups of 5 plots
plot_chunks <- split_into_chunks(plots, 5)

# Combining each chunk into its own figure
figures <- lapply(plot_chunks, function(chunk) {
  wrap_plots(chunk, ncol = 5) 
})

# Printing each figure
plot <- figures[[1]]
plot
ggsave("JT_1.png", plot, dpi = 600, limitsize = FALSE, bg = "transparent")

# Combining into one figure (grid layout)
combined_plot <- wrap_plots(plots, ncol = 7, nrow = 2) + 
  plot_annotation(title = "Per-donor Proportion by Diagnosis",
                  subtitle = "Separate plots for each cell type with Jonckheere–Terpstra p-values")

combined_plot

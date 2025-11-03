# Loading the packages

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)

# Reading CSV
df <- read_csv("OXPHOS_complexes(All_complexes).csv")

# Removing any row where any column contains MT-
df_no_mt <- df %>%
  filter(!if_any(everything(), ~ grepl("^MT-", .)))

# Saving cleaned file
write.csv(df_no_mt, "nuclear_encoded_complex_genes_noMT.csv", row.names = FALSE)

# Deg file
All_deg_cell_types_AsymAD_vs_AD <- All_deg_cell_types_AsymAD_vs_AD %>% 
  filter(p_val_adj < 0.01) %>% 
  filter(!grepl("MT-", gene))
#=====
# Converting to long format: two columns gene + complex
complex_map <- df_no_mt %>%
  pivot_longer(
    cols = everything(),               # all columns are complexes
    names_to = "complex",               # new column for complex names
    values_to = "gene"                   # new column for gene names
  ) %>%
  filter(!is.na(gene) & gene != "")      # remove NAs or empty cells

# Saving as mapping file
write_csv(complex_map, "complex_mapping.csv")
#======

# Loading DEG file
deg <- All_deg_cell_types_AsymAD_vs_AD
# Mapping the DEGs to complexes and tagging the direction

deg_dir <- deg %>%
  mutate(direction = ifelse(avg_log2FC > 0, "Up",
                            ifelse(avg_log2FC < 0, "Down", NA))) %>%
  filter(!is.na(direction)) %>%                
  inner_join(complex_map, by = "gene")          # add complex

# Complex totals (denominator for % coverage)
complex_totals <- complex_map %>%
  distinct(complex, gene) %>%
  count(complex, name = "total_genes")

# Coverage by direction (so Up% + Down% = total coverage %)
coverage_dir <- deg_dir %>%
  distinct(cell_type, complex, direction, gene) %>%
  count(cell_type, complex, direction, name = "n_dir") %>%
  left_join(complex_totals, by = "complex") %>%
  mutate(coverage_pct = 100 * n_dir / total_genes)

# ordering factors
coverage_dir$complex <- factor(
  coverage_dir$complex,
  levels = c("Complex I","Complex II","Complex III","Complex IV","Complex V")
)

coverage_dir <- coverage_dir %>%
  mutate(cell_type = dplyr::recode(as.character(cell_type),
                                   Ast_CHI3L1 = "Ast CHI3L1",
                      
                                   Ast_GRM3   = "Ast GRM3",
                                   `Deep-layer neuron`   = "Exc deep-layer neuron",
                                   `Mid-layer neuron`    = "Exc mid-layer neuron",
                                   `Superficial neuron`  = "Exc superficial neuron",
                                   Inh_LAMP5 = "Inh LAMP5",
                                   Inh_PVALB = "Inh PVALB",
                                   Inh_SST   = "Inh SST",
                                   Inh_VIP   = "Inh VIP",
                                  
                                   Oligodendrocyte = "Oligodendrocyte",
                                   OPC = "OPC",
                                   .default = as.character(cell_type)   
  ))

# Setting order
coverage_dir$cell_type <- factor(
  coverage_dir$cell_type,
  levels = c(
    
    "Exc superficial neuron","Exc mid-layer neuron", "Exc deep-layer neuron",
    "Inh PVALB","Inh SST","Inh VIP", "Inh LAMP5",
    "Ast CHI3L1","Ast GRM3",
    "Oligodendrocyte",
    "OPC"
  )
)


# Plot: stacked by direction and the bar height is total complex coverage
p_facet <- ggplot(
  coverage_dir,
  aes(x = fct_rev(cell_type), y = coverage_pct, fill = direction)
) +
  geom_col(width = 0.8, color = "black") +
  geom_text(
    data = subset(coverage_dir, coverage_pct > 5),
    aes(label = paste0(round(coverage_pct, 1))),
    position = position_stack(vjust = 0),  # end of each stacked piece
    hjust = 0,                         # push a bit past the bar edge
    size = 5, fontface = "bold"
  )+
  # label segments only if they are large enough
  facet_wrap(~ complex, nrow = 1) +
  coord_flip(clip = "off") +
  scale_y_continuous(labels = percent_format(scale = 1),
                     limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.06))) +
  scale_fill_manual(values = c(Up = "#ff4d4d", Down = "#7570B3")) +  # red/orange & blue
  labs(x = NULL, y = "Coverage of complex (%)",
       title = "Coverage of ETC Complexes by Cell Type (Up vs Down)") +
  theme_classic(base_size = 16) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    axis.text.y  = element_text(face = "bold", color = "black"),
    axis.title.x  = element_text(face = "bold", color = "black"),
    axis.text.x  = element_text(face = "bold", colour = "black", size = 9),
    legend.title = element_blank(),
    panel.spacing.x = unit(10, "pt"),
    plot.margin = margin(10, 20, 10, 10)
  )

p_facet
ggsave("complex_facet.png", p_facet, dpi = 600, limitsize = FALSE, bg = "transparent")

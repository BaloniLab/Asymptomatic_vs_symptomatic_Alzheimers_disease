# Load all packages
library(dplyr)
library(ggplot2)
library(forcats)
library(scales)

# Setting factor based on cell type subtype
prop_df$cell_type <- factor(prop_df$cell_type, levels = c("Exc", "Inh", "Astrocytes", "Microglia", "Oli", "OPC"))
prop_df <- metadata_final_cell_types %>%
  count(cell_type, name = "n") %>%
  mutate(prop = n / sum(n))

# Setting colors for each cell type
celltype_colors <- c(
  "Exc"       = "#ff4d4d",
  "Inh" = "#004c99",
  "Astrocytes"   = "#cc6600",
  "Microglia" = "darkgreen",
  "Oli" = "plum3",
  "OPC" = "pink3"
)
# Plotting the data
p <- ggplot(prop_df, aes(x = "All cells", y = prop, fill = cell_type)) +
  geom_col(width = 0.6, color = "black") +
  geom_text(
    aes(label = scales::percent(prop, accuracy = 0.1)),
    position = position_stack(vjust = 0.5),   # center of each segment
    size = 10, fontface = "bold", color = "black"
  ) +
  scale_y_continuous(labels = percent, limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = celltype_colors) + 
  labs(x = NULL, y = "Proportion", title = "Cell class composition") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.text  = element_text(face = "bold", size = 16, color = "black")
  ) 
p
# Saving the data
ggsave("cell_proportion.png", p, dpi = 600, limitsize = FALSE, bg = "transparent")

#=================

# Plotting now based on diagnosis
prop_df <- all_md
all_md$diagnosis <- factor(all_md$diagnosis, levels = c("AD", "AsymAD", "Control"))

# Percentage covered by each cell type
plot_df <- all_md %>%
  count(cell_type, cell_type_high_resolution, name = "n_cells") %>%
  group_by(cell_type) %>%
  mutate(pct = 100 * n_cells / sum(n_cells)) %>%
  ungroup()

# Setting factor order
plot_df$cell_type <- factor(plot_df$cell_type, 
                            levels = c("Exc", "Inh", "Astrocyte", "Microglia", "Oli", "OPC" ))

# Reversing order in the plot (top to bottom):
plot_df$cell_type <- forcats::fct_rev(plot_df$cell_type)


# Defining the colo
celltype_colors <- c(
  # Excitatory neurons - red shades
  "Exc deep-layer neuron"       = "#ff9999",
  "Exc mid-layer neuron"        = "#ff4d4d",
  "Exc superficial neuron"      = "#b30000",
  
  # Inhibitory neurons - blue shades
  "Inh LAMP5" = "#99ccff",
  "Inh PVALB" = "#004c99",
  "Inh SST"   = "#3399ff",
  "Inh VIP"   = "#0066cc",
  
  # Astrocytes - orange shades
  "Ast CHI3L1" = "#ffb366",
  "Ast DPP10"  = "#ff9933",
  "Ast GRM3"   = "#cc6600",
  # Microglia - green shades
  "Mic P2RY12" = "darkgreen",
  "Mic TPT1"   = "seagreen3",
  
  # Oligodendrocytes - purple shades
  "Oli"                     = "plum3",
  
  # OPC - pink
  "OPC" = "pink3"
)
# Factoring each cell type
plot_df$cell_type_high_resolution <- factor(
  plot_df$cell_type_high_resolution,
  levels = c(
    "Exc superficial neuron", "Exc mid-layer neuron", "Exc deep-layer neuron",
    "Inh PVALB","Inh SST","Inh VIP","Inh LAMP5",
    "Ast CHI3L1","Ast DPP10","Ast GRM3",
    "Mic P2RY12","Mic TPT1",
    "Oli",
    "OPC"
  )
)
# Plotting the data
plot <- ggplot(plot_df, aes(
  y = (cell_type), 
  x = pct, 
  fill = cell_type_high_resolution
)) +
  geom_col(width = 1, color = "black") +
  scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 101), expand = c(0,0)) +
  labs(x = "Percent of cells", y = "Cell type", fill = "cell_type_high_resolution"
  ) +
  scale_fill_manual(values = celltype_colors) +  # <-- set custom colors
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "bold", size = 18, color = "black"),
    axis.text.x = element_text(face = "bold", size = 18, color = "black"),
    axis.title.x  = element_text(face = "bold", size = 18, color = "black"),
    axis.title.y  = element_blank(),
    legend.title = element_blank(),
    legend.text  = element_text(face = "bold", size = 16, color = "black")
  )

plot
ggsave("cell_type.png", plot, dpi = 600, limitsize = FALSE, bg = "transparent")

#==================================
# Plotting data based on sex
plot_df <- all_md %>%
  count(cell_type, sex, name = "n_cells") %>%
  group_by(cell_type) %>%
  mutate(pct = 100 * n_cells / sum(n_cells)) %>%
  ungroup()

plot <- ggplot(plot_df, aes(y = fct_rev(cell_type), x = pct, fill = sex)) +
  geom_col(width = 1, color = "black") +
  scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100), expand = c(0,0)) +
  labs(x = "Percent of cells", y = "Cell type", fill = "Sex") +
  theme_classic(base_size = 14) +
  theme(axis.text.y = element_text(face = "bold", size = 18, color = "black"),
        axis.text.x = element_text(face = "bold", size = 18, color = "black"),
        axis.title.x  = element_text(face = "bold", size = 18, color = "black"),
        axis.title.y  = element_text(face = "bold", size = 18, color = "black"),
        legend.title = element_text(face = "bold"),
        legend.text  = element_text(face = "bold", size = 16, color = "black"))
plot
ggsave("Sex.png", plot, dpi = 600, limitsize = FALSE, bg = "transparent")

#=======================

# Plotting data based on diagnosis
plot_df <- all_md %>%
  count(cell_type, diagnosis, name = "n_cells") %>%
  group_by(cell_type) %>%
  mutate(pct = 100 * n_cells / sum(n_cells)) %>%
  ungroup()
# Diagnosis color
diagnosis_colors <- c("Control" = "#99ccff",
                      "AD" = "brown",
                      "AsymAD" = "plum3")
# Plotting the data
plot <- ggplot(plot_df, aes(y = fct_rev(cell_type), x = pct, fill = diagnosis)) +
  geom_col(width = 1, color = "black") +
  scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100), expand = c(0,0)) +
  scale_fill_manual(values = diagnosis_colors) + 
  labs(x = "Percent of cells", y = "Cell type", fill = "diagnosis") +
  theme_classic(base_size = 14) +
  theme(axis.text.y = element_text(face = "bold", size = 18, color = "black"),
        axis.text.x = element_text(face = "bold", size = 18, color = "black"),
        axis.title.x  = element_text(face = "bold", size = 18, color = "black"),
        axis.title.y  = element_text(face = "bold", size = 18, color = "black"),
        legend.title = element_text(face = "bold"),
        legend.text  = element_text(face = "bold", size = 16, color = "black"))
plot
ggsave("diagnosis.png", plot, dpi = 600, limitsize = FALSE, bg = "transparent")


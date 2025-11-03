# DE Limma with cell_type_correction

# To evaluate DEG without have bias of cell type proportions...

# This code performs differential expression analysis (DEA) on proteomics data while adjusting for cell type proportions.
# Step 0: Setup
# Step 1: Load cell type marker genes
# Step 2: Load expression matrix
# Step 3–4: Aggregate by gene name
# Step 5: Score cell types per sample
# Step 6: Add diagnosis metadata
# Step 7: Differential expression with cell type adjustment

#==================================
# Load packages
library(tidyr)
library(dplyr)
library(limma)
# Step 1: Load the CSV, skipping the first row
setwd("/depot/pbaloni/data/Lab_members/Purba_Mandal/Proteomics")
marker_df <- read.csv("Cell_type_markers_5_types.csv")

# Reshape marker_df into long format: gene | cell_type
marker_long <- pivot_longer(marker_df, cols = everything(),
                            names_to = "cell_type", values_to = "gene") %>%
  filter(!is.na(gene)) %>% distinct()
#===============================

# Step 2: Load the Exp_matrix 

# BA9
setwd("/depot/pbaloni/data/Lab_members/Purba_Mandal/Proteomics/Brain_regions/BA9(Banner+ROSMAP)")
exp_mat_sub <- read.csv("3.cleanDat-Consensus-MainArm.csv",row.names = 1) 
rownames(exp_mat_sub)

#Extract protein names from rownames
gene_names <- sapply(strsplit(rownames(exp_mat_sub), "\\|"), function(x) x[1])

# Convert expression matrix to data frame and add gene as a proper column
exp_df <- as.data.frame(exp_mat_sub)
exp_df$gene <- gene_names  

#======================
# Step 3: Group by protein isoforms, aggregate with mean

exp_df_agg <- exp_df %>%
  group_by(gene) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

# Step 4: Convert back to matrix, set rownames
exp_mat_sub_unique <- as.matrix(exp_df_agg[,-1])  # remove gene column
rownames(exp_mat_sub_unique) <- exp_df_agg$gene
write.csv(exp_mat_sub_unique, "Aggregated_expression_matrix.csv", row.names = TRUE)
#=======================

# Step 5: For each cell type, calculate average expression of its marker genes per sample
# defines a custom function named get_marker_score, 
# which calculates the average expression of marker genes for a given cell type across all samples.
get_marker_score <- function(cell_type, mat, marker_table) {
  genes <- marker_table %>% filter(cell_type == !!cell_type) %>% pull(gene)
  genes <- intersect(genes, rownames(mat))  # only keep genes present in the dataset
  if (length(genes) < 3) return(rep(NA, ncol(mat)))  # skip if too few markers
  colMeans(mat[genes, , drop = FALSE], na.rm = TRUE)
}

# Step 6: Build matrix of scores: samples x cell types
cell_types <- unique(marker_long$cell_type) # This creates a list of all unique cell types in marker table
marker_scores <- sapply(cell_types, get_marker_score, mat = exp_mat_sub_unique, marker_table = marker_long)
# Applies the get_marker_score() function to each cell type.
# Returns a numeric vector of marker scores for that cell type across all samples.

# Convert to data.frame and add diagnosis info
marker_scores_df <- as.data.frame(marker_scores)
marker_scores_df$id <- colnames(exp_mat_sub_unique)

# Read metadata
diagnosis_sub <- read.csv("Final_metadata.csv")

# Reorder diagnosis to match expression matrix columns

marker_scores_df$id == diagnosis_sub$sample_id
marker_scores_df$diagnosis <- diagnosis_sub$diagnosis

#==========
# Step 7: Calculate DEA

# Only keep AsymAD and AD
marker_scores <- marker_scores_df
marker_scores_df <- marker_scores_df %>%
  filter(diagnosis %in% c("AsymAD", "AD"))

# Diagnosis is a factor with correct levels
marker_scores_df$diagnosis <- factor(marker_scores_df$diagnosis, levels = c("AD", "AsymAD"))
# Subset expression matrix to match sample IDs
exp_mat_subset_only_2 <- exp_mat_sub_unique[, marker_scores_df$id]
write.csv(exp_mat_subset_only_2, "Aggregated_expression_subset_for_AsymAD_AD.csv", row.names = TRUE)
# Design matrix


design <- model.matrix(~ diagnosis + Astrocytes + Microglia + Neuron + Oligodendrocytes + Endothelia , data = marker_scores_df)


# Fit model
fit <- lmFit(exp_mat_subset_only_2, design)
fit <- eBayes(fit)

# Get DE results for AD vs AsymAD
res_table <- topTable(fit, coef = "diagnosisAsymAD", adjust.method = "BH", number = Inf)

# After fitting your linear model with lmFit(), eBayes():
# Stabilizes the variance estimates (especially helpful with small sample sizes).
# Increases statistical power by shrinking noisy variance estimates toward a common value.
# Helps avoid false positives and false negatives in differential expression analysis.
#=============

# Step 8: Now significant mito_proteins
filtered_data <- res_table %>% 
  filter(adj.P.Val < 0.05) 
filtered_data$protein <- rownames(filtered_data)
mito <- readLines("/depot/pbaloni/data/Lab_members/Purba_Mandal/PFC_MIT_data/DEG_results/mito.txt")
filtered_data_mito <-filtered_data %>% 
  filter( protein%in% mito)
write.csv(filtered_data, "DE_limma_proteins_AsymAD_vs_AD_cell_correction.csv")
write.csv(filtered_data_mito, "DE_limma_mito_proteins_AsymAD_vs_AD_cell_correction.csv")
write.table(filtered_data_mito$protein, file = "DE_limma_mito_proteins_AsymAD_vs_AD_cell_correction.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)


# Ad vs Control
#==========
# Step 7: Calculate DEA

# Only keep AD and Control
marker_scores_df <- marker_scores
marker_scores_df <- marker_scores_df %>%
  filter(diagnosis %in% c("Control", "AD"))

# Diagnosis is a factor with correct levels
marker_scores_df$diagnosis <- factor(marker_scores_df$diagnosis, levels = c("Control", "AD"))
# Subset expression matrix to match sample IDs
exp_mat_subset_only_2 <- exp_mat_sub_unique[, marker_scores_df$id]
write.csv(exp_mat_subset_only_2, "Aggregated_expression_subset_for_AD_vs_Control.csv", row.names = TRUE)
# Design matrix


design <- model.matrix(~ diagnosis + Astrocytes + Microglia + Neuron + Oligodendrocytes + Endothelia , data = marker_scores_df)


# Fit model
fit <- lmFit(exp_mat_subset_only_2, design)
fit <- eBayes(fit)

# Get DE results for AD vs AsymAD
res_table <- topTable(fit, coef = "diagnosisAD", adjust.method = "BH", number = Inf)


#=============

# Step 8: Now significant mito_proteins
filtered_data <- res_table %>% 
  filter(adj.P.Val < 0.05) 
filtered_data$protein <- rownames(filtered_data)
mito <- readLines("/depot/pbaloni/data/Lab_members/Purba_Mandal/PFC_MIT_data/DEG_results/mito.txt")
filtered_data_mito <-filtered_data %>% 
  filter( protein%in% mito)
write.csv(filtered_data, "DE_limma_proteins_AD_vs_Control_cell_correction.csv")
write.csv(filtered_data_mito, "DE_limma_mito_proteins_AD_vs_Control_cell_correction.csv")
write.table(filtered_data_mito$protein, file = "DE_limma_mito_proteins_AD_vs_Control_cell_correction.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

#========================
# AsymAD vs Control
#==========
# Step 7: Calculate DEA

# Only keep AD and Control
marker_scores_df <- marker_scores
marker_scores_df <- marker_scores_df %>%
  filter(diagnosis %in% c("Control", "AsymAD"))

# Diagnosis is a factor with correct levels
marker_scores_df$diagnosis <- factor(marker_scores_df$diagnosis, levels = c("Control", "AsymAD"))
# Subset expression matrix to match sample IDs
exp_mat_subset_only_2 <- exp_mat_sub_unique[, marker_scores_df$id]
write.csv(exp_mat_subset_only_2, "Aggregated_expression_subset_for_AsymAD_vs_Control.csv", row.names = TRUE)
# Design matrix


design <- model.matrix(~ diagnosis + Astrocytes + Microglia + Neuron + Oligodendrocytes + Endothelia , data = marker_scores_df)


# Fit model
fit <- lmFit(exp_mat_subset_only_2, design)
fit <- eBayes(fit)

# Get DE results for AD vs AsymAD
res_table <- topTable(fit, coef = "diagnosisAsymAD", adjust.method = "BH", number = Inf)


#=============

# Step 8: Now significant mito_proteins
filtered_data <- res_table %>% 
  filter(adj.P.Val < 0.05) 
filtered_data$protein <- rownames(filtered_data)
mito <- readLines("/depot/pbaloni/data/Lab_members/Purba_Mandal/PFC_MIT_data/DEG_results/mito.txt")
filtered_data_mito <-filtered_data %>% 
  filter( protein%in% mito)
write.csv(filtered_data, "DE_limma_proteins_AsymAD_vs_Control_cell_correction.csv")
write.csv(filtered_data_mito, "DE_limma_mito_proteins_AsymAD_vs_Control_cell_correction.csv")
write.table(filtered_data_mito$protein, file = "DE_limma_mito_proteins_AsymAD_vs_Control_cell_correction.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)





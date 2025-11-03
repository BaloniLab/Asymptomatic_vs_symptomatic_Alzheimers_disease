# Code generated for applying DE_limma 
# Region BA9
#==================================
# Load packages
library(tidyr)
library(dplyr)
library(limma)

# Read the expr data for BA9
setwd("/depot/pbaloni/data/Lab_members/Purba_Mandal/Proteomics/Brain_regions/BA9(Banner+ROSMAP)")
exp_mat_sub <- read.csv("3.cleanDat-Consensus-MainArm.csv",row.names = 1) 
# Read the metadata file
metadata <- read.csv("Final_metadata.csv") # 488 samples
#===========
#Extract protein names from rownames
gene_names <- sapply(strsplit(rownames(exp_mat_sub), "\\|"), function(x) x[1])
# Convert expression matrix to data frame and add gene as a proper column
exp_df <- as.data.frame(exp_mat_sub)
exp_df$gene <- gene_names  
# Group by protein isoforms, aggregate with mean
exp_df_agg <- exp_df %>%
  group_by(gene) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
# Convert back to matrix, set rownames
exp_mat_sub_unique <- as.matrix(exp_df_agg[,-1])  # remove gene column
rownames(exp_mat_sub_unique) <- exp_df_agg$gene
write.csv(exp_mat_sub_unique, "Aggregated_expression_matrix.csv", row.names = TRUE) # Final aggregated expression data
#=======================

setwd("/depot/pbaloni/data/Lab_members/Purba_Mandal/Proteomics/Brain_regions/BA9(Banner+ROSMAP)")
expr_data <- read.csv("Aggregated_expression_matrix.csv", row.names = 1)

# AD vs AsymAD
metadata <- read.csv("Final_metadata.csv") # 488 samples
metadata <- metadata %>% 
  filter(diagnosis %in% c("AD", "AsymAD"))
expr_data <- expr_data[, colnames(expr_data) %in% metadata$sample_id]

# Match column order of expr_data to metadata
expr_data <- expr_data[, metadata$sample_id]
# Diagnosis must be a factor with AD as reference
metadata$diagnosis <- factor(metadata$diagnosis, levels = c("AD", "AsymAD"))

# DEA
design <- model.matrix(~ diagnosis, data = metadata)
fit <- lmFit(expr_data, design)
fit <- eBayes(fit)

# Differential expression results
results <- topTable(fit, coef = "diagnosisAsymAD", number = Inf)
# Results includes logFC, P.Value, adj.P.Val, t-statistics, etc.
head(results)
# Filtering significant proteins (FDR < 0.05)
sig_proteins <- results %>% filter(adj.P.Val < 0.05)
sig_proteins$proteins <- rownames(sig_proteins)
mito <- readLines("/depot/pbaloni/data/Lab_members/Purba_Mandal/PFC_MIT_data/DEG_results/mito.txt")
sig_proteins <- sig_proteins %>% 
  filter(proteins %in% mito)
# Write to CSV
write.csv(sig_proteins, "DE_limma_proteins_AD_vs_AsymAD_no_cell_correction.csv")
write.table(sig_proteins$proteins, file = "DE_limma_proteins_AD_vs_AsymAD_no_cell_correction.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
#==================

# AD vs Control
metadata <- read.csv("Final_metadata.csv") # 488 samples
metadata <- metadata %>% 
  filter(diagnosis %in% c("AD", "Control"))
expr_data <- read.csv("Aggregated_expression_matrix.csv", row.names = 1)
expr_data <- expr_data[, colnames(expr_data) %in% metadata$sample_id]

# Match column order of expr_data to metadata
expr_data <- expr_data[, metadata$sample_id]
# Diagnosis must be a factor with Control as reference
metadata$diagnosis <- factor(metadata$diagnosis, levels = c("Control", "AD"))

# DEA
design <- model.matrix(~ diagnosis, data = metadata)
fit <- lmFit(expr_data, design)
fit <- eBayes(fit)

# Differential expression results
results <- topTable(fit, coef = "diagnosisAD", number = Inf)
# Results includes logFC, P.Value, adj.P.Val, t-statistics, etc.
head(results)
# Filtering significant proteins (FDR < 0.05)
sig_proteins <- results %>% filter(adj.P.Val < 0.05)
sig_proteins$proteins <- rownames(sig_proteins)
mito <- readLines("/depot/pbaloni/data/Lab_members/Purba_Mandal/PFC_MIT_data/DEG_results/mito.txt")
sig_proteins <- sig_proteins %>% 
  filter(proteins %in% mito)
# Write to CSV
write.csv(sig_proteins, "DE_limma_proteins_AD_vs_Control_no_cell_correction.csv")
write.table(sig_proteins$proteins, file = "DE_limma_proteins_AD_vs_Control_no_cell_correction.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
#==================

# AsymAD vs Control
metadata <- read.csv("Final_metadata.csv") # 488 samples
metadata <- metadata %>% 
  filter(diagnosis %in% c("Control", "AsymAD"))
expr_data <- read.csv("Aggregated_expression_matrix.csv", row.names = 1)
expr_data <- expr_data[, colnames(expr_data) %in% metadata$sample_id]

# Match column order of expr_data to metadata
expr_data <- expr_data[, metadata$sample_id]
# Diagnosis must be a factor with Control as reference
metadata$diagnosis <- factor(metadata$diagnosis, levels = c("Control", "AsymAD"))

# DEA
design <- model.matrix(~ diagnosis, data = metadata)
fit <- lmFit(expr_data, design)
fit <- eBayes(fit)

# Differential expression results
results <- topTable(fit, coef = "diagnosisAsymAD", number = Inf)
# Results includes logFC, P.Value, adj.P.Val, t-statistics, etc.
head(results)
# Filtering significant proteins (FDR < 0.05)
sig_proteins <- results %>% filter(adj.P.Val < 0.05)
sig_proteins$proteins <- rownames(sig_proteins)
mito <- readLines("/depot/pbaloni/data/Lab_members/Purba_Mandal/PFC_MIT_data/DEG_results/mito.txt")
sig_proteins <- sig_proteins %>% 
  filter(proteins %in% mito)
# Write to CSV
write.csv(sig_proteins, "DE_limma_proteins_AsymAD_vs_Control_no_cell_correction.csv")
write.table(sig_proteins$proteins, file = "DE_limma_proteins_AsymADAD_vs_Control_no_cell_correction.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
#==================



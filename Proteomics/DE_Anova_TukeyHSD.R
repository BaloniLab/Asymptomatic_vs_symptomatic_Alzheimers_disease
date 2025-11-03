# Code generated for applying anova followed by Tukeys HSD 
# Region BA9
# BA9

# Read the expression data taken from synapse id(syn25006659) (ROSMAP+Banner cohort)
setwd("/depot/pbaloni/data/Lab_members/Purba_Mandal/Proteomics/Brain_regions/BA9(Banner+ROSMAP)")
expr_data <- read.csv("3.cleanDat-Consensus-MainArm.csv", row.names = 1) #488 samples
# Read the metadata file
metadata <- read.csv("Final_metadata.csv") # 488 samples
diagnosis_sub <- metadata
# Reorder diagnosis to match expression matrix columns
diagnosis_sub <- diagnosis_sub$diagnosis[match(colnames(expr_data), diagnosis_sub$sample_id)]

#===================
# Testing normality
data <- expr_data
data <- na.omit(data)
# Running normality test for a single protein (example: first row)
data_numeric <- as.matrix(data)  
shapiro.test(data_numeric[1, ])

# Applying Shapiro-Wilk test to all proteins
normality_p_values <- apply(data_numeric, 1, function(protein) shapiro.test(protein)$p.value)

# Number of proteins that pass the normality test (p > 0.05 means normal)
sum(normality_p_values > 0.05) / length(normality_p_values) * 100  # Percentage of normally distributed proteins
hist(data_numeric[488, ], main = "Histogram of Protein 1", xlab = "Expression Level", col = "lightblue", border = "black")

#=====================================

# Now we perform pairwise Tukey HSD tests for each protein
# Running Tukey HSD and ensuring correct indexing for p-values
tukey_results <- apply(expr_data, 1, function(x) {
  aov_res <- aov(x ~ diagnosis_sub)
  tukey_res <- TukeyHSD(aov_res)
  return(tukey_res$diagnosis_sub[, "p adj"])  # Get adjusted p-values
})

tukey_results <-data.frame(t(tukey_results)) 
# Spliting row names into 'protein'
tukey_results$protein <- sapply(strsplit(rownames(tukey_results), "\\|"), function(x) x[1])
write.csv(tukey_results, "Tukey_results_proteins.csv", row.names = TRUE)

#=======================

# Now adding logFC for each pair "Control",  "AD", "AsymAD" 

volcano_data <- tukey_results

volcano_data$logFC_AsymAD_vs_Control <- apply(expr_data, 1, function(x) {
  mean(x[diagnosis_sub == "AsymAD"], na.rm = TRUE) - mean(x[diagnosis_sub == "Control"], na.rm = TRUE)
})

volcano_data$logFC_AD_vs_Control <- apply(expr_data, 1, function(x) {
  mean(x[diagnosis_sub == "AD"], na.rm = TRUE) - mean(x[diagnosis_sub == "Control"], na.rm = TRUE)
})

volcano_data$logFC_AsymAD_vs_AD <- apply(expr_data, 1, function(x) {
  mean(x[diagnosis_sub == "AsymAD"], na.rm = TRUE) - mean(x[diagnosis_sub == "AD"], na.rm = TRUE)
})

write.csv(volcano_data, "Tukey_results_proteins_for_volcano_plot.csv", row.names = TRUE)

#=========================================

# All mitochondrial proteins downloaded from Mitocarta3.0
mito <- readLines("/depot/pbaloni/data/Lab_members/Purba_Mandal/PFC_MIT_data/DEG_results/mito.txt")

# AD vs Control
AD_Control <- volcano_data %>% 
  select(Control.AD, logFC_AD_vs_Control, protein) %>% 
  filter(Control.AD < 0.05)
AD_Control_filtered <- na.omit(AD_Control)
Protein <- AD_Control_filtered$protein
Protein <- intersect(mito, Protein)
write.table(Protein, file = "All_mito_proteins_AD_vs_Control.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.csv(AD_Control_filtered, "AD_Control_filtered_all_proteins.csv", row.names = TRUE)

# AsymAD vs Control
AsymAD_Control <- volcano_data %>% 
  select(Control.AsymAD, logFC_AsymAD_vs_Control, protein) %>% 
  filter(Control.AsymAD < 0.05)
AsymAD_Control_filtered <- na.omit(AsymAD_Control)
Protein <- AsymAD_Control_filtered$protein
Protein <- intersect(mito, Protein)
write.table(Protein, file = "All_mito_proteins_AsymAD_Control.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.csv(AsymAD_Control_filtered, "AsymAD_Control_filtered_all_proteins.csv", row.names = TRUE)

# AsymAD vs AD
AsymAD_AD <- volcano_data %>% 
  select(AsymAD.AD, logFC_AsymAD_vs_AD, protein) %>% 
  filter(AsymAD.AD < 0.05)
AsymAD_AD_filtered <- na.omit(AsymAD_AD)
Protein <- AsymAD_AD_filtered$protein
Protein <- intersect(mito, Protein)
write.table(Protein, file = "All_mito_proteins_AsymAD_AD.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.csv(AsymAD_AD_filtered, "AsymAD_AD_filtered_all_proteins.csv", row.names = TRUE)
#=======================================================================








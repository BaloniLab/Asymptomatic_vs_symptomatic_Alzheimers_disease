# Making a heatmap by taking cell_types into account
# exp_mat_sub_unique and marker_scores_df obtained from DE_limma_cell_type_correction.R
# Initialize an empty matrix to store residuals
adjusted_mat <- sapply(rownames(exp_mat_sub_unique), function(gene) {
  expr <- exp_mat_sub_unique[gene, ]
  model_data <- cbind(expr = expr, marker_scores_df[, c("Astrocytes", "Microglia", "Neuron", "Oligodendrocytes", "Endothelia")])
  fit <- try(lm(expr ~ ., data = model_data), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(rep(NA, length(expr)))
  } else {
    return(residuals(fit))
  }
})
# Convert list to matrix: genes x samples
adjusted_mat <- do.call(rbind, adjusted_mat)  

# calculating group medians
group_medians <- sapply(unique(marker_scores_df$diagnosis), function(group) {
  samples_in_group <- marker_scores_df$id[marker_scores_df$diagnosis == group]
  apply(adjusted_mat[, samples_in_group, drop = FALSE], 1, median, na.rm = TRUE)
})

# Z-score and heatmap
library(pheatmap)
complex <- read.csv("Complex_I_IV.csv")
complex <- complex[1:45, ] 
complex <- complex$elements

group_median_mito <- group_medians[rownames(group_medians) %in% complex, ]

group_median_mito <- t(group_median_mito)
new_order = c("Control","AD","AsymAD")
group_median_mito <- group_median_mito[new_order,]
df <- group_median_mito
df <- df[, order(factor(colnames(df), levels = c(
  "NDUFA10", "NDUFA12", "NDUFA13", "NDUFA2", "NDUFA6", "NDUFA7", "NDUFA8", "NDUFA9", "NDUFAB1",
  "NDUFB11", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFC2",
  "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2",
  "MT-ND3", "MT-ND4", "SDHC", "UQCR10", "UQCRB", "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRH", "UQCRQ",
  "COX4I1", "COX5A", "COX5B", "COX6B1", "COX7C", "CYC1", "MT-CO3"
)))]

plot <- pheatmap(df,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 20, 
         main = "Cell-type corrected median mitochondrial protein expression",
         color = colorRampPalette(c("#00008B", "white", "#8B0000"))(50))  
plot
ggsave("1_Complex_I_IV_heatmap.png", plot = plot , dpi = 600, bg = "transparent", limitsize = F)

# Loading mitochondrial gene list

library(Seurat)
setwd("/depot/pbaloni/data/Lab_members/Purba_Mandal/snRNA_seq_analysis/Dammer_classification/Objects")
mito <- readLines("/depot/pbaloni/data/Lab_members/Purba_Mandal/PFC_MIT_data/DEG_results/mito.txt")

# Directory containing .rds objects
obj_dir <- "/depot/pbaloni/data/Lab_members/Purba_Mandal/snRNA_seq_analysis/Dammer_classification/Objects"

# A directory to store results 
results_dir <- file.path("/depot/pbaloni/data/Lab_members/Purba_Mandal/snRNA_seq_analysis/Dammer_classification/Objects", "DEG_results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Listing all .rds files in the objects directory
obj_files <- list.files(obj_dir, pattern = "\\.rds$", full.names = TRUE)

# Looping over each .rds file
for (obj_file in obj_files) {
  obj <- readRDS(obj_file)
  
  # Extracting object name for file naming
  obj_name <- tools::file_path_sans_ext(basename(obj_file))
  
  # Identifying unique cell subtypes
  Mic_subtypes <- unique(obj$cell_type_high_resolution)
  # Sub-setting to mitochondrial genes
  obj <- subset(obj, features = mito)
  
  # Looping over each cell subtype
  for (subtype in Mic_subtypes) {
    message("Running DEG for object: ", obj_name, " - subtype: ", subtype)
    
    # Subsetting the object for this subtype
    subtype_obj <- subset(obj, subset = cell_type_high_resolution == subtype)
    Idents(subtype_obj) <- "diagnosis"
    
    # differential expression
    deg <- tryCatch({
      Seurat::FindMarkers(
        subtype_obj,
        ident.1 = "AsymAD",
        ident.2 = "AD",
        test.use = "MAST",
        logfc.threshold = 0.1,
        min.pct = 0.1,
        return.thresh = 0.05,
        latent.vars = c("pmi", "age_death", "sex")
      )
    }, error = function(e) {
      message("Skipping ", subtype, " in object ", obj_name, " due to error: ", e$message)
      return(NULL)
    })
    
    if (!is.null(deg) && nrow(deg) > 0) {
      # Adding columns for gene and subtype
      deg$gene <- rownames(deg)
      deg$subtype <- subtype
      
      # subtype names 
      clean_subtype <- gsub("[/ ]", "_", subtype)
      
      # Saving all DEGs for this subtype
      output_name_all <- paste0(obj_name, "_", clean_subtype, "_DEG_MAST.csv")
      write.csv(deg, file.path(results_dir, output_name_all), row.names = FALSE)
    }
  }
}


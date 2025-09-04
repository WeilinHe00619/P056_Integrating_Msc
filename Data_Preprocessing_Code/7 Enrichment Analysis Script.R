# ===============================================================================
# Enrichment Analysis Script - Tailored to Real Data Structure
# Based on actual data analysis from uploaded files
# ===============================================================================

cat("=== ENRICHMENT ANALYSIS - REAL DATA VERSION ===\n")

# Load required packages
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(readxl)
  library(dplyr)
})

# Fix potential namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# Define paths
INPUT_PATH <- "/Users/heweilin/Desktop/P056_Code_2/Raw_Data"
OUTPUT_PATH <- "/Users/heweilin/Desktop/P056_Code_2/Processed_Data"

# Create output directory
dir.create(OUTPUT_PATH, recursive = TRUE, showWarnings = FALSE)
cat("Input directory:", INPUT_PATH, "\n")
cat("Output directory:", OUTPUT_PATH, "\n")

# ===============================================================================
# Read data using correct file paths
# ===============================================================================

cat("\n=== READING DATA ===\n")

# Define input file paths
dna_mrna_file <- file.path(INPUT_PATH, "9Integrated_DNA_mRNA.xlsx")
mirna_mrna_file <- file.path(INPUT_PATH, "10Integrated_miRNA_mRNA.xlsx")

# Check if files exist
if (!file.exists(dna_mrna_file)) {
  stop("DNA methylation file not found: ", dna_mrna_file)
}
if (!file.exists(mirna_mrna_file)) {
  stop("miRNA file not found: ", mirna_mrna_file)
}

# Read DNA methylation + mRNA data
cat("Reading DNA methylation data from:", basename(dna_mrna_file), "\n")
# Sheet: Hypo_Up (147 rows) - HYPO methylation + UP expression
hypo_up_data <- read_excel(dna_mrna_file, sheet = "Hypo_Up")
cat("  Hypo_Up:", nrow(hypo_up_data), "rows\n")

# Sheet: Hyper_Down (106 rows) - HYPER methylation + DOWN expression  
hyper_down_data <- read_excel(dna_mrna_file, sheet = "Hyper_Down")
cat("  Hyper_Down:", nrow(hyper_down_data), "rows\n")

# Read miRNA + mRNA data
cat("Reading miRNA data from:", basename(mirna_mrna_file), "\n")
# Sheet: mirUP_mrnaDown (537 rows) - miRNA UP + mRNA DOWN
mir_up_mrna_down_data <- read_excel(mirna_mrna_file, sheet = "mirUP_mrnaDown")
cat("  mirUP_mrnaDown:", nrow(mir_up_mrna_down_data), "rows\n")

# Sheet: mirDown_mrnaUP (544 rows) - miRNA DOWN + mRNA UP
mir_down_mrna_up_data <- read_excel(mirna_mrna_file, sheet = "mirDown_mrnaUP")
cat("  mirDown_mrnaUP:", nrow(mir_down_mrna_up_data), "rows\n")

# ===============================================================================
# Extract gene symbols and convert to Entrez IDs
# ===============================================================================

cat("\n=== EXTRACTING GENE LISTS ===\n")

# Function to extract unique gene symbols and convert to Entrez IDs
extract_entrez_ids <- function(data, symbol_column, group_name) {
  cat("Processing", group_name, "...\n")
  
  # Check if column exists
  if (!symbol_column %in% names(data)) {
    cat("  ERROR: Column", symbol_column, "not found\n")
    cat("  Available columns:", paste(names(data)[1:10], collapse = ", "), "\n")
    return(character(0))
  }
  
  # Extract unique symbols
  symbols <- unique(data[[symbol_column]][!is.na(data[[symbol_column]])])
  symbols <- symbols[symbols != ""]  # Remove empty strings
  cat("  Found", length(symbols), "unique gene symbols\n")
  
  if (length(symbols) == 0) {
    cat("  No valid symbols found\n")
    return(character(0))
  }
  
  # Convert to Entrez IDs
  tryCatch({
    conversion_result <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    entrez_ids <- unique(conversion_result$ENTREZID)
    cat("  Converted", length(symbols), "symbols to", length(entrez_ids), "Entrez IDs\n")
    
    # Show conversion rate
    conversion_rate <- round(length(entrez_ids) / length(symbols) * 100, 1)
    cat("  Conversion rate:", conversion_rate, "%\n")
    
    return(entrez_ids)
  }, error = function(e) {
    cat("  ERROR in symbol conversion:", e$message, "\n")
    return(character(0))
  })
}

# Extract gene lists using the correct column names from real data
gene_groups <- list()

# DNA methylation data - use "SYMBOL" column (index 37 in headers)
gene_groups[["Hypo_Up"]] <- extract_entrez_ids(hypo_up_data, "SYMBOL", "Hypo_Up")
gene_groups[["Hyper_Down"]] <- extract_entrez_ids(hyper_down_data, "SYMBOL", "Hyper_Down")

# miRNA data - use "Target.Gene" column (index 2 in headers)
gene_groups[["mirUP_mrnaDown"]] <- extract_entrez_ids(mir_up_mrna_down_data, "Target.Gene", "mirUP_mrnaDown")
gene_groups[["mirDown_mrnaUP"]] <- extract_entrez_ids(mir_down_mrna_up_data, "Target.Gene", "mirDown_mrnaUP")

# Remove empty groups
gene_groups <- gene_groups[sapply(gene_groups, length) > 0]

cat("\nSummary of gene groups:\n")
for (group_name in names(gene_groups)) {
  cat("  ", group_name, ":", length(gene_groups[[group_name]]), "genes\n")
}

if (length(gene_groups) == 0) {
  cat("ERROR: No valid gene groups found. Check column names and data.\n")
  stop("Cannot proceed without gene lists")
}

# ===============================================================================
# Perform enrichment analysis
# ===============================================================================

cat("\n=== ENRICHMENT ANALYSIS ===\n")

# Function to perform GO and KEGG enrichment
perform_enrichment <- function(entrez_ids, group_name) {
  cat("\nAnalyzing", group_name, "(", length(entrez_ids), "genes)...\n")
  
  if (length(entrez_ids) < 3) {
    cat("  Too few genes for enrichment analysis (minimum 3 required)\n")
    return(list())
  }
  
  results <- list()
  
  # GO enrichment analysis
  cat("  Running GO enrichment...\n")
  tryCatch({
    go_result <- enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "ALL",  # BP, MF, CC
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      minGSSize = 3,
      maxGSSize = 500
    )
    
    if (!is.null(go_result) && nrow(go_result@result) > 0) {
      # Convert gene IDs to symbols for readability
      go_readable <- setReadable(go_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      go_df <- as.data.frame(go_readable)
      
      # Add metadata
      go_df$Group <- group_name
      go_df$Analysis <- "GO"
      
      results$GO <- go_df
      cat("    GO results:", nrow(go_df), "significant terms\n")
    } else {
      cat("    GO: No significant terms found\n")
    }
  }, error = function(e) {
    cat("    GO analysis failed:", e$message, "\n")
  })
  
  # KEGG pathway enrichment
  cat("  Running KEGG enrichment...\n")
  tryCatch({
    kegg_result <- enrichKEGG(
      gene = entrez_ids,
      organism = "hsa",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 3,
      maxGSSize = 500
    )
    
    if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
      # Convert gene IDs to symbols for readability
      kegg_readable <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      kegg_df <- as.data.frame(kegg_readable)
      
      # Add ONTOLOGY column for KEGG (to match GO format)
      kegg_df$ONTOLOGY <- "KEGG"
      
      # Add metadata
      kegg_df$Group <- group_name
      kegg_df$Analysis <- "KEGG"
      
      results$KEGG <- kegg_df
      cat("    KEGG results:", nrow(kegg_df), "significant pathways\n")
    } else {
      cat("    KEGG: No significant pathways found\n")
    }
  }, error = function(e) {
    cat("    KEGG analysis failed:", e$message, "\n")
  })
  
  return(results)
}

# Run enrichment for all gene groups
all_results <- list()
total_terms <- 0

for (group_name in names(gene_groups)) {
  group_results <- perform_enrichment(gene_groups[[group_name]], group_name)
  all_results[[group_name]] <- group_results
  
  # Count terms for this group
  group_terms <- sum(sapply(group_results, function(x) if(is.data.frame(x)) nrow(x) else 0))
  total_terms <- total_terms + group_terms
}

cat("\nTotal enrichment terms found:", total_terms, "\n")

# ===============================================================================
# Combine and save results
# ===============================================================================

cat("\n=== SAVING RESULTS ===\n")

# Function to standardize column structure
standardize_columns <- function(df, analysis_type) {
  # Convert tibble to data.frame to avoid method dispatch issues
  if (inherits(df, "tbl_df")) {
    df <- as.data.frame(df)
  }
  
  # Define required columns
  required_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", 
                     "p.adjust", "qvalue", "geneID", "Count", "ONTOLOGY", 
                     "Group", "Analysis")
  
  # Add missing columns with appropriate default values
  for (col in required_cols) {
    if (!col %in% names(df)) {
      if (col == "ONTOLOGY" && analysis_type == "KEGG") {
        df[[col]] <- "KEGG"
      } else if (col == "ONTOLOGY" && analysis_type == "GO") {
        # ONTOLOGY should already exist for GO results
        df[[col]] <- ifelse(is.na(df[[col]]), "GO", df[[col]])
      } else if (col %in% c("Group", "Analysis")) {
        # These should already be added
        next
      } else {
        # Add default values for missing columns
        df[[col]] <- NA
      }
    }
  }
  
  # Ensure columns are in the right order and convert to data.frame
  available_cols <- intersect(required_cols, names(df))
  df <- df[, available_cols, drop = FALSE]
  df <- as.data.frame(df)
  
  return(df)
}

# Combine all results into one data frame
combined_results <- data.frame()

for (group_name in names(all_results)) {
  group_data <- all_results[[group_name]]
  
  for (analysis_type in names(group_data)) {
    df <- group_data[[analysis_type]]
    if (is.data.frame(df) && nrow(df) > 0) {
      # Standardize column structure
      df_standardized <- standardize_columns(df, analysis_type)
      
      # Use dplyr::bind_rows for safer combining
      if (nrow(combined_results) == 0) {
        combined_results <- df_standardized
      } else {
        combined_results <- bind_rows(combined_results, df_standardized)
      }
    }
  }
}

cat("Combined results:", nrow(combined_results), "terms\n")

# Save complete results
complete_file <- file.path(OUTPUT_PATH, "2_Enrichment_Complete_Results.csv")
if (nrow(combined_results) > 0) {
  write.csv(combined_results, complete_file, row.names = FALSE)
  cat("✓ Complete results saved:", basename(complete_file), "\n")
} else {
  # Create empty file with proper structure
  empty_df <- data.frame(
    ID = character(0), Description = character(0), GeneRatio = character(0),
    BgRatio = character(0), pvalue = numeric(0), p.adjust = numeric(0),
    qvalue = numeric(0), geneID = character(0), Count = integer(0),
    ONTOLOGY = character(0), Group = character(0), Analysis = character(0)
  )
  write.csv(empty_df, complete_file, row.names = FALSE)
  cat("✓ Empty complete results file created\n")
}

# Create visualization table (use_pathway)
viz_file <- file.path(OUTPUT_PATH, "2_Enrichment_use_pathway.csv")
if (nrow(combined_results) > 0) {
  # Select top pathways for visualization
  use_pathway <- combined_results %>%
    dplyr::group_by(ONTOLOGY) %>%
    dplyr::arrange(p.adjust, desc(Count)) %>%
    dplyr::slice_head(n = 5) %>%  # Top 5 per ontology
    dplyr::group_by(p.adjust) %>%
    dplyr::slice_max(Count, n = 1, with_ties = FALSE) %>%  # If same p-value, take highest count
    dplyr::ungroup() %>%
    dplyr::arrange(ONTOLOGY, p.adjust) %>%
    dplyr::mutate(index = dplyr::row_number()) %>%
    dplyr::select(Description, p.adjust, geneID, Count, ONTOLOGY, index, Group, Analysis)
  
  write.csv(use_pathway, viz_file, row.names = FALSE)
  cat("✓ Visualization table saved:", basename(viz_file), "(", nrow(use_pathway), "pathways)\n")
} else {
  # Create empty visualization table
  empty_viz <- data.frame(
    Description = character(0), p.adjust = numeric(0), geneID = character(0),
    Count = integer(0), ONTOLOGY = character(0), index = integer(0),
    Group = character(0), Analysis = character(0)
  )
  write.csv(empty_viz, viz_file, row.names = FALSE)
  cat("✓ Empty visualization table created\n")
}

# Save individual group results
for (group_name in names(all_results)) {
  group_data <- all_results[[group_name]]
  
  if (length(group_data) > 0) {
    # Combine GO and KEGG results for this group
    group_combined <- data.frame()
    
    for (analysis_type in names(group_data)) {
      df <- group_data[[analysis_type]]
      if (is.data.frame(df) && nrow(df) > 0) {
        # Standardize column structure for individual group files too
        df_standardized <- standardize_columns(df, analysis_type)
        
        if (nrow(group_combined) == 0) {
          group_combined <- df_standardized
        } else {
          group_combined <- bind_rows(group_combined, df_standardized)
        }
      }
    }
    
    if (nrow(group_combined) > 0) {
      group_file <- file.path(OUTPUT_PATH, paste0("2_Enrichment_", group_name, ".csv"))
      write.csv(group_combined, group_file, row.names = FALSE)
      cat("✓ Group", group_name, "saved:", basename(group_file), "(", nrow(group_combined), "terms)\n")
    }
  }
}

# ===============================================================================
# Final verification and summary
# ===============================================================================

cat("\n=== FINAL VERIFICATION ===\n")

# Check all created files
created_files <- list.files(OUTPUT_PATH, pattern = "^2_Enrichment.*\\.csv$", full.names = TRUE)

if (length(created_files) > 0) {
  cat("Created files:\n")
  for (file_path in created_files) {
    file_name <- basename(file_path)
    file_size <- round(file.info(file_path)$size / 1024, 1)
    cat(sprintf("  ✓ %s (%.1f KB)\n", file_name, file_size))
  }
} else {
  cat("⚠️  No files found in output directory\n")
}

# Summary report
cat("\n" + paste(rep("=", 70), collapse = "") + "\n")
cat("ENRICHMENT ANALYSIS SUMMARY\n")
cat(paste(rep("=", 70), collapse = "") + "\n")
cat("Data processed:\n")
cat("  - Hypo_Up:", length(gene_groups[["Hypo_Up"]] %||% character(0)), "genes\n")
cat("  - Hyper_Down:", length(gene_groups[["Hyper_Down"]] %||% character(0)), "genes\n")
cat("  - mirUP_mrnaDown:", length(gene_groups[["mirUP_mrnaDown"]] %||% character(0)), "genes\n")
cat("  - mirDown_mrnaUP:", length(gene_groups[["mirDown_mrnaUP"]] %||% character(0)), "genes\n")
cat("Total enrichment terms:", nrow(combined_results), "\n")
cat("Files created:", length(created_files), "\n")
cat("Input directory:", INPUT_PATH, "\n")
cat("Output directory:", OUTPUT_PATH, "\n")
cat(paste(rep("=", 70), collapse = "") + "\n")

if (length(created_files) > 0) {
  cat("✅ SUCCESS: Enrichment analysis completed successfully!\n")
} else {
  cat("❌ WARNING: No output files were created\n")
}

cat("=== ANALYSIS COMPLETE ===\n")
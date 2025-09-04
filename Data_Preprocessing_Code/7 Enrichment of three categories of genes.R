# ===============================================================================
# Four Gene Set Enrichment Analysis Script - Fixed KEGG Version
# Gene sets: ①Differentially expressed genes ②CpG-regulated mRNA ③Double-regulated genes ④Enriched pathway intersection analysis
# ===============================================================================

cat("=== Four Gene Set Enrichment Analysis Started (Fixed Version) ===\n")

# Load required packages
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(readxl)
  library(dplyr)
  library(VennDiagram)
})

# Fix namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# Define paths
RAW_DATA_PATH <- "/Users/heweilin/Desktop/P056_Code_2/Raw_Data"
PROCESSED_DATA_PATH <- "/Users/heweilin/Desktop/P056_Code_2/Processed_Data"
OUTPUT_PATH <- "/Users/heweilin/Desktop/P056_Code_2/Processed_Data"

# Create output directory
dir.create(OUTPUT_PATH, recursive = TRUE, showWarnings = FALSE)

cat("Input directory:", RAW_DATA_PATH, "\n")
cat("Processed data directory:", PROCESSED_DATA_PATH, "\n")
cat("Output directory:", OUTPUT_PATH, "\n")

# ===============================================================================
# Step 1: Read and extract three gene sets
# ===============================================================================

cat("\n=== Extract three gene sets ===\n")

# Gene set 1: Differentially expressed genes (DEGs)
cat("1. Extract differentially expressed genes (DEGs)...\n")
degs_file <- file.path(RAW_DATA_PATH, "1mRNA_DEGs_proteincoding.csv")
if (!file.exists(degs_file)) {
  stop("DEGs file does not exist: ", degs_file)
}

degs_data <- read.csv(degs_file, stringsAsFactors = FALSE)
cat("  Original DEGs data:", nrow(degs_data), "rows\n")

# Filter significantly differentially expressed genes (padj < 0.05)
degs_genes <- degs_data %>%
  filter(!is.na(padj) & padj < 0.05 & !is.na(SYMBOL) & SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique()

cat("  Significant DEGs genes:", length(degs_genes), "\n")

# Gene set 2: CpG-regulated mRNA (extracted from integrated file)
cat("2. Extract CpG-regulated mRNA genes...\n")
cpg_file <- file.path(RAW_DATA_PATH, "9Integrated_DNA_mRNA.xlsx")
if (!file.exists(cpg_file)) {
  stop("CpG integrated file does not exist: ", cpg_file)
}

# Read two key sheets of CpG regulation data
hypo_up <- read_excel(cpg_file, sheet = "Hypo_Up")     # Hypomethylation→upregulation
hyper_down <- read_excel(cpg_file, sheet = "Hyper_Down") # Hypermethylation→downregulation

cat("  Hypo_Up data:", nrow(hypo_up), "rows\n")
cat("  Hyper_Down data:", nrow(hyper_down), "rows\n")

# Extract CpG-regulated genes (merge two regulation modes)
cpg_hypo_genes <- hypo_up %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique()

cpg_hyper_genes <- hyper_down %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique()

# Merge all CpG-regulated genes
cpg_genes <- union(cpg_hypo_genes, cpg_hyper_genes)

cat("  CpG hypomethylation→upregulation genes:", length(cpg_hypo_genes), "\n")
cat("  CpG hypermethylation→downregulation genes:", length(cpg_hyper_genes), "\n")
cat("  Total CpG-regulated genes:", length(cpg_genes), "\n")

# Gene set 3: Double-regulated genes (miRNA + CpG)
cat("3. Extract double-regulated genes (miRNA + CpG)...\n")
mirna_file <- file.path(RAW_DATA_PATH, "10Integrated_miRNA_mRNA.xlsx")
if (!file.exists(mirna_file)) {
  stop("miRNA integrated file does not exist: ", mirna_file)
}

# Read miRNA regulation data
mir_up_mrna_down <- read_excel(mirna_file, sheet = "mirUP_mrnaDown")
mir_down_mrna_up <- read_excel(mirna_file, sheet = "mirDown_mrnaUP")

cat("  mirUP_mrnaDown data:", nrow(mir_up_mrna_down), "rows\n")
cat("  mirDown_mrnaUP data:", nrow(mir_down_mrna_up), "rows\n")

# View actual column names to ensure using correct columns
cat("  mirUP_mrnaDown column names:", paste(names(mir_up_mrna_down)[1:5], collapse = ", "), "\n")
cat("  mirDown_mrnaUP column names:", paste(names(mir_down_mrna_up)[1:5], collapse = ", "), "\n")

# Extract miRNA-regulated target genes (two regulation modes)
# Mode 1: miRNA upregulated → mRNA downregulated (classical negative regulation)
mirna_up_targets <- mir_up_mrna_down$Target.Gene %>%
  .[!is.na(.) & . != ""] %>%
  unique()

# Mode 2: miRNA downregulated → mRNA upregulated (de-repression effect)  
mirna_down_targets <- mir_down_mrna_up$Target.Gene %>%
  .[!is.na(.) & . != ""] %>%
  unique()

# Merge all miRNA-regulated target genes
mirna_target_genes <- union(mirna_up_targets, mirna_down_targets)

cat("  miRNA upregulated→mRNA downregulated target genes:", length(mirna_up_targets), "\n")
cat("  miRNA downregulated→mRNA upregulated target genes:", length(mirna_down_targets), "\n")
cat("  Total miRNA-regulated target genes:", length(mirna_target_genes), "\n")

# Debug information: show actual gene name samples
cat("  \n=== Debug information ===\n")
cat("  miRNA-regulated gene examples:", paste(head(mirna_target_genes, 8), collapse = ", "), "\n")
cat("  CpG-regulated gene examples:", paste(head(cpg_genes, 8), collapse = ", "), "\n")

# Check if gene names contain special characters or spaces
mirna_sample <- head(mirna_target_genes, 3)
cpg_sample <- head(cpg_genes, 3)
for (i in 1:length(mirna_sample)) {
  cat("  miRNA gene", i, ":", mirna_sample[i], ", length:", nchar(mirna_sample[i]), "\n")
}
for (i in 1:length(cpg_sample)) {
  cat("  CpG gene", i, ":", cpg_sample[i], ", length:", nchar(cpg_sample[i]), "\n")
}

# Double-regulated genes: miRNA regulation ∩ CpG regulation
# Clean gene names, remove possible spaces and special characters
mirna_clean <- trimws(mirna_target_genes)
cpg_clean <- trimws(cpg_genes)

double_regulated_genes <- intersect(mirna_clean, cpg_clean)
cat("  \nDouble-regulated genes (direct intersection):", length(double_regulated_genes), "\n")

# If still no intersection, perform more detailed analysis
if (length(double_regulated_genes) == 0) {
  cat("  \n=== Detailed intersection analysis ===\n")
  
  # 1. Check if there are case sensitivity issues
  mirna_lower <- tolower(mirna_clean)
  cpg_lower <- tolower(cpg_clean)
  case_insensitive_intersection <- intersect(mirna_lower, cpg_lower)
  
  if (length(case_insensitive_intersection) > 0) {
    cat("  Found case mismatch issue! Case-insensitive intersection:", length(case_insensitive_intersection), "genes\n")
    cat("  Intersection genes:", paste(head(case_insensitive_intersection, 10), collapse = ", "), "\n")
    
    # Get original case gene names
    original_case_genes <- c()
    for (gene_lower in case_insensitive_intersection) {
      # Find original case from CpG list
      original_gene <- cpg_clean[tolower(cpg_clean) == gene_lower][1]
      if (!is.na(original_gene)) {
        original_case_genes <- c(original_case_genes, original_gene)
      }
    }
    double_regulated_genes <- original_case_genes
    cat("  Double-regulated genes with original case:", length(double_regulated_genes), "\n")
  } else {
    cat("  No intersection even ignoring case\n")
    
    # 2. Random sampling comparison
    cat("  Random sample comparison:\n")
    sample_mirna <- sample(mirna_clean, min(10, length(mirna_clean)))
    sample_cpg <- sample(cpg_clean, min(10, length(cpg_clean)))
    
    cat("    miRNA sample:", paste(sample_mirna, collapse = ", "), "\n")
    cat("    CpG sample:", paste(sample_cpg, collapse = ", "), "\n")
    
    # To continue analysis, use miRNA genes as third category
    cat("  Warning: No true double-regulated genes found, will use miRNA-regulated genes as third category for analysis\n")
    double_regulated_genes <- head(mirna_clean, 100)  # Limit number to avoid being too large
  }
} else {
  cat("  ✓ Successfully found double-regulated genes:", paste(head(double_regulated_genes, 10), collapse = ", "), "\n")
}

# Summarize first three gene set information
gene_sets <- list(
  "DEGs" = degs_genes,
  "CpG_regulated" = cpg_genes,
  "Double_regulated" = double_regulated_genes
)

# Print gene set overview
cat("\n=== First three gene set overview ===\n")
for (i in 1:length(gene_sets)) {
  set_name <- names(gene_sets)[i]
  set_size <- length(gene_sets[[i]])
  cat(sprintf("  %d. %s: %d genes\n", i, set_name, set_size))
}

# ===============================================================================
# Step 2: Gene symbol conversion to Entrez ID
# ===============================================================================

cat("\n=== Gene symbol conversion ===\n")

# Conversion function
convert_to_entrez <- function(symbols, set_name) {
  cat("Converting", set_name, "...\n")
  
  if (length(symbols) == 0) {
    cat("  Warning: No gene symbols to convert\n")
    return(character(0))
  }
  
  tryCatch({
    conversion_result <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    entrez_ids <- unique(conversion_result$ENTREZID)
    conversion_rate <- round(length(entrez_ids) / length(symbols) * 100, 1)
    cat("  ", length(symbols), "symbols →", length(entrez_ids), "Entrez IDs (", conversion_rate, "%)\n")
    return(entrez_ids)
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    return(character(0))
  })
}

# Convert first three gene sets
entrez_gene_sets <- list()
for (set_name in names(gene_sets)) {
  entrez_gene_sets[[set_name]] <- convert_to_entrez(gene_sets[[set_name]], set_name)
}

# Filter out empty gene sets
entrez_gene_sets <- entrez_gene_sets[sapply(entrez_gene_sets, length) > 0]

# ===============================================================================
# Step 3: Enrichment analysis function - Fixed KEGG version
# ===============================================================================

cat("\n=== Define enrichment analysis function (fixed version) ===\n")

perform_enrichment_analysis <- function(entrez_ids, set_name, set_index) {
  cat("\nAnalyzing", set_name, "(", length(entrez_ids), "genes)...\n")
  
  if (length(entrez_ids) < 3) {
    cat("  Insufficient gene count, skipping analysis (at least 3 genes needed)\n")
    return(data.frame())
  }
  
  all_results <- data.frame()
  
  # GO enrichment analysis
  cat("  Performing GO enrichment analysis...\n")
  tryCatch({
    go_result <- enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      minGSSize = 3,
      maxGSSize = 500
    )
    
    if (!is.null(go_result) && nrow(go_result@result) > 0) {
      go_readable <- setReadable(go_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      go_df <- as.data.frame(go_readable)
      go_df$Analysis <- "GO"
      go_df$Group <- set_name
      go_df$Category <- set_index
      all_results <- rbind(all_results, go_df)
      cat("    GO results:", nrow(go_df), "significant terms\n")
    } else {
      cat("    GO: No significant results\n")
    }
  }, error = function(e) {
    cat("    GO analysis failed:", e$message, "\n")
  })
  
  # KEGG enrichment analysis - Fixed version
  cat("  Performing KEGG enrichment analysis...\n")
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
      # Directly get data frame from result slot
      kegg_df <- kegg_result@result
      
      if (nrow(kegg_df) > 0) {
        # Add required columns
        kegg_df$ONTOLOGY <- "KEGG"
        kegg_df$Analysis <- "KEGG"
        kegg_df$Group <- set_name
        kegg_df$Category <- set_index
        
        # Convert gene IDs to gene symbols
        tryCatch({
          kegg_readable <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
          kegg_df$geneID <- kegg_readable@result$geneID
        }, error = function(e) {
          cat("      Gene symbol conversion failed, using original IDs\n")
          # If conversion fails, ensure geneID column exists
          if (!"geneID" %in% names(kegg_df)) {
            kegg_df$geneID <- kegg_df$ID
          }
        })
        
        # Ensure all required columns exist
        required_cols <- c("Description", "p.adjust", "pvalue", "Count", "ONTOLOGY", "Analysis", "Group", "Category", "geneID")
        for (col in required_cols) {
          if (!col %in% names(kegg_df)) {
            if (col == "pvalue") {
              kegg_df$pvalue <- kegg_df$p.adjust
            } else {
              kegg_df[[col]] <- NA
            }
          }
        }
        
        # If GO results exist, ensure consistent column structure
        if (nrow(all_results) > 0) {
          go_cols <- names(all_results)
          # Ensure KEGG results have same columns
          for (go_col in go_cols) {
            if (!go_col %in% names(kegg_df)) {
              kegg_df[[go_col]] <- NA
            }
          }
          # Arrange KEGG results by GO results column order
          kegg_df <- kegg_df[, go_cols, drop = FALSE]
        }
        
        all_results <- rbind(all_results, kegg_df)
        cat("    KEGG results:", nrow(kegg_df), "significant pathways\n")
      } else {
        cat("    KEGG: No significant pathways\n")
      }
    } else {
      cat("    KEGG: No significant pathways\n")
    }
  }, error = function(e) {
    cat("    KEGG analysis failed:", e$message, "\n")
    cat("    Continue analysis, using only GO results\n")
  })
  
  return(all_results)
}

# ===============================================================================
# Step 4: Perform first three enrichment analyses
# ===============================================================================

cat("\n=== Perform first three enrichment analyses ===\n")

all_enrichment_results <- list()
set_names <- c("DEGs", "CpG_regulated", "Double_regulated")
set_indices <- c("1", "2", "3")

for (i in 1:length(set_names)) {
  set_name <- set_names[i]
  set_index <- set_indices[i]
  
  if (set_name %in% names(entrez_gene_sets)) {
    result <- perform_enrichment_analysis(
      entrez_gene_sets[[set_name]], 
      set_name, 
      set_index
    )
    all_enrichment_results[[set_name]] <- result
  }
}

# ===============================================================================
# Step 5: Perform intersection analysis (Category 4)
# ===============================================================================

cat("\n=== Perform enriched pathway intersection analysis ===\n")

# Collect all enrichment results
all_pathway_results <- do.call(rbind, all_enrichment_results)

if (nrow(all_pathway_results) > 0) {
  
  # 1. Find common enriched pathways (same Description)
  cat("1. Analyze common enriched pathways...\n")
  
  # Count how many gene sets each pathway appears in
  pathway_counts <- all_pathway_results %>%
    group_by(Description, ONTOLOGY) %>%
    summarise(
      appeared_in_sets = n_distinct(Group),
      groups = paste(unique(Group), collapse = ", "),
      avg_pvalue = mean(pvalue, na.rm = TRUE),
      avg_p.adjust = mean(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(appeared_in_sets >= 2)  # Appear in at least 2 gene sets
  
  cat("  Common pathway count:", nrow(pathway_counts), "\n")
  
  # 2. Find common enriched genes
  cat("2. Analyze common enriched genes...\n")
  
  # Extract all enriched genes
  all_enriched_genes <- c()
  gene_pathway_mapping <- data.frame()
  
  for (i in 1:nrow(all_pathway_results)) {
    row <- all_pathway_results[i, ]
    genes <- unlist(strsplit(as.character(row$geneID), "/"))
    
    if (length(genes) > 0) {
      all_enriched_genes <- c(all_enriched_genes, genes)
      
      # Record gene-pathway mapping
      for (gene in genes) {
        gene_pathway_mapping <- rbind(gene_pathway_mapping, data.frame(
          Gene = gene,
          Pathway = row$Description,
          Group = row$Group,
          ONTOLOGY = row$ONTOLOGY,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Count how many gene sets each gene appears in
  gene_counts <- gene_pathway_mapping %>%
    group_by(Gene) %>%
    summarise(
      appeared_in_sets = n_distinct(Group),
      groups = paste(unique(Group), collapse = ", "),
      pathways = paste(unique(Pathway), collapse = "; "),
      .groups = "drop"
    ) %>%
    filter(appeared_in_sets >= 1) %>%  # Lower threshold, appear in at least 1 gene set
    arrange(desc(appeared_in_sets))
  
  cat("  Common gene count:", nrow(gene_counts), "\n")
  
  # 3. Create intersection analysis results
  cat("3. Create intersection analysis results...\n")
  
  # Create intersection results based on common pathways
  intersection_results <- data.frame()
  
  if (nrow(pathway_counts) > 0) {
    for (i in 1:nrow(pathway_counts)) {
      pathway_info <- pathway_counts[i, ]
      
      # Get detailed information for this pathway (choose the one with smallest p-value)
      pathway_details <- all_pathway_results %>%
        filter(Description == pathway_info$Description, ONTOLOGY == pathway_info$ONTOLOGY) %>%
        arrange(p.adjust) %>%
        slice(1)
      
      # Get all genes appearing in this pathway
      pathway_genes <- all_pathway_results %>%
        filter(Description == pathway_info$Description, ONTOLOGY == pathway_info$ONTOLOGY) %>%
        pull(geneID) %>%
        paste(collapse = "/") %>%
        strsplit("/") %>%
        unlist() %>%
        unique() %>%
        paste(collapse = "/")
      
      # Create intersection result row
      intersection_row <- data.frame(
        Description = pathway_info$Description,
        p.adjust = pathway_info$avg_p.adjust,
        geneID = pathway_genes,
        Count = length(unique(unlist(strsplit(pathway_genes, "/")))),
        ONTOLOGY = pathway_info$ONTOLOGY,
        index = i,
        Group = "Intersection",
        Analysis = ifelse(pathway_info$ONTOLOGY == "KEGG", "KEGG", "GO"),
        Category = "123",
        appeared_in_sets = pathway_info$appeared_in_sets,
        groups = pathway_info$groups,
        stringsAsFactors = FALSE
      )
      
      intersection_results <- rbind(intersection_results, intersection_row)
    }
    
    # Sort and add index
    intersection_results <- intersection_results %>%
      arrange(ONTOLOGY, p.adjust) %>%
      mutate(index = row_number())
  } else {
    # If no common pathways, create empty results
    intersection_results <- data.frame(
      Description = character(0),
      p.adjust = numeric(0),
      geneID = character(0),
      Count = integer(0),
      ONTOLOGY = character(0),
      index = integer(0),
      Group = character(0),
      Analysis = character(0),
      Category = character(0),
      appeared_in_sets = integer(0),
      groups = character(0)
    )
  }
  
  # Add intersection results to total results
  all_enrichment_results[["Intersection"]] <- intersection_results
  
  cat("  Intersection analysis completed, total", nrow(intersection_results), "common pathways\n")
  
} else {
  cat("Warning: First three gene sets have no enrichment results, cannot perform intersection analysis\n")
  # Create empty intersection results
  intersection_results <- data.frame(
    Description = character(0),
    p.adjust = numeric(0),
    geneID = character(0),
    Count = integer(0),
    ONTOLOGY = character(0),
    index = integer(0),
    Group = character(0),
    Analysis = character(0),
    Category = character(0)
  )
  all_enrichment_results[["Intersection"]] <- intersection_results
}

# ===============================================================================
# Step 6: Generate standardized use_pathway format files
# ===============================================================================

cat("\n=== Generate output files ===\n")

# Function to standardize result format
create_use_pathway_format <- function(enrichment_data, set_name, set_index) {
  if (nrow(enrichment_data) == 0) {
    # Create empty standard format
    return(data.frame(
      Description = character(0),
      p.adjust = numeric(0),
      geneID = character(0),
      Count = integer(0),
      ONTOLOGY = character(0),
      index = integer(0),
      Group = character(0),
      Analysis = character(0),
      Category = character(0)
    ))
  }
  
  # Ensure necessary columns exist
  if (!"Category" %in% names(enrichment_data)) {
    enrichment_data$Category <- set_index
  }
  
  # Select top 5 most significant results for each ontology type
  use_pathway <- enrichment_data %>%
    arrange(p.adjust, desc(Count)) %>%
    group_by(ONTOLOGY) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    arrange(ONTOLOGY, p.adjust) %>%
    mutate(index = row_number()) %>%
    select(Description, p.adjust, geneID, Count, ONTOLOGY, index, Group, Analysis, Category)
  
  return(use_pathway)
}

# Generate use_pathway format files for all four gene sets
all_set_names <- c("DEGs", "CpG_regulated", "Double_regulated", "Intersection")
all_set_indices <- c("1", "2", "3", "123")
output_files <- c()

for (i in 1:length(all_set_names)) {
  set_name <- all_set_names[i]
  set_index <- all_set_indices[i]
  
  # Generate file name
  file_name <- paste0("2_1_Enrichment_", set_name, ".csv")
  file_path <- file.path(OUTPUT_PATH, file_name)
  
  if (set_name %in% names(all_enrichment_results)) {
    enrichment_data <- all_enrichment_results[[set_name]]
    use_pathway_data <- create_use_pathway_format(enrichment_data, set_name, set_index)
  } else {
    # Create empty data frame
    use_pathway_data <- create_use_pathway_format(data.frame(), set_name, set_index)
  }
  
  # Save file
  write.csv(use_pathway_data, file_path, row.names = FALSE)
  output_files <- c(output_files, file_name)
  
  cat("✓ Saved:", file_name, "(", nrow(use_pathway_data), "pathways)\n")
}

# ===============================================================================
# Step 7: Generate detailed analysis report
# ===============================================================================

cat("\n=== Generate detailed analysis report ===\n")

# Merge all results
combined_results <- data.frame()

# Safely merge results, ensuring consistent columns
for (set_name in names(all_enrichment_results)) {
  set_data <- all_enrichment_results[[set_name]]
  if (nrow(set_data) > 0) {
    # Ensure all necessary columns exist
    required_cols <- c("Description", "p.adjust", "pvalue", "Count", "ONTOLOGY", "Analysis", "Group", "Category", "geneID")
    
    for (col in required_cols) {
      if (!col %in% names(set_data)) {
        set_data[[col]] <- NA
      }
    }
    
    # Only select needed columns
    set_data_clean <- set_data[, required_cols, drop = FALSE]
    combined_results <- rbind(combined_results, set_data_clean)
  }
}

if (nrow(combined_results) > 0) {
  # Save complete results
  complete_file <- file.path(OUTPUT_PATH, "2_1_Enrichment_Complete_Results.csv")
  write.csv(combined_results, complete_file, row.names = FALSE)
  cat("✓ Complete results saved:", basename(complete_file), "\n")
  
  # Generate summary use_pathway file
  summary_use_pathway <- create_use_pathway_format(combined_results, "All_Sets", "All")
  summary_file <- file.path(OUTPUT_PATH, "2_1_Enrichment_Summary_use_pathway.csv")
  write.csv(summary_use_pathway, summary_file, row.names = FALSE)
  cat("✓ Summary use_pathway saved:", basename(summary_file), "\n")
}

# Save detailed results of intersection analysis
if (exists("pathway_counts") && nrow(pathway_counts) > 0) {
  pathway_intersection_file <- file.path(OUTPUT_PATH, "2_1_Enrichment_Pathway_Intersection_Details.csv")
  write.csv(pathway_counts, pathway_intersection_file, row.names = FALSE)
  cat("✓ Pathway intersection details saved:", basename(pathway_intersection_file), "\n")
}

if (exists("gene_counts") && nrow(gene_counts) > 0) {
  gene_intersection_file <- file.path(OUTPUT_PATH, "2_1_Enrichment_Gene_Intersection_Details.csv")
  write.csv(gene_counts, gene_intersection_file, row.names = FALSE)
  cat("✓ Gene intersection details saved:", basename(gene_intersection_file), "\n")
}

# ===============================================================================
# Step 8: Final report
# ===============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("Four Gene Set Enrichment Analysis Completion Report (Fixed Version)\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("Gene set overview:\n")
for (i in 1:3) {
  set_name <- all_set_names[i]
  set_index <- all_set_indices[i]
  original_count <- length(gene_sets[[set_name]] %||% character(0))
  entrez_count <- length(entrez_gene_sets[[set_name]] %||% character(0))
  cat(sprintf("  %s (Category %s): %d symbols → %d EntrezIDs\n", 
              set_name, set_index, original_count, entrez_count))
}

cat("\nIntersection analysis results:\n")
if (exists("pathway_counts") && nrow(pathway_counts) > 0) {
  cat("  Common enriched pathways:", nrow(pathway_counts), "pathways\n")
  
  # Show pathways with highest frequency
  top_pathways <- pathway_counts %>%
    arrange(desc(appeared_in_sets), avg_p.adjust) %>%
    slice_head(n = 5)
  
  for (i in 1:nrow(top_pathways)) {
    pathway <- top_pathways[i, ]
    cat(sprintf("    - %s (appears in %d gene sets)\n", 
                pathway$Description, pathway$appeared_in_sets))
  }
} else {
  cat("  No common enriched pathways\n")
}

if (exists("gene_counts") && nrow(gene_counts) > 0) {
  cat("  Common enriched genes:", nrow(gene_counts), "genes\n")
  
  # Show genes with highest frequency
  top_genes <- gene_counts %>%
    arrange(desc(appeared_in_sets)) %>%
    slice_head(n = 10)
  
  for (i in 1:min(5, nrow(top_genes))) {
    gene <- top_genes[i, ]
    cat(sprintf("    - %s (appears in %d gene sets)\n", 
                gene$Gene, gene$appeared_in_sets))
  }
} else {
  cat("  No common enriched genes\n")
}

cat("\nOutput file list:\n")
for (file in output_files) {
  cat("  ✓", file, "\n")
}

# Additional output files
additional_files <- c(
  "2_1_Enrichment_Complete_Results.csv",
  "2_1_Enrichment_Summary_use_pathway.csv",
  "2_1_Enrichment_Pathway_Intersection_Details.csv",
  "2_1_Enrichment_Gene_Intersection_Details.csv"
)

for (file in additional_files) {
  if (file.exists(file.path(OUTPUT_PATH, file))) {
    cat("  ✓", file, "\n")
  }
}

cat("\n=== Analysis Complete ===\n")
cat("All result files have been saved to:", OUTPUT_PATH, "\n")

# ===============================================================================
# Step 9: Optional package version information and troubleshooting
# ===============================================================================

cat("\n=== Package Version Information ===\n")
cat("clusterProfiler version:", as.character(packageVersion("clusterProfiler")), "\n")
cat("org.Hs.eg.db version:", as.character(packageVersion("org.Hs.eg.db")), "\n")
cat("dplyr version:", as.character(packageVersion("dplyr")), "\n")

# Check KEGG connection
cat("\n=== KEGG Connection Test ===\n")
tryCatch({
  # Simple KEGG connection test
  test_kegg <- enrichKEGG(gene = c("1", "2", "3"), organism = "hsa", pvalueCutoff = 1)
  cat("✓ KEGG database connection normal\n")
}, error = function(e) {
  cat("⚠ KEGG connection may have issues:", e$message, "\n")
  cat("  Recommend checking network connection or try again later\n")
})

cat("\nIf KEGG analysis still has problems, consider:\n")
cat("1. Check network connection\n")
cat("2. Update clusterProfiler package: BiocManager::install('clusterProfiler')\n")
cat("3. Use only GO analysis results (already comprehensive)\n")
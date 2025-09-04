# Load required libraries
library(tidyverse)
library(ComplexHeatmap)
library(corrplot)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(circlize)

# Set file paths
cpg_file <- "/Users/heweilin/Desktop/P056_Code_2/Processed_Data/1_OR_CpG_expr_raw.csv"
mirna_file <- "/Users/heweilin/Desktop/P056_Code_2/Processed_Data/1_OR_miRNA_expr_raw.csv"
mrna_file <- "/Users/heweilin/Desktop/P056_Code_2/Processed_Data/1_OR_mRNA_expr_raw.csv"
adjusted_cpg_file <- "/Users/heweilin/Desktop/P056_Code_2/Processed_Data/3_CpG_expr_adjusted.csv"
adjusted_mirna_file <- "/Users/heweilin/Desktop/P056_Code_2/Processed_Data/3_miRNA_expr_adjusted.csv"
adjusted_mrna_file <- "/Users/heweilin/Desktop/P056_Code_2/Processed_Data/3_mRNA_expr_adjusted.csv"
output_dir <- "/Users/heweilin/Desktop/P056_Code_3/Figure"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Loading data files...\n")

# ================== PART 1: Load Raw Data for Correlation Heatmap ==================

# Read CpG data (Gene_Symbol column contains row names)
cpg_raw <- read.csv(cpg_file, check.names = FALSE, stringsAsFactors = FALSE)
cpg_data <- cpg_raw[, -1]  # Remove Gene_Symbol column
rownames(cpg_data) <- cpg_raw$Gene_Symbol
# Convert all columns to numeric
cpg_data <- data.frame(lapply(cpg_data, function(x) as.numeric(as.character(x))))

# Read miRNA data (sRNA column contains row names)
mirna_raw <- read.csv(mirna_file, check.names = FALSE, stringsAsFactors = FALSE)
mirna_data <- mirna_raw[, -1]  # Remove sRNA column
rownames(mirna_data) <- mirna_raw$sRNA
# Convert all columns to numeric
mirna_data <- data.frame(lapply(mirna_data, function(x) as.numeric(as.character(x))))

# Read mRNA data (first column contains row names)
mrna_raw <- read.csv(mrna_file, check.names = FALSE, stringsAsFactors = FALSE)
mrna_data <- mrna_raw[, -1]  # Remove first column
rownames(mrna_data) <- mrna_raw[, 1]  # Set first column as row names
# Convert all columns to numeric
mrna_data <- data.frame(lapply(mrna_data, function(x) as.numeric(as.character(x))))

# Print data dimensions and value ranges
cat(sprintf("CpG methylation data (raw): %d features x %d samples\n", nrow(cpg_data), ncol(cpg_data)))
cat(sprintf("miRNA expression data (raw): %d features x %d samples\n", nrow(mirna_data), ncol(mirna_data)))
cat(sprintf("mRNA expression data (raw): %d features x %d samples\n", nrow(mrna_data), ncol(mrna_data)))

# Extract base patient IDs (remove suffix letters)
cpg_samples <- colnames(cpg_data)
mirna_samples <- colnames(mirna_data) 
mrna_samples <- colnames(mrna_data)

cpg_base <- gsub("[a-zA-Z]+$", "", cpg_samples)
mirna_base <- gsub("[a-zA-Z]+$", "", mirna_samples)
mrna_base <- gsub("[a-zA-Z]+$", "", mrna_samples)

# Find common base patient IDs
common_base_ids <- Reduce(intersect, list(cpg_base, mirna_base, mrna_base))
cat(sprintf("Common base patient IDs: %d\n", length(common_base_ids)))

# Create matching sample vectors
cpg_matched <- cpg_samples[cpg_base %in% common_base_ids]
mirna_matched <- mirna_samples[mirna_base %in% common_base_ids]
mrna_matched <- mrna_samples[mrna_base %in% common_base_ids]

# Sort by base ID to ensure proper matching
cpg_order <- order(gsub("[a-zA-Z]+$", "", cpg_matched))
mirna_order <- order(gsub("[a-zA-Z]+$", "", mirna_matched))
mrna_order <- order(gsub("[a-zA-Z]+$", "", mrna_matched))

cpg_matched <- cpg_matched[cpg_order]
mirna_matched <- mirna_matched[mirna_order]
mrna_matched <- mrna_matched[mrna_order]

# Subset to matched samples
cpg_common <- cpg_data[, cpg_matched]
mirna_common <- mirna_data[, mirna_matched]
mrna_common <- mrna_data[, mrna_matched]

# Rename columns to have consistent sample names for correlation
common_base_sorted <- common_base_ids[order(common_base_ids)]
colnames(cpg_common) <- common_base_sorted
colnames(mirna_common) <- common_base_sorted
colnames(mrna_common) <- common_base_sorted

cat(sprintf("Matched samples for analysis: %d\n", length(common_base_sorted)))

# Feature filtering: Remove low expression/variance features
cat("\nFiltering low-expression and low-variance features...\n")

# Filter CpG: Remove features with very low variance (< 0.001)
rownames(cpg_common) <- make.names(rownames(cpg_common), unique = TRUE)
cpg_clean <- cpg_common[complete.cases(cpg_common), ]
cpg_vars <- tryCatch({
  apply(cpg_clean, 1, function(x) {
    x_numeric <- as.numeric(x)
    if (all(is.na(x_numeric)) || length(x_numeric) < 2) return(0)
    var(x_numeric, na.rm = TRUE)
  })
}, error = function(e) {
  rep(1, nrow(cpg_clean))  # Fallback: keep all features
})
cpg_filtered <- cpg_clean[cpg_vars > 0.001 & !is.na(cpg_vars), ]

# Filter miRNA: Remove very low expressed miRNAs (mean < 1 TPM)
rownames(mirna_common) <- make.names(rownames(mirna_common), unique = TRUE)
mirna_clean <- mirna_common[complete.cases(mirna_common), ]
mirna_means <- tryCatch({
  apply(mirna_clean, 1, function(x) {
    x_numeric <- as.numeric(x)
    if (all(is.na(x_numeric))) return(0)
    mean(x_numeric, na.rm = TRUE)
  })
}, error = function(e) {
  rowMeans(mirna_clean, na.rm = TRUE)
})
mirna_filtered <- mirna_clean[mirna_means > 1, ]

# Filter mRNA: Remove very low expressed genes (mean < 10 TPM)
rownames(mrna_common) <- make.names(rownames(mrna_common), unique = TRUE)
mrna_clean <- mrna_common[complete.cases(mrna_common), ]
mrna_means <- tryCatch({
  apply(mrna_clean, 1, function(x) {
    x_numeric <- as.numeric(x)
    if (all(is.na(x_numeric))) return(0)
    mean(x_numeric, na.rm = TRUE)
  })
}, error = function(e) {
  rowMeans(mrna_clean, na.rm = TRUE)
})
mrna_filtered <- mrna_clean[mrna_means > 10, ]

cat(sprintf("CpG features after filtering: %d, miRNA: %d, mRNA: %d\n", 
            nrow(cpg_filtered), nrow(mirna_filtered), nrow(mrna_filtered)))

# Select top variable features for correlation analysis
select_top_variable <- function(data, n_features = 500) {
  if (nrow(data) <= n_features) return(data)
  
  cat(sprintf("Selecting top %d features from %d total features\n", n_features, nrow(data)))
  
  # Calculate coefficient of variation safely
  feature_cv <- tryCatch({
    apply(data, 1, function(x) {
      x_numeric <- as.numeric(x)
      x_clean <- x_numeric[!is.na(x_numeric) & is.finite(x_numeric)]
      
      if (length(x_clean) < 2) return(0)
      
      mean_x <- mean(x_clean)
      sd_x <- sd(x_clean)
      
      if (mean_x == 0 || is.na(mean_x) || is.na(sd_x) || !is.finite(mean_x) || !is.finite(sd_x)) {
        return(0)
      }
      
      cv <- sd_x / mean_x
      if (!is.finite(cv)) return(0)
      return(cv)
    })
  }, error = function(e) {
    apply(data, 1, function(x) {
      x_numeric <- as.numeric(x)
      x_clean <- x_numeric[!is.na(x_numeric) & is.finite(x_numeric)]
      if (length(x_clean) < 2) return(0)
      var(x_clean)
    })
  })
  
  # Remove features with zero, NA, or infinite CV
  valid_features <- !is.na(feature_cv) & is.finite(feature_cv) & feature_cv > 0
  
  if (sum(valid_features) == 0) {
    return(data[1:min(n_features, nrow(data)), ])
  }
  
  feature_cv_valid <- feature_cv[valid_features]
  data_valid <- data[valid_features, ]
  
  # Select top variable features
  if (length(feature_cv_valid) > n_features) {
    top_indices <- order(feature_cv_valid, decreasing = TRUE)[1:n_features]
    return(data_valid[top_indices, ])
  } else {
    return(data_valid)
  }
}

# Select top variable features from each omics layer
cat("\nSelecting top variable features for analysis...\n")
cpg_subset <- select_top_variable(cpg_filtered, 300)
mirna_subset <- select_top_variable(mirna_filtered, 200)
mrna_subset <- select_top_variable(mrna_filtered, 300)

cat(sprintf("Final feature counts - CpG: %d, miRNA: %d, mRNA: %d\n", 
            nrow(cpg_subset), nrow(mirna_subset), nrow(mrna_subset)))

# Standardize the data for correlation analysis
standardize_data <- function(data) {
  cat(sprintf("Standardizing %d features...\n", nrow(data)))
  
  standardized <- tryCatch({
    t(apply(data, 1, function(x) {
      x_numeric <- as.numeric(x)
      x_clean <- x_numeric[!is.na(x_numeric) & is.finite(x_numeric)]
      
      if (length(x_clean) < 2) return(x_numeric)
      
      mean_x <- mean(x_clean)
      sd_x <- sd(x_clean)
      
      if (sd_x == 0 || is.na(sd_x) || !is.finite(sd_x)) {
        return(x_numeric)
      }
      
      standardized_values <- (x_numeric - mean_x) / sd_x
      standardized_values[!is.finite(standardized_values)] <- 0
      return(standardized_values)
    }))
  }, error = function(e) {
    as.matrix(data)
  })
  
  return(standardized)
}

cpg_std <- standardize_data(cpg_subset)
mirna_std <- standardize_data(mirna_subset)
mrna_std <- standardize_data(mrna_subset)

# Create combined dataframe for correlation analysis
rownames(cpg_std) <- paste0("CpG_", rownames(cpg_std))
rownames(mirna_std) <- paste0("miRNA_", rownames(mirna_std))
rownames(mrna_std) <- paste0("mRNA_", rownames(mrna_std))

# Combine all omics data
omics_combined <- rbind(cpg_std, mirna_std, mrna_std)

# Transpose for correlation (features as columns)
omics_matrix <- t(omics_combined)

# Create metadata for annotations
meta_df <- data.frame(
  ftr_name = colnames(omics_matrix),
  omic_layer = case_when(
    str_starts(colnames(omics_matrix), "CpG_") ~ "CpG",
    str_starts(colnames(omics_matrix), "miRNA_") ~ "miRNA", 
    str_starts(colnames(omics_matrix), "mRNA_") ~ "mRNA"
  )
)

# Calculate correlation matrix
cat("Calculating correlation matrix...\n")
cormat <- cor(omics_matrix, method = "pearson", use = "pairwise.complete.obs")

# Create annotations for heatmap
annotation <- data.frame(
  ftr_name = colnames(cormat),
  index = 1:ncol(cormat)
) %>%
  left_join(meta_df, by = "ftr_name") %>%
  mutate(omic_layer = str_to_title(omic_layer))

# ================== FIGURE 1: Generate correlation heatmap ==================
cat("Generating correlation heatmap from true raw data...\n")
jpeg(file.path(output_dir, "correlation_heatmap_true_raw.jpg"), width = 12, height = 10, units = "in", res = 300, quality = 95)
ht <- Heatmap(cormat,
              row_split = annotation$omic_layer,
              column_split = annotation$omic_layer,
              show_row_names = FALSE,
              show_column_names = FALSE,
              column_title_gp = gpar(fontsize = 12),
              row_title_gp = gpar(fontsize = 12),
              heatmap_legend_param = list(title = "Correlation"),
              col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
              name = "Correlation")
draw(ht)
dev.off()

# ================== PART 2: Load Adjusted Data for Sample Distribution ==================

cat("\nLoading adjusted data files for sample distribution...\n")

# Read adjusted CpG data (first column is empty, row names in first column after comma)
adj_cpg_data <- read.csv(adjusted_cpg_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
adj_cpg_data <- data.frame(lapply(adj_cpg_data, function(x) as.numeric(as.character(x))))
rownames(adj_cpg_data) <- gsub("CpG_", "", rownames(adj_cpg_data))

# Read adjusted miRNA data (sRNA column contains row names)
adj_mirna_raw <- read.csv(adjusted_mirna_file, check.names = FALSE, stringsAsFactors = FALSE)
adj_mirna_data <- adj_mirna_raw[, -1]  # Remove sRNA column
rownames(adj_mirna_data) <- adj_mirna_raw$sRNA
adj_mirna_data <- data.frame(lapply(adj_mirna_data, function(x) as.numeric(as.character(x))))
rownames(adj_mirna_data) <- gsub("hsa-", "", rownames(adj_mirna_data))

# Read adjusted mRNA data (first column is empty, row names in first column after comma)
adj_mrna_data <- read.csv(adjusted_mrna_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
adj_mrna_data <- data.frame(lapply(adj_mrna_data, function(x) as.numeric(as.character(x))))

cat(sprintf("Adjusted data - CpG: %d x %d, miRNA: %d x %d, mRNA: %d x %d\n",
            nrow(adj_cpg_data), ncol(adj_cpg_data),
            nrow(adj_mirna_data), ncol(adj_mirna_data),
            nrow(adj_mrna_data), ncol(adj_mrna_data)))

# Check for common samples in adjusted data
adj_common_samples <- Reduce(intersect, list(colnames(adj_cpg_data), 
                                             colnames(adj_mirna_data), 
                                             colnames(adj_mrna_data)))
cat(sprintf("Common samples in adjusted data: %d\n", length(adj_common_samples)))

# Subset to common samples
adj_cpg_common <- adj_cpg_data[, adj_common_samples]
adj_mirna_common <- adj_mirna_data[, adj_common_samples]
adj_mrna_common <- adj_mrna_data[, adj_common_samples]

# ================== FIGURE 2: Generate sample distribution and boxplot ==================
cat("Generating sample distribution and boxplot figure...\n")
jpeg(file.path(output_dir, "sample_distribution_boxplot.jpg"), width = 16, height = 10, units = "in", res = 300, quality = 95)

# Set up 2x3 layout
par(mfrow = c(2, 3))

# Top row: Sample distribution plots (1x3)
plot_distribution <- function(data, title) {
  sample_means <- colMeans(data, na.rm = TRUE)
  hist(sample_means, main = paste("Sample Distribution -", title),
       xlab = "Mean Expression", col = "lightblue", breaks = 10,
       cex.main = 1.2, cex.lab = 1.1)
}

plot_distribution(adj_cpg_common, "CpG")
plot_distribution(adj_mirna_common, "miRNA") 
plot_distribution(adj_mrna_common, "mRNA")

# Bottom row: Boxplot for expression distribution (1x3)
plot_boxplot <- function(data, title) {
  # Convert data to long format for boxplot
  plot_data <- data %>%
    as.data.frame() %>%
    rownames_to_column("Feature") %>%
    pivot_longer(cols = -Feature, names_to = "Sample", values_to = "Expression") %>%
    filter(!is.na(Expression))
  
  boxplot(plot_data$Expression ~ plot_data$Sample, 
          main = paste("Expression Distribution -", title),
          xlab = "Samples", ylab = "Expression Level",
          col = "lightblue", las = 2, cex.axis = 0.7,
          cex.main = 1.2, cex.lab = 1.1)
}

plot_boxplot(adj_cpg_common, "CpG")
plot_boxplot(adj_mirna_common, "miRNA")
plot_boxplot(adj_mrna_common, "mRNA")

dev.off()

cat("Analysis complete! Generated 2 figures:\n")
cat("1. correlation_heatmap_true_raw.jpg - Multi-omics correlation heatmap (from raw data)\n")
cat("2. sample_distribution_boxplot.jpg - Sample distributions and boxplots (from adjusted data)\n")

# Print session info
cat("\nSession Info:\n")
print(sessionInfo())
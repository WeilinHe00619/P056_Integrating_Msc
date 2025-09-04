
rm(list=ls())

# Load required packages
library(ggplot2)
library(tibble)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(patchwork)
library(RColorBrewer)

# Check if RColorBrewer is available, if not use base colors
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  print("RColorBrewer not available, using base R colors")
}

##### Configuration
INPUT_DIR <- "/Users/heweilin/Desktop/P056"
OUTPUT_DIR <- "/Users/heweilin/Desktop/P056_Code_3/Figure"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

##### 01. Load Data
print("Loading differential expression data...")

# Load DEGs (mRNA)
degs_df <- read.csv(file.path(INPUT_DIR, "1mRNA_DEGs_proteincoding.csv"))
cat("Loaded", nrow(degs_df), "mRNA records\n")

# Load DEmiRs
demirs_df <- read.csv(file.path(INPUT_DIR, "2miRNA_DEmirs.csv"))
cat("Loaded", nrow(demirs_df), "miRNA records\n")

# Load DMRs
dmrs_df <- read.csv(file.path(INPUT_DIR, "4DNA_DMRs.csv"))
cat("Loaded", nrow(dmrs_df), "DMR records\n")

##### 02. Data Processing and Classification
print("Classifying differential features...")

# Classify DEGs
degs_df$significance <- "normal"
degs_df$significance[degs_df$log2FoldChange > 0 & degs_df$padj < 0.05] <- "up"
degs_df$significance[degs_df$log2FoldChange < 0 & degs_df$padj < 0.05] <- "down"
degs_df$neg_log10_padj <- -log10(degs_df$padj)

# Replace infinite values and handle missing data
degs_df$neg_log10_padj[is.infinite(degs_df$neg_log10_padj) | is.na(degs_df$neg_log10_padj)] <- 
  max(degs_df$neg_log10_padj[!is.infinite(degs_df$neg_log10_padj) & !is.na(degs_df$neg_log10_padj)], na.rm = TRUE) + 1

# Count DEGs
deg_counts <- table(degs_df$significance)
cat("DEGs:", sum(deg_counts[c("up", "down")]), "significant out of", nrow(degs_df), "total\n")
cat("  - Upregulated:", deg_counts["up"], "\n")
cat("  - Downregulated:", deg_counts["down"], "\n")

# Classify DEmiRs
demirs_df$significance <- "normal"
demirs_df$significance[demirs_df$log2FoldChange > 0 & demirs_df$pvalue < 0.05] <- "up"
demirs_df$significance[demirs_df$log2FoldChange < 0 & demirs_df$pvalue < 0.05] <- "down"
demirs_df$neg_log10_pvalue <- -log10(demirs_df$pvalue)

# Replace infinite values and handle missing data
demirs_df$neg_log10_pvalue[is.infinite(demirs_df$neg_log10_pvalue) | is.na(demirs_df$neg_log10_pvalue)] <- 
  max(demirs_df$neg_log10_pvalue[!is.infinite(demirs_df$neg_log10_pvalue) & !is.na(demirs_df$neg_log10_pvalue)], na.rm = TRUE) + 1

# Count DEmiRs
demir_counts <- table(demirs_df$significance)
cat("DEmiRs:", sum(demir_counts[c("up", "down")]), "significant out of", nrow(demirs_df), "total\n")
cat("  - Upregulated:", demir_counts["up"], "\n")
cat("  - Downregulated:", demir_counts["down"], "\n")

# DMR counts
dmr_counts <- table(dmrs_df$DiffMethylated)
cat("DMRs:", nrow(dmrs_df), "significant regions\n")
cat("  - Hypermethylated:", dmr_counts["HYPER"], "\n")
cat("  - Hypomethylated:", dmr_counts["HYPO"], "\n")

##### 03. Enhanced Volcano Plots (Style similar to reference)
print("Creating enhanced volcano plots...")

# Color schemes - using the new provided color palette (excluding middle white)
# RGB to HEX conversion for the new palette
new_palette <- c(
  "#1B3B70",  # R:027, G:059, B:112 (最蓝)
  "#276FB3",  # R:039, G:111, B:175 (第二蓝)
  "#4D9AC7",  # R:077, G:154, B:199 (第三蓝)
  "#99C8E0",  # R:153, G:200, B:224 (第四蓝)
  "#D4E6EF",  # R:212, G:230, B:239 (第五蓝)
  # 跳过中间白色 R:248, G:244, B:242
  "#F6D8C3",  # R:251, G:216, B:195 (第一橙)
  "#F2A281",  # R:242, G:164, B:129 (第二橙)
  "#D6564D",  # R:214, G:096, B:077 (第三红)
  "#B5202E",  # R:181, G:032, B:046 (第二红)
  "#70221C"   # R:112, G:012, B:034 (最红)
)

# For volcano plots and two-color schemes (second blue and second red)
volcano_colors <- c("up" = "#B5202E", "down" = "#276FB3", "normal" = "#474747")

# For methylation pie chart (second blue and second red)
meth_colors <- c("HYPER" = "#B5202E", "HYPO" = "#276FB3")

# DEGs volcano plot data preparation
data_bg_degs <- degs_df[degs_df$significance == "normal", ]
data_sig_degs <- degs_df[degs_df$significance != "normal", ]

# Select top genes for labeling (most significant by p-value and fold change)
top_up_degs <- degs_df[degs_df$significance == "up" & !is.na(degs_df$log2FoldChange) & !is.na(degs_df$SYMBOL), ] %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  head(5)

top_down_degs <- degs_df[degs_df$significance == "down" & !is.na(degs_df$log2FoldChange) & !is.na(degs_df$SYMBOL), ] %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  head(5)

data_label_degs <- rbind(top_up_degs, top_down_degs)

# DEGs Volcano Plot
p1 <- ggplot(degs_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(data = data_bg_degs, shape = 21, fill = "#474747", color = "black", 
             alpha = 0.6, size = 1.2, stroke = 0.2) +
  geom_point(data = data_sig_degs, aes(fill = significance), shape = 21, 
             color = "black", size = 1.8, stroke = 0.3, alpha = 0.8) +
  scale_fill_manual(values = volcano_colors,
                    labels = c("up" = paste0("Upregulated (", deg_counts["up"], ")"),
                               "down" = paste0("Downregulated (", deg_counts["down"], ")")),
                    name = "DEGs") +
  geom_text_repel(data = data_label_degs, 
                  aes(x = log2FoldChange, y = neg_log10_padj, label = SYMBOL),
                  size = 3, force_pull = 0, point.padding = 0, box.padding = 0.8,
                  min.segment.length = 0, vjust = 0.5, segment.color = "grey20",
                  segment.size = 0.4, segment.alpha = 0.8, max.overlaps = Inf,
                  seed = 123, nudge_y = 0.3) +
  geom_hline(yintercept = -log10(0.05), color = "#474747", 
             linewidth = 0.6, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), color = "#474747", 
             linewidth = 0.4, linetype = "dotted", alpha = 0.7) +
  scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, 2)) +
  scale_y_continuous(limits = c(0, 5), expand = c(0.02, 0)) +
  labs(title = "Differentially Expressed Genes", 
       x = expression(log[2](Fold~Change)), 
       y = expression(-log[10](adj.P.value))) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(color = "grey90", size = 0.3),
        legend.box.background = element_rect(color = "black", size = 0.5))

# DEmiRs volcano plot data preparation
data_bg_demirs <- demirs_df[demirs_df$significance == "normal", ]
data_sig_demirs <- demirs_df[demirs_df$significance != "normal", ]

# Create miRNA identifier column
if ("Unnamed..0" %in% colnames(demirs_df)) {
  demirs_df$miRNA_name <- demirs_df$Unnamed..0
} else if ("X" %in% colnames(demirs_df)) {
  demirs_df$miRNA_name <- demirs_df$X
} else {
  demirs_df$miRNA_name <- rownames(demirs_df)
}

# Select top miRNAs for labeling (most significant by p-value and fold change)
top_up_demirs <- demirs_df[demirs_df$significance == "up" & !is.na(demirs_df$log2FoldChange), ] %>%
  arrange(pvalue, desc(abs(log2FoldChange))) %>%
  head(3)

top_down_demirs <- demirs_df[demirs_df$significance == "down" & !is.na(demirs_df$log2FoldChange), ] %>%
  arrange(pvalue, desc(abs(log2FoldChange))) %>%
  head(3)

data_label_demirs <- rbind(top_up_demirs, top_down_demirs)

# DEmiRs Volcano Plot
p2 <- ggplot(demirs_df, aes(x = log2FoldChange, y = neg_log10_pvalue)) +
  geom_point(data = data_bg_demirs, shape = 21, fill = "#474747", color = "black", 
             alpha = 0.6, size = 1.2, stroke = 0.2) +
  geom_point(data = data_sig_demirs, aes(fill = significance), shape = 21, 
             color = "black", size = 1.8, stroke = 0.3, alpha = 0.8) +
  scale_fill_manual(values = volcano_colors,
                    labels = c("up" = paste0("Upregulated (", demir_counts["up"], ")"),
                               "down" = paste0("Downregulated (", demir_counts["down"], ")")),
                    name = "DEmiRs") +
  geom_text_repel(data = data_label_demirs, 
                  aes(x = log2FoldChange, y = neg_log10_pvalue, label = miRNA_name),
                  size = 3, force_pull = 0, point.padding = 0, box.padding = 0.8,
                  min.segment.length = 0, vjust = 0.5, segment.color = "grey20",
                  segment.size = 0.4, segment.alpha = 0.8, max.overlaps = Inf,
                  seed = 123, nudge_y = 0.3) +
  geom_hline(yintercept = -log10(0.05), color = "#474747", 
             linewidth = 0.6, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), color = "#474747", 
             linewidth = 0.4, linetype = "dotted", alpha = 0.7) +
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-4, 4, 2)) +
  scale_y_continuous(limits = c(0, 5.5), expand = c(0.02, 0)) +
  labs(title = "Differentially Expressed miRNAs", 
       x = expression(log[2](Fold~Change)), 
       y = expression(-log[10](P.value))) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(color = "grey90", size = 0.3),
        legend.box.background = element_rect(color = "black", size = 0.5))

# Combine volcano plots
volcano_combined <- p1 + p2 + plot_layout(ncol = 2)

ggsave(filename = file.path(OUTPUT_DIR, "01_volcano_plots_enhanced.png"), 
       plot = volcano_combined, width = 14, height = 6, dpi = 300, bg = "white")

##### 04. Fold Change Analysis (1×3 layout with statistics)
print("Creating fold change analysis...")

# Get significant data
degs_significant <- degs_df[degs_df$significance != "normal", ]
demirs_significant <- demirs_df[demirs_df$significance != "normal", ]

# Calculate statistics
deg_stats <- degs_significant %>%
  summarise(
    Mean = round(mean(log2FoldChange, na.rm = TRUE), 3),
    Median = round(median(log2FoldChange, na.rm = TRUE), 3),
    SD = round(sd(log2FoldChange, na.rm = TRUE), 3),
    Min = round(min(log2FoldChange, na.rm = TRUE), 3),
    Max = round(max(log2FoldChange, na.rm = TRUE), 3)
  )

demir_stats <- demirs_significant %>%
  summarise(
    Mean = round(mean(log2FoldChange, na.rm = TRUE), 3),
    Median = round(median(log2FoldChange, na.rm = TRUE), 3),
    SD = round(sd(log2FoldChange, na.rm = TRUE), 3),
    Min = round(min(log2FoldChange, na.rm = TRUE), 3),
    Max = round(max(log2FoldChange, na.rm = TRUE), 3)
  )

# DEGs fold change histogram with statistics
stats_text_deg <- paste0("DEGs Statistics:\n",
                         "Mean: ", deg_stats$Mean, "\n",
                         "Median: ", deg_stats$Median, "\n",
                         "SD: ", deg_stats$SD, "\n",
                         "Range: ", deg_stats$Min, " to ", deg_stats$Max, "\n",
                         "n = ", nrow(degs_significant))

p3 <- ggplot(degs_significant, aes(x = log2FoldChange)) +
  geom_histogram(bins = 25, fill = "#276FB3", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = deg_stats$Mean, color = "#B5202E", linetype = "dashed", size = 1) +
  annotate("text", x = min(degs_significant$log2FoldChange, na.rm = TRUE) + 0.5, 
           y = max(hist(degs_significant$log2FoldChange, plot = FALSE)$counts) * 0.8,
           label = stats_text_deg, hjust = 0, vjust = 1, size = 3,
           fontfamily = "mono", color = "black") +
  labs(title = "DEGs Fold Change Distribution",
       x = expression(log[2](Fold~Change)), y = "Frequency") +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10))

# DEmiRs fold change histogram with statistics
stats_text_demir <- paste0("DEmiRs Statistics:\n",
                           "Mean: ", demir_stats$Mean, "\n",
                           "Median: ", demir_stats$Median, "\n",
                           "SD: ", demir_stats$SD, "\n",
                           "Range: ", demir_stats$Min, " to ", demir_stats$Max, "\n",
                           "n = ", nrow(demirs_significant))

p4 <- ggplot(demirs_significant, aes(x = log2FoldChange)) +
  geom_histogram(bins = 15, fill = "#B5202E", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = demir_stats$Mean, color = "#276FB3", linetype = "dashed", size = 1) +
  annotate("text", x = min(demirs_significant$log2FoldChange, na.rm = TRUE) + 0.5, 
           y = max(hist(demirs_significant$log2FoldChange, plot = FALSE)$counts) * 0.8,
           label = stats_text_demir, hjust = 0, vjust = 1, size = 3,
           fontfamily = "mono", color = "black") +
  labs(title = "DEmiRs Fold Change Distribution",
       x = expression(log[2](Fold~Change)), y = "Frequency") +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10))

# Box plot comparison
data_combined <- rbind(
  data.frame(Type = "DEGs", log2FoldChange = degs_significant$log2FoldChange),
  data.frame(Type = "DEmiRs", log2FoldChange = demirs_significant$log2FoldChange)
)

p5 <- ggplot(data_combined, aes(x = Type, y = log2FoldChange, fill = Type)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  scale_fill_manual(values = c("DEGs" = "#276FB3", "DEmiRs" = "#B5202E")) +
  labs(title = "Fold Change Comparison",
       x = "Feature Type", y = expression(log[2](Fold~Change))) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10))

# Combine fold change plots
fc_combined <- p3 + p4 + p5 + plot_layout(ncol = 3)

ggsave(filename = file.path(OUTPUT_DIR, "02_fold_change_analysis.png"), 
       plot = fc_combined, width = 18, height = 5, dpi = 300, bg = "white")

# DMR Genomic Region Distribution (horizontal bar chart with wider layout)
region_counts <- dmrs_df %>%
  count(annot.type) %>%
  mutate(percentage = round(n/sum(n)*100, 1))

# Generate colors using the new palette (excluding middle white)
n_regions <- nrow(region_counts)
if(n_regions <= length(new_palette)) {
  region_colors <- new_palette[1:n_regions]
} else {
  # If more regions than colors, extend with gradient
  region_colors <- c(new_palette, 
                     colorRampPalette(c("#1B3B70", "#70221C"))(n_regions - length(new_palette)))
}
names(region_colors) <- region_counts$annot.type

p6 <- ggplot(region_counts, aes(x = reorder(annot.type, n), y = n, fill = annot.type)) +
  geom_col(alpha = 0.8, color = "white", size = 0.5) +
  geom_text(aes(label = paste0(n, " (", percentage, "%)")), 
            hjust = -0.05, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = region_colors) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "DMR Genomic Region Distribution",
       x = "Genomic Region", y = "Count") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 20, 5, 5),
        panel.grid.major.x = element_line(color = "grey90", size = 0.3))

# DMR Methylation Status Pie Chart
meth_counts <- dmrs_df %>%
  count(DiffMethylated) %>%
  mutate(percentage = round(n/sum(n)*100, 1),
         label = paste0(DiffMethylated, "\n", n, " (", percentage, "%)"))

meth_colors <- c("HYPER" = "#B5202E", "HYPO" = "#276FB3")

p7 <- ggplot(meth_counts, aes(x = "", y = n, fill = DiffMethylated)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = meth_colors, labels = meth_counts$label) +
  labs(title = "DMR Methylation Status", fill = "Status") +
  theme_void() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10))

# Combine DMR plots (2×1 layout instead of 1×3)
dmr_combined <- p6 + p7 + plot_layout(ncol = 2, widths = c(2, 1))

ggsave(filename = file.path(OUTPUT_DIR, "03_dmr_analysis.png"), 
       plot = dmr_combined, width = 14, height = 6, dpi = 300, bg = "white")

##### 06. Chromosomal Distribution Analysis (2×1 layout)
print("Creating chromosomal distribution plots...")

# DEGs chromosomal distribution
degs_sig <- degs_df[degs_df$significance != "normal", ]
deg_chr_counts <- degs_sig %>%
  count(Chromosome, name = "DEG_count") %>%
  arrange(as.numeric(ifelse(Chromosome %in% c("X", "Y"), 
                            ifelse(Chromosome == "X", 23, 24), 
                            Chromosome)))

# Sort chromosomes properly
chr_order <- c(as.character(1:22), "X", "Y")
deg_chr_counts <- deg_chr_counts %>%
  filter(Chromosome %in% chr_order) %>%
  mutate(Chromosome = factor(Chromosome, levels = chr_order)) %>%
  arrange(Chromosome)

# DMRs chromosomal distribution
dmrs_df$chromosome <- gsub("chr", "", dmrs_df$seqnames)
dmr_chr_counts <- dmrs_df %>%
  count(chromosome, name = "DMR_count") %>%
  filter(chromosome %in% chr_order) %>%
  mutate(chromosome = factor(chromosome, levels = chr_order)) %>%
  arrange(chromosome)

# Create 2×1 plot
fig_chr <- plot_grid(
  # DEGs plot
  ggplot(deg_chr_counts, aes(x = Chromosome, y = DEG_count)) +
    geom_col(fill = "#276FB3", alpha = 0.8, color = "white", size = 0.3) +
    geom_text(aes(label = DEG_count), vjust = -0.3, size = 3, angle = 0) +
    labs(title = paste0("Chromosomal Distribution of DEGs (n=", sum(deg_chr_counts$DEG_count), ")"),
         x = "Chromosome", y = "Number of DEGs") +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 11),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.y = element_line(color = "grey90", size = 0.3)),
  
  # DMRs plot
  ggplot(dmr_chr_counts, aes(x = chromosome, y = DMR_count)) +
    geom_col(fill = "#B5202E", alpha = 0.8, color = "white", size = 0.3) +
    geom_text(aes(label = DMR_count), vjust = -0.3, size = 3, angle = 0) +
    labs(title = paste0("Chromosomal Distribution of DMRs (n=", sum(dmr_chr_counts$DMR_count), ")"),
         x = "Chromosome", y = "Number of DMRs") +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 11),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.y = element_line(color = "grey90", size = 0.3)),
  
  ncol = 1
)

ggsave(filename = file.path(OUTPUT_DIR, "04_chromosomal_distribution.png"), 
       plot = fig_chr, width = 14, height = 10, dpi = 300, bg = "white")

suppressMessages({
  library(gground)
  library(ggprism)
  library(tidyverse)
  library(org.Hs.eg.db)
  library(clusterProfiler)
})

INPUT_PATH <- "/Users/heweilin/Desktop/P056_Code_2/Processed_Data"
OUTPUT_PATH <- "/Users/heweilin/Desktop/P056_Code_2/Figure"

dir.create(OUTPUT_PATH, recursive = TRUE, showWarnings = FALSE)

pathway_file <- file.path(INPUT_PATH, "2_Enrichment_use_pathway.csv")

if (!file.exists(pathway_file)) {
  stop("File not found: ", pathway_file)
}

raw_data <- read.csv(pathway_file, stringsAsFactors = FALSE)

use_pathway <- raw_data %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c('BP', 'CC', 'MF', 'KEGG')))) %>%
  dplyr::arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  select(-any_of("index")) %>%
  tibble::rowid_to_column('index')

pal <- c('#eaa052', '#b74147', '#90ad5b', '#23929c')
names(pal) <- c('BP', 'CC', 'MF', 'KEGG')

width <- 0.5
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1

rect.data <- use_pathway %>%
  group_by(ONTOLOGY) %>%
  summarize(n = n(), .groups = 'drop') %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )

plot_enrichment <- function() {
  p <- use_pathway %>%
    ggplot(aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
    geom_round_col(
      aes(y = Description), 
      width = 0.6, 
      alpha = 0.8
    ) +
    geom_text(
      aes(x = 0.05, label = Description),
      hjust = 0, 
      size = 4,
      color = "black"
    ) +
    geom_text(
      aes(x = 0.1, label = geneID, colour = ONTOLOGY),
      hjust = 0, 
      vjust = 2.6, 
      size = 3, 
      fontface = 'italic',
      show.legend = FALSE
    ) +
    geom_point(
      aes(x = -width, size = Count),
      shape = 21,
      color = "black",
      fill = "white",
      stroke = 1
    ) +
    geom_text(
      aes(x = -width, label = Count),
      size = 3,
      color = "black"
    ) +
    scale_size_continuous(
      name = 'Gene Count', 
      range = c(4, 10),
      breaks = c(5, 10, 15, 20)
    ) +
    geom_round_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ONTOLOGY),
      data = rect.data,
      radius = unit(2, 'mm'),
      inherit.aes = FALSE,
      alpha = 0.9
    ) +
    geom_text(
      aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
      data = rect.data,
      inherit.aes = FALSE,
      color = "white",
      fontface = "bold",
      size = 4
    ) +
    annotate(
      "segment",
      x = 0, y = 0, 
      xend = xaxis_max, yend = 0,
      linewidth = 1.5,
      color = "black"
    ) +
    labs(
      x = expression(-log[10](p.adjust)),
      y = NULL,
      title = "GO/KEGG Enrichment Analysis",
      subtitle = paste("Top enriched pathways (n =", nrow(use_pathway), ")")
    ) +
    scale_fill_manual(
      name = 'Category', 
      values = pal
    ) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(
      breaks = seq(0, ceiling(xaxis_max), 1),
      expand = expansion(c(0, 0.02))
    ) +
    theme_prism() +
    theme(
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title.x = element_text(size = 12, face = "bold"),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      legend.position = "right",
      plot.margin = margin(20, 20, 20, 20)
    )
  
  return(p)
}

tryCatch({
  enrichment_plot <- plot_enrichment()
}, error = function(e) {
  plot_enrichment_simple <- function() {
    p <- use_pathway %>%
      ggplot(aes(x = -log10(p.adjust), y = reorder(Description, -index), fill = ONTOLOGY)) +
      geom_col(width = 0.6, alpha = 0.8) +
      geom_text(
        aes(label = Count, x = -0.1),
        hjust = 1,
        size = 3,
        color = "black"
      ) +
      labs(
        x = expression(-log[10](p.adjust)),
        y = NULL,
        title = "GO/KEGG Enrichment Analysis",
        subtitle = paste("Top enriched pathways (n =", nrow(use_pathway), ")")
      ) +
      scale_fill_manual(
        name = 'Category', 
        values = pal
      ) +
      theme_prism() +
      theme(
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    
    return(p)
  }
  
  enrichment_plot <- plot_enrichment_simple()
})

output_file <- file.path(OUTPUT_PATH, "GO_KEGG_Enrichment.jpg")

ggsave(
  filename = output_file, 
  plot = enrichment_plot, 
  width = 14, 
  height = 10, 
  dpi = 300, 
  units = "in",
  bg = "white"
)

print(enrichment_plot)
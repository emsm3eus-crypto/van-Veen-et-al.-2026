need <- c("tidyverse")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, Ncpus = 2)
invisible(lapply(need, library, character.only = TRUE))


setwd("C:/Users/u0171728/Downloads")

df <- read.csv(
  "fe-wp-dataset-124 (2).csv",
  sep = ",",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Data dimensions: ", nrow(df), " rows x ", ncol(df), " columns\n", sep = "")
cat("First few column names:\n")
print(colnames(df)[1:min(10, ncol(df))])


rownames(df) <- make.unique(as.character(df[[1]]))

df_expr <- df[, -1, drop = FALSE]


df_expr[] <- lapply(df_expr, as.numeric)


aggregate_celltypes <- function(prefix, df) {
  cols <- grep(paste0("^", prefix), colnames(df), value = TRUE)
  
  if (length(cols) == 0) return(NULL)
  
  if (length(cols) == 1) {
    return(df[[cols]])
  } else {
    return(rowMeans(df[, cols, drop = FALSE], na.rm = TRUE))
  }
}


cell_types <- c(
  "astrocytes_fetal",
  "astrocytes_mature",
  "neurons",
  "oligodendrocytes",
  "microglla",
  "endothelial"
)


agg_list <- lapply(cell_types, aggregate_celltypes, df = df_expr)
valid_idx <- !sapply(agg_list, is.null)

agg_list <- agg_list[valid_idx]
valid_cell_types <- cell_types[valid_idx]

agg_df <- as.data.frame(agg_list)
colnames(agg_df) <- valid_cell_types
rownames(agg_df) <- rownames(df)

cat("\nDetected cell-type prefixes:\n")
print(valid_cell_types)


nice_labels <- c(
  astrocytes_fetal = "Fetal astrocytes",
  astrocytes_mature = "Mature astrocytes",
  neurons = "Neurons",
  oligodendrocytes = "Oligodendrocytes",
  microglla = "Microglia/macrophage",
  endothelial = "Endothelial"
)

color_map <- c(
  astrocytes_fetal = "red",
  astrocytes_mature = "red",
  neurons = "blue",
  oligodendrocytes = "blue",
  microglla = "blue",
  endothelial = "red"
)


genes_to_plot <- c("ATP13A2", "ATP13A3", "ATP13A4", "ATP13A5")


dir.create("zhang2016_atp13_plots", showWarnings = FALSE)

plot_gene <- function(gene_symbol, agg_df, valid_cell_types, nice_labels, color_map) {
  
  # Search gene symbol inside the long row label
  gene_hits <- grep(gene_symbol, rownames(agg_df), ignore.case = TRUE, value = TRUE)
  
  if (length(gene_hits) == 0) {
    cat("Gene not found: ", gene_symbol, "\n", sep = "")
    return(NULL)
  }
  
  gene_name <- gene_hits[1]
  
  cat("\nMatched row for ", gene_symbol, ":\n", sep = "")
  print(gene_name)
  
  gene_expr <- agg_df[gene_name, , drop = FALSE]
  gene_expr <- as.data.frame(t(gene_expr))
  gene_expr$CellType <- rownames(gene_expr)
  colnames(gene_expr)[1] <- "Expression"
  
  gene_expr <- gene_expr %>%
    mutate(
      Expression = as.numeric(Expression),
      CellType = factor(CellType, levels = rev(valid_cell_types)),
      CellLabel = nice_labels[as.character(CellType)]
    )
  
  p <- ggplot(gene_expr, aes(x = CellType, y = Expression, fill = CellType)) +
    geom_col(show.legend = FALSE, width = 0.7, color = "black", linewidth = 0.3) +
    coord_flip() +
    scale_fill_manual(values = color_map) +
    scale_x_discrete(labels = rev(nice_labels[valid_cell_types])) +
    theme_classic() +
    labs(
      title = paste0("h", gene_symbol, " mRNA"),
      subtitle = "Zhang et al., 2016",
      x = NULL,
      y = "FPKM"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank()
    )
  
  print(p)
  
  ggsave(
    filename = file.path("zhang2016_atp13_plots", paste0(gene_symbol, "_Zhang2016.png")),
    plot = p,
    width = 6,
    height = 4.5,
    dpi = 300
  )
  
  cat("\n====================\n")
  cat(gene_symbol, "\n")
  cat("====================\n")
  print(gene_expr[, c("CellLabel", "Expression")])
  
  invisible(gene_expr)
}


results <- lapply(
  genes_to_plot,
  plot_gene,
  agg_df = agg_df,
  valid_cell_types = valid_cell_types,
  nice_labels = nice_labels,
  color_map = color_map
)

names(results) <- genes_to_plot

cat("\nDone. Plots saved in folder: zhang2016_atp13_plots\n")

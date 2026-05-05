library(tidyverse)
library(pheatmap)

setwd("C:/Users/u0171728/Downloads")
df <- read.csv("fe-wp-dataset-120.csv", check.names = FALSE)
rownames(df) <- make.unique(df$gene_id)
df_expr <- df %>% select(-gene_id, -id)
aggregate_celltypes <- function(prefix, df) {
  cols <- grep(paste0("^", prefix), colnames(df), value = TRUE)
  if (length(cols) == 0) return(NULL)
  if (length(cols) == 1) {
    return(df[[cols]])  # single column vector
  } else {
    return(rowMeans(df[, cols]))  # average of multiple columns
  }
}

cell_types <- c(
  "astrocytes", 
  "endothelial", 
  "microglia_macrophage", 
  "myelinating_oligodendrocyte", 
  "neurons", 
  "newly_formed_oligodendrocyte", 
  "opc"
)


agg_list <- lapply(cell_types, aggregate_celltypes, df = df_expr)
valid_idx <- !sapply(agg_list, is.null)
agg_list <- agg_list[valid_idx]
agg_df <- as.data.frame(agg_list)
valid_cell_types <- cell_types[valid_idx]
colnames(agg_df) <- valid_cell_types
nlgn1_names <- grep("ATP13A2", rownames(agg_df), ignore.case = TRUE, value = TRUE)
print(nlgn1_names)  # just to confirm
gene_name <- nlgn1_names[1]
gene_expr <- agg_df[gene_name, , drop = FALSE]
gene_expr <- as.data.frame(t(gene_expr))
gene_expr$CellType <- valid_cell_types
colnames(gene_expr)[1] <- "Expression"
gene_expr <- gene_expr[, c("CellType", "Expression")]

gene_expr$CellType <- factor(gene_expr$CellType, levels = valid_cell_types)
color_map <- c(
  "neurons" = "forestgreen",
  "astrocytes" = "dodgerblue",
  "microglia_macrophage" = "firebrick",
  "myelinating_oligodendrocyte" = "darkorange",
  "newly_formed_oligodendrocyte" = "darkorange",
  "opc" = "darkorange",
  "endothelial" = "goldenrod"
)

ggplot(gene_expr, aes(x = CellType, y = Expression, fill = CellType)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = color_map) +
  theme_minimal() +
  labs(
    title = paste("ATP13A2 Expression Across Cell Types"),
    subtitle = "Data extracted from Brainrnaseq.org (Zhang et al.,2014)",
    x = "Cell Type",
    y = "FPKM Values"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12)
  )

cat("CellType\tExpression\n")
apply(gene_expr, 1, function(r) cat(r["CellType"], "\t", r["Expression"], "\n", sep=""))

cat("cell_types = ", paste(shQuote(as.character(gene_expr$CellType)), collapse = ", "), "\n", sep = "")
cat("values     = ", paste(signif(as.numeric(gene_expr$Expression), 6), collapse = ", "), "\n", sep = "")


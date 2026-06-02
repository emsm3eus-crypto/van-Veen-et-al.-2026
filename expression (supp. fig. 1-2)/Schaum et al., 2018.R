
library(data.table)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

setwd("C:/Users/u0171728/Downloads/FACS")

gene_for_both <- "Atp13a5"

gene_non_myeloid <- NULL   
gene_myeloid     <- NULL  

gene_non_myeloid <- if (is.null(gene_non_myeloid)) gene_for_both else gene_non_myeloid
gene_myeloid     <- if (is.null(gene_myeloid))     gene_for_both else gene_myeloid

metadata          <- fread("metadata_FACS.csv") %>% as_tibble()     # optional
annotations       <- fread("annotations_facs.csv") %>% as_tibble()
myeloid_counts    <- fread("Brain_Myeloid-counts.csv")
non_myeloid_counts<- fread("Brain_Non-Myeloid-counts.csv")

get_gene_counts_row <- function(counts_df, gene_name) {
  idx <- which(toupper(counts_df[[1]]) == toupper(gene_name))
  if (length(idx) == 0) {
    stop(sprintf("Gene '%s' not found in the first column (V1). Check spelling/case.", gene_name))
  }
  counts_df[idx, , drop = FALSE]
}

compute_tp100k <- function(counts_df, gene_name) {
  gene_row <- get_gene_counts_row(counts_df, gene_name)
  gene_counts <- gene_row %>% select(-1)
  total_counts <- colSums(counts_df %>% select(-1), na.rm = TRUE)
  tp100k <- (as.numeric(gene_counts) / total_counts) * 100000
  tibble(
    sample  = colnames(counts_df)[-1],
    TP100K  = tp100k,
    gene    = gene_name
  )
}

df_myeloid <- compute_tp100k(myeloid_counts, gene_myeloid) %>%
  left_join(annotations, by = c("sample" = "cell")) %>%
  mutate(cell_group = "Myeloid")

df_non_myeloid <- compute_tp100k(non_myeloid_counts, gene_non_myeloid) %>%
  left_join(annotations, by = c("sample" = "cell")) %>%
  mutate(cell_group = "Non-Myeloid")

combined_df <- bind_rows(df_non_myeloid, df_myeloid) %>%
  filter(!is.na(cell_ontology_class), !is.na(TP100K))

combined_df <- combined_df %>%
  mutate(Category = case_when(
    grepl("neuron",          cell_ontology_class, ignore.case = TRUE) ~ "Neuron",
    grepl("astro",           cell_ontology_class, ignore.case = TRUE) ~ "Astrocytes",
    grepl("oligodendrocyte", cell_ontology_class, ignore.case = TRUE) ~ "Oligodendrocytes",
    grepl("\\bopc\\b",       cell_ontology_class, ignore.case = TRUE) ~ "OPCs",
    grepl("microglia",       cell_ontology_class, ignore.case = TRUE) ~ "Microglia",
    TRUE ~ "Other"
  ))

colors <- c(
  "Neuron" = "forestgreen",
  "Astrocytes" = "dodgerblue",
  "Oligodendrocytes" = "darkorange",
  "OPCs" = "darkorange",
  "Microglia" = "firebrick",
  "Other" = "goldenrod"
)

plot_title <- if (tolower(gene_myeloid) == tolower(gene_non_myeloid)) {
  sprintf("%s Expression Across Cell Types", gene_for_both)
} else {
  sprintf("Expression Across Cell Types (Non-Myeloid: %s, Myeloid: %s)", gene_non_myeloid, gene_myeloid)
}

ggplot(combined_df, aes(x = cell_ontology_class, y = TP100K, fill = Category)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_wrap(~ cell_group, scales = "free_x") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(
    title = plot_title,
    subtitle = "Single-Cell RNA-seq of FACS sorted cells",
    x = "Cell Type",
    y = "Transcripts per 100,000 (TP100K)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )


mean_out <- combined_df %>%
  group_by(cell_group, cell_ontology_class) %>%
  summarise(mean_TP100K = mean(TP100K, na.rm = TRUE), .groups = "drop") %>%
  arrange(cell_group, cell_ontology_class)

cat("cell_group\tcell_ontology_class\tmean_TP100K\n")
write.table(mean_out, file = "", sep = "\t", row.names = FALSE, quote = FALSE)


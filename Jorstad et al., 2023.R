WORKDIR <- "C:/Users/u0171728/Downloads/JorstadNL"  
setwd(WORKDIR)

genes_file   <- "human_MTG_2018-06-14_genes-rows.csv"
samples_file <- "human_MTG_2018-06-14_samples-columns.csv"
exon_file    <- "human_MTG_2018-06-14_exon-matrix.csv"
intron_file  <- "human_MTG_2018-06-14_intron-matrix.csv"

genes_of_interest <- c("ATP13A2","ATP13A3","ATP13A4","ATP13A5")

need <- c("data.table","dplyr","tidyr","stringr","readr")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, Ncpus = 2)
invisible(lapply(need, library, character.only = TRUE))

genes   <- data.table::fread(genes_file)
samples <- data.table::fread(samples_file)

genes <- genes %>%
  mutate(
    gene_symbol = toupper(gene),
    row_index   = dplyr::row_number()
  )

samples <- samples %>%
  mutate(
    column_index = dplyr::row_number(),
    class   = ifelse(is.na(class),   "Unknown", class),
    cluster = ifelse(is.na(cluster), "Unknown", cluster)
  )

ng <- nrow(genes)
nc <- nrow(samples)
cat("Expect matrices (genes x cells): ", ng, " x ", nc, "\n", sep="")

read_matrix_gxc <- function(path, expect_genes, expect_cells) {
  dt <- suppressWarnings(data.table::fread(path, showProgress = FALSE))
 
  if (nrow(dt) == expect_genes && ncol(dt) == expect_cells + 1) {
    dt <- dt[, -1, with = FALSE]
  }
  if (nrow(dt) != expect_genes || ncol(dt) != expect_cells) {
    stop(sprintf("Matrix dims (%d x %d) != expected (%d x %d) for %s",
                 nrow(dt), ncol(dt), expect_genes, expect_cells, basename(path)))
  }
  
  for (j in seq_len(ncol(dt))) if (!is.numeric(dt[[j]])) dt[[j]] <- as.numeric(dt[[j]])
  dt
}


cat("Reading exon matrix...\n")
exon   <- read_matrix_gxc(exon_file,   ng, nc)
cat("Reading intron matrix...\n")
intron <- read_matrix_gxc(intron_file, ng, nc)

counts <- exon + intron
rm(exon, intron); gc()


libsize <- colSums(counts)
libsize[libsize == 0] <- NA_real_


subtype_order <- c(
  "Excitatory neuron","Inhibitory neuron",
  "Astrocyte","Oligodendrocyte","OPC",
  "Microglia","Endothelial","VLMC/Pericyte",
  "Other/Unknown"
)

samples_sub <- samples %>%
  mutate(
    subtype = case_when(
      class == "Glutamatergic" ~ "Excitatory neuron",
      class == "GABAergic"     ~ "Inhibitory neuron",
      
      class %in% c("Non-neuronal","no class","Unknown") & str_detect(cluster, regex("astro",  TRUE)) ~ "Astrocyte",
      class %in% c("Non-neuronal","no class","Unknown") & str_detect(cluster, regex("\\boligo\\b|oligodend", TRUE)) ~ "Oligodendrocyte",
      class %in% c("Non-neuronal","no class","Unknown") & str_detect(cluster, regex("\\bopc\\b", TRUE)) ~ "OPC",
      class %in% c("Non-neuronal","no class","Unknown") & str_detect(cluster, regex("micro",  TRUE)) ~ "Microglia",
      class %in% c("Non-neuronal","no class","Unknown") & str_detect(cluster, regex("endo",   TRUE)) ~ "Endothelial",
      class %in% c("Non-neuronal","no class","Unknown") & str_detect(cluster, regex("vlmc|peric", TRUE)) ~ "VLMC/Pericyte",
      
      TRUE ~ "Other/Unknown"
    ),
    subtype = factor(subtype, levels = subtype_order)
  )


gene_idx <- which(genes$gene_symbol %in% toupper(genes_of_interest))
if (length(gene_idx) == 0) {
  stop("None of the requested genes found. Check genes_of_interest spelling.")
}
found_genes <- genes$gene_symbol[gene_idx]
cat("Found genes: ", paste(found_genes, collapse=", "), "\n", sep="")


expr_long <- lapply(gene_idx, function(i) {
  raw <- as.numeric(counts[i, ])
  cpm <- 1e6 * raw / libsize
  tibble(
    gene_symbol  = genes$gene_symbol[i],
    column_index = seq_len(nc),
    raw_count    = raw,
    cpm          = cpm,
    logCPM       = log1p(cpm)
  )
}) %>% bind_rows() %>%
  left_join(samples_sub %>% select(column_index, sample_id, class, cluster, subtype),
            by = "column_index")

summ_sub <- expr_long %>%
  group_by(gene_symbol, subtype) %>%
  summarize(
    n_cells       = dplyr::n(),
    mean_logCPM   = mean(logCPM, na.rm = TRUE),
    median_logCPM = median(logCPM, na.rm = TRUE),
    frac_detected = mean(raw_count > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(gene_symbol, subtype)

dir.create("mtg_atp_output", showWarnings = FALSE)
readr::write_csv(summ_sub, "mtg_atp_output/ATP13A2_A3_A4_A5_summary_by_subtype.csv")


subtype_names <- levels(samples_sub$subtype)
subtype_names <- subtype_names[!is.na(subtype_names)]

cat("\n--- GraphPad block (MEAN log1p(CPM)) ---\n")
cat(paste(subtype_names, collapse=","), "\n")

for (g in genes_of_interest) {
  g_up <- toupper(g)
  row <- summ_sub %>%
    filter(gene_symbol == g_up) %>%
    
    right_join(tibble(subtype = factor(subtype_names, levels = subtype_names)), by = "subtype") %>%
    mutate(gene_symbol = g_up) %>%
    arrange(subtype)
  
  vals <- row$mean_logCPM
  vals_txt <- ifelse(is.na(vals), "", format(vals, digits = 6, scientific = FALSE, trim = TRUE))
  
  cat("\n# ", g_up, "\n", sep="")
  cat(paste(vals_txt, collapse=","), "\n")
}


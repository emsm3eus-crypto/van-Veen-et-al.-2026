getwd()
setwd("C:/Users/u0171728/Downloads/comparative")
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

genes_of_interest <- c("ATP13A2","ATP13A3","ATP13A4","ATP13A5")
drop_low_quality  <- TRUE

species_files <- list(
  human   = list(genes="human_LGN_2021_genes-rows.csv",
                 exon ="human_LGN_2021_exon-matrix.csv",
                 meta ="human_LGN_2021_metadata.csv"),
  macaque = list(genes="macaque_LGN_2021_genes-rows.csv",
                 exon ="macaque_LGN_2021_exon-matrix.csv",
                 meta ="macaque_LGN_2021_metadata.csv"),
  mouse   = list(genes="mouse_LGN_2021_genes-rows.csv",
                 exon ="mouse_LGN_2021_exon-matrix.csv",
                 meta ="mouse_LGN_2021_metadata.csv")
)

normalize_id <- function(x){
  x <- as.character(x)
  x <- tolower(x)
  x <- gsub("\\s+", "", x)
  x <- gsub("[\\._-]", "", x)
  x <- gsub("-\\d+$", "", x)  # drop "-1" suffixes
  x
}

auto_log2_if_needed <- function(vec){
  q95 <- suppressWarnings(quantile(vec, 0.95, na.rm = TRUE))
  if (is.finite(q95) && q95 > 40) list(x = log2(vec + 1), label = "log2(CPM + 1) [auto]")
  else                            list(x = vec,            label = "log2(CPM + 1)")
}

pick_gene_symbol_col <- function(gr){
  intersect(c("gene_symbol","gene","symbol","Gene","GeneSymbol","name"), names(gr))[1]
}

pick_subtype_col <- function(meta){
  intersect(c("cluster_label","cluster","subclass_label","subclass","cell_type",
              "celltype","cell_class","class","subtype","cluster_name"), names(meta))[1]
}

best_id_column <- function(meta, expr_ids){
  expr_norm <- normalize_id(expr_ids)
  scores <- sapply(names(meta), function(cn){
    v <- meta[[cn]]
    if (is.list(v)) return(0L)
    sum(normalize_id(v) %in% expr_norm)
  })
  best <- names(scores)[which.max(scores)]
  list(col = best, matches = max(scores), total = length(expr_ids))
}


process_species <- function(sp, files){
 
  if (!all(file.exists(files))) {
    cat("\n[", sp, "] Missing files → skipped.\n", sep = "")
    return(invisible(NULL))
  }
  
  cat("\n================ ", toupper(sp), " ================\n", sep = "")
  

  gr <- fread(files["genes"])
  gene_col <- pick_gene_symbol_col(gr)
  if (is.na(gene_col)) { cat("[", sp, "] No gene symbol column in genes-rows.\n", sep=""); return(invisible(NULL)) }
  gr[[gene_col]] <- toupper(as.character(gr[[gene_col]]))
  present <- intersect(genes_of_interest, unique(gr[[gene_col]]))
  if (!length(present)) { cat("[", sp, "] None of ATP13A2–A5 present.\n", sep=""); return(invisible(NULL)) }
  sel_rows <- which(gr[[gene_col]] %in% present)
  gene_map <- tibble(.row_index = seq_len(nrow(gr)), gene_symbol = gr[[gene_col]])
  

  expr_dt <- fread(files["exon"])

  setnames(expr_dt, 1, "first_col_tag")
  expr_sel <- expr_dt[sel_rows, ]
  
 
  num_cols <- names(expr_sel)[vapply(expr_sel, is.numeric, logical(1))]
  if (!length(num_cols)) {

    cand <- setdiff(names(expr_sel), "first_col_tag")
    for (j in cand) suppressWarnings(expr_sel[[j]] <- as.numeric(expr_sel[[j]]))
    num_cols <- names(expr_sel)[vapply(expr_sel, is.numeric, logical(1))]
  }
  if (!length(num_cols)) { cat("[", sp, "] No numeric expression columns found.\n", sep=""); return(invisible(NULL)) }
  

  auto <- auto_log2_if_needed(as.numeric(as.matrix(expr_sel[, ..num_cols])))
  expr_mat <- matrix(auto$x, nrow = nrow(expr_sel), byrow = FALSE)
  colnames(expr_mat) <- num_cols
  expr_vals <- as.data.table(expr_mat)
  expr_vals[, .row_index := sel_rows]
  

  long <- melt(
    expr_vals,
    id.vars = ".row_index",
    measure.vars = num_cols,
    variable.name = "cell_id",
    value.name = "expr",
    variable.factor = FALSE
  ) |>
    as_tibble() |>
    left_join(gene_map, by = ".row_index") |>
    select(gene_symbol, cell_id, expr)
  

  meta <- fread(files["meta"]) |> as_tibble()
  id_best <- best_id_column(meta, unique(long$cell_id))
  id_col  <- id_best$col
  subtype_col <- pick_subtype_col(meta)
  
  if (is.null(id_col) || is.na(id_col) || is.na(subtype_col)) {
    cat("[", sp, "] Could not detect ID or subtype columns.\n", sep = "")
    return(invisible(NULL))
  }
  
  meta <- meta %>%
    mutate(cell_key_norm = normalize_id(.data[[id_col]]))
  long <- long %>%
    mutate(cell_key_norm = normalize_id(cell_id))
  
  meta_min <- meta %>%
    transmute(cell_key_norm,
              subtype = as.character(.data[[subtype_col]]))
  
  if (drop_low_quality) {
    meta_min <- meta_min %>% filter(!str_detect(subtype, regex("low\\s*quality", TRUE)))
  }
  
  joined <- inner_join(long, meta_min, by = "cell_key_norm")
  

  summary_sp <- joined %>%
    group_by(subtype, gene_symbol) %>%
    summarise(mean_expr = mean(expr, na.rm = TRUE),
              .groups = "drop")
  

  if (!nrow(summary_sp)) {
    cat("[", sp, "] No joined data after matching IDs.\n", sep = "")
    return(invisible(NULL))
  }
  

  sub_order <- summary_sp %>%
    group_by(subtype) %>%
    summarise(overall = mean(mean_expr, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(overall)) %>%
    pull(subtype)
  

  cat(sprintf("[Info] join on '%s' | subtype: '%s' | scale: %s\n",
              id_col, subtype_col, auto$label))
  cat("\n# Copy this header once per gene:\n", paste(sub_order, collapse = ","), "\n", sep = "")
  
  for (g in genes_of_interest) {
    vals <- summary_sp %>%
      filter(gene_symbol == g) %>%
      select(subtype, mean_expr) %>%
      right_join(tibble(subtype = sub_order), by = "subtype") %>%
      arrange(match(subtype, sub_order)) %>%
      pull(mean_expr)
    vals_fmt <- ifelse(is.na(vals), "", formatC(vals, digits = 3, format = "f"))
    cat("\n# ", sp, " — ", g, "\n", paste(vals_fmt, collapse = ","), "\n", sep = "")
  }
  
  invisible(NULL)
}

invisible(lapply(names(species_files), function(sp){
  files <- unlist(species_files[[sp]], use.names = FALSE)
  names(files) <- names(species_files[[sp]])
  process_species(sp, files)
}))

cat("\n=== Done. (GraphPad rows printed above for each species) ===\n")


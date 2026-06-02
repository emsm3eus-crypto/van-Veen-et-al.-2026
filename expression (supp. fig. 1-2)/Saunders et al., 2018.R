
setwd("C:/Users/u0171728/Downloads")
library(ggplot2)

metacell_expr <- readRDS("metacells.BrainCellAtlas_Saunders_version_2018.04.01.rds")
annotations <- readRDS("annotation.BrainCellAtlas_Saunders_version_2018.04.01.rds")

groups <- list(
  `Exc neurons (sbc)` = annotations$tissue_subcluster[
    annotations$class == "NEURON" & grepl("sbc", annotations$common_name, ignore.case = TRUE)
  ],
  `Exc neurons (ca2:ca3)` = annotations$tissue_subcluster[
    annotations$class == "NEURON" & grepl("ca2|ca3", annotations$common_name, ignore.case = TRUE)
  ],
  `Exc neurons (sbc/etcx)` = annotations$tissue_subcluster[
    annotations$class == "NEURON" & grepl("etcx", annotations$common_name, ignore.case = TRUE)
  ],
  `Exc neurons (ca1:sbc:psbc)` = annotations$tissue_subcluster[
    annotations$class == "NEURON" & grepl("ca1|psbc", annotations$common_name, ignore.case = TRUE)
  ],
  `Exc neurons (dg)` = annotations$tissue_subcluster[
    annotations$class == "NEURON" & grepl("dg", annotations$common_name, ignore.case = TRUE)
  ],
  `Inhibitory neurons` = annotations$tissue_subcluster[
    annotations$class == "NEURON" & grepl("inhibitory|gaba", annotations$common_name, ignore.case = TRUE)
  ],
  `Cajal Retzius cells` = annotations$tissue_subcluster[
    grepl("cajal", annotations$common_name, ignore.case = TRUE)
  ],
  `Astrocytes` = annotations$tissue_subcluster[
    annotations$class == "ASTROCYTE"
  ],
  `Oligodendrocytes` = annotations$tissue_subcluster[
    annotations$class == "OLIGODENDROCYTE"
  ],
  `OPCs` = annotations$tissue_subcluster[
    annotations$class == "POLYDENDROCYTE"
  ],
  `Microglia` = annotations$tissue_subcluster[
    annotations$class == "MICROGLIA"
  ],
  `Neural progenitors` = annotations$tissue_subcluster[
    grepl("progenitor", annotations$common_name, ignore.case = TRUE)
  ],
  `Choroid plexus` = annotations$tissue_subcluster[
    grepl("choroid", annotations$common_name, ignore.case = TRUE)
  ],
  `Ependymal cells` = annotations$tissue_subcluster[
    grepl("ependymal", annotations$common_name, ignore.case = TRUE)
  ],
  `Endothelial cells` = annotations$tissue_subcluster[
    grepl("endothelial", annotations$common_name, ignore.case = TRUE)
  ],
  `Fibroblast-like cells` = annotations$tissue_subcluster[
    grepl("fibro", annotations$common_name, ignore.case = TRUE)
  ],
  `Mural` = annotations$tissue_subcluster[
    annotations$class == "MURAL"
  ]
)


neuron_groups <- c(
  "Exc neurons (sbc)", "Exc neurons (ca2:ca3)", "Exc neurons (sbc/etcx)",
  "Exc neurons (ca1:sbc:psbc)", "Exc neurons (dg)", "Inhibitory neurons",
  "Cajal Retzius cells", "Neural progenitors"
)

colors <- c(
  "Neuron" = "forestgreen",
  "Astrocytes" = "dodgerblue",
  "Oligodendrocytes" = "darkorange",
  "OPCs" = "darkorange",
  "Microglia" = "firebrick",
  "Other" = "goldenrod"
)


gene <- "Atp13a5"

df <- data.frame()
for (celltype in names(groups)) {
  mcells <- groups[[celltype]]
  valid_mcells <- mcells[mcells %in% colnames(metacell_expr)]
  if (length(valid_mcells) == 0) next
  
  expr_vals <- metacell_expr[gene, valid_mcells]
  total_umis <- sum(metacell_expr[, valid_mcells])
  tp100k <- 100000 * sum(expr_vals) / total_umis
  
  df <- rbind(df, data.frame(
    CellType = celltype,
    TP100K = tp100k
  ))
}

df$Category <- ifelse(df$CellType %in% neuron_groups, "Neuron",
                      ifelse(df$CellType == "Astrocytes", "Astrocytes",
                             ifelse(df$CellType %in% c("Oligodendrocytes", "OPCs"), "Oligodendrocytes",
                                    ifelse(df$CellType == "Microglia", "Microglia", "Other"))))

df$CellType <- factor(df$CellType, levels = c(
  neuron_groups,
  "Astrocytes",
  "Oligodendrocytes",
  "OPCs",
  "Microglia",
  setdiff(df$CellType, c(neuron_groups, "Astrocytes", "Oligodendrocytes", "OPCs", "Microglia"))
))


ggplot(df, aes(x = CellType, y = TP100K, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(
    title = paste(gene, "Expression Across Cell Types"),
    subtitle = "Data extracted from DropViz.org (Saunders et al., 2018)",
    x = NULL, y = "Transcripts per 100,000 (TP100K)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(face = "plain", hjust = 0.5),
    legend.position = "none"   # remove legend
  )

cat("CellType,TP100K,Category\n")
for (i in 1:nrow(df)) {
  cat(as.character(df$CellType[i]), ",", df$TP100K[i], ",", df$Category[i], "\n", sep = "")
}



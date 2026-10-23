library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


# Load previously saved DESeq2 results
results_df <- read.csv(
  "GSE251911_DESeq2_all_genes.csv",
  stringsAsFactors = FALSE
)

deg_df <- read.csv(
  "GSE251911_DESeq2_DEGs_FDR0.05_log2FC0.58.csv",
  stringsAsFactors = FALSE
)

# Separate upregulated and downregulated genes
up_df <- subset(
  deg_df,
  regulation == "Upregulated"
)

down_df <- subset(
  deg_df,
  regulation == "Downregulated"
)

# Verify imported data
stopifnot(
  all(c("Ensembl_ID", "padj") %in% names(results_df)),
  all(c("Ensembl_ID", "regulation") %in% names(deg_df))
)

cat("Upregulated:", nrow(up_df), "\n")
cat("Downregulated:", nrow(down_df), "\n")


# Define tested genes as the enrichment background
background_genes <- unique(
  results_df$Ensembl_ID[!is.na(results_df$padj)]
)

up_genes <- unique(up_df$Ensembl_ID)
down_genes <- unique(down_df$Ensembl_ID)


# Perform GO Biological Process enrichment
go_up <- enrichGO(
  gene = up_genes,
  universe = background_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 1
)

go_down <- enrichGO(
  gene = down_genes,
  universe = background_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 1
)


# Export complete enrichment results
write.csv(
  as.data.frame(go_up),
  "GSE251911_GO_BP_Upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(go_down),
  "GSE251911_GO_BP_Downregulated.csv",
  row.names = FALSE
)


# Create publication-style GO enrichment dot plots
plot_go <- function(go_result, title, colors, filename) {
  
  df <- as.data.frame(go_result)
  
  if (nrow(df) == 0) {
    message("No enriched terms found: ", title)
    return(invisible(NULL))
  }
  
  df <- head(df[order(df$p.adjust), ], 12)
  
  df$ratio <- vapply(
    strsplit(df$GeneRatio, "/"),
    function(x) as.numeric(x[1]) / as.numeric(x[2]),
    numeric(1)
  )
  
  df$term <- vapply(
    df$Description,
    function(x) paste(strwrap(x, width = 42), collapse = "\n"),
    character(1)
  )
  
  df$term <- factor(df$term, levels = rev(df$term))
  df$significance <- -log10(pmax(df$p.adjust, 1e-300))
  
  p <- ggplot(
    df,
    aes(x = ratio, y = term)
  ) +
    geom_point(
      aes(size = Count, color = significance),
      alpha = 0.95
    ) +
    scale_color_gradient(
      low = colors[1],
      high = colors[2],
      name = "-log10(FDR)"
    ) +
    scale_size_continuous(
      range = c(3, 8),
      name = "Gene count"
    ) +
    scale_x_continuous(
      labels = function(x) paste0(round(x * 100, 1), "%")
    ) +
    labs(
      title = title,
      x = "Gene ratio",
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 17,
        margin = margin(b = 18)
      ),
      axis.text.y = element_text(
        size = 11,
        color = "#253244"
      ),
      axis.title.x = element_text(
        face = "bold"
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  print(p)
  
  ggsave(
    filename,
    plot = p,
    width = 11,
    height = 7.5,
    dpi = 320,
    bg = "white"
  )
}


# GO enrichment of upregulated genes
plot_go(
  go_up,
  title = "GO Biological Processes | Upregulated",
  colors = c("#F2B4A8", "#B83235"),
  filename = "GSE251911_GO_BP_Upregulated.png"
)


# GO enrichment of downregulated genes
plot_go(
  go_down,
  title = "GO Biological Processes | Downregulated",
  colors = c("#A6CDE2", "#16588A"),
  filename = "GSE251911_GO_BP_Downregulated.png"
)


# Summarize significant GO terms
cat("\nGO ENRICHMENT COMPLETED\n")
cat("Upregulated GO terms:", nrow(as.data.frame(go_up)), "\n")
cat("Downregulated GO terms:", nrow(as.data.frame(go_down)), "\n")

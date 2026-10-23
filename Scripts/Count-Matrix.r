library(DESeq2)
library(apeglm)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(pheatmap)


# Download raw RNA-seq count matrix

url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE251nnn/GSE251911/suppl/",
  "GSE251911_CYT_RNA-seq_Raw_Count_Matrix.txt.gz"
)

file <- "GSE251911_CYT_RNA-seq_Raw_Count_Matrix.txt.gz"

download.file(url, destfile = file, mode = "wb")

# Load into data frame
counts_df <- read.delim(
  gzfile(file),
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Inspect dataset
dim(counts_df)
colnames(counts_df)
head(counts_df[, 1:min(5, ncol(counts_df))])


# Load cleaned metadata
metadata_clean <- read.csv(
  "GSE251911_clean_metadata.csv",
  stringsAsFactors = FALSE
)

# Create matching sample names
sample_names <- paste0(
  metadata_clean$islet_id,
  "_24Hr_",
  metadata_clean$condition
)

# Verify sample matching
stopifnot(
  nrow(metadata_clean) == 38,
  !anyDuplicated(sample_names),
  setequal(sample_names, colnames(counts_df)[-1]),
  all(table(
    metadata_clean$donor_id,
    metadata_clean$condition
  ) == 1)
)

# Align count matrix with metadata
counts <- as.matrix(
  counts_df[, sample_names, drop = FALSE]
)

# Assign Ensembl gene IDs
rownames(counts) <- counts_df$Ensembl_ID

# Assign GEO sample IDs
colnames(counts) <- metadata_clean$sample_id

# Prepare metadata
metadata_clean$donor_id <- factor(
  metadata_clean$donor_id
)

metadata_clean$condition <- factor(
  metadata_clean$condition,
  levels = c("Untreated", "Cytokines")
)

rownames(metadata_clean) <- metadata_clean$sample_id

# Final alignment check
stopifnot(
  identical(
    colnames(counts),
    rownames(metadata_clean)
  )
)

# Save aligned count matrix
write.csv(
  data.frame(
    Ensembl_ID = rownames(counts),
    counts,
    check.names = FALSE
  ),
  "GSE251911_aligned_counts.csv",
  row.names = FALSE
)

# Save aligned metadata
write.csv(
  metadata_clean,
  "GSE251911_aligned_metadata.csv",
  row.names = FALSE
)

# Inspect
dim(counts)
dim(metadata_clean)

table(metadata_clean$condition)

identical(
  colnames(counts),
  rownames(metadata_clean)
)




# GSE251911 | Hierarchical Clustering and PCA

# Verify count matrix and metadata alignment
stopifnot(
  identical(colnames(counts), rownames(metadata_clean))
)


# Construct DESeq2 dataset using paired experimental design
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata_clean,
  design = ~ donor_id + condition
)


# Exclude all-zero genes for exploratory analysis
dds_qc <- dds[
  rowSums(DESeq2::counts(dds)) > 0,
]


# Apply variance stabilizing transformation
vsd <- vst(
  dds_qc,
  blind = FALSE
)


# Calculate sample-to-sample expression correlations
sample_cor <- cor(
  SummarizedExperiment::assay(vsd),
  method = "pearson"
)


# Perform hierarchical clustering using correlation distance
sample_hclust <- hclust(
  as.dist(1 - sample_cor),
  method = "complete"
)


# Prepare sample labels with treatment information
sample_labels <- paste0(
  metadata_clean$islet_id,
  "_",
  ifelse(
    metadata_clean$condition == "Cytokines",
    "CYT",
    "UNT"
  )
)


# Visualize hierarchical clustering dendrogram
png(
  "GSE251911_hierarchical_clustering.png",
  width = 3000,
  height = 1800,
  res = 300
)

plot(
  sample_hclust,
  labels = sample_labels,
  main = "GSE251911 | Hierarchical Clustering",
  xlab = "",
  sub = "",
  ylab = "1 - Pearson Correlation",
  hang = -1,
  cex = 0.7
)

dev.off()


# Perform PCA using the 500 most variable genes
pca_df <- plotPCA(
  vsd,
  intgroup = "condition",
  ntop = 500,
  returnData = TRUE
)

percent_var <- round(
  100 * attr(pca_df, "percentVar"),
  1
)


# Visualize principal component analysis
p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = condition
  )
) +
  geom_point(
    size = 3.5,
    alpha = 0.9
  ) +
  labs(
    title = "GSE251911 | Principal Component Analysis",
    x = paste0(
      "PC1: ", percent_var[1], "% variance"
    ),
    y = paste0(
      "PC2: ", percent_var[2], "% variance"
    ),
    color = "Condition"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank()
  )

print(p_pca)

ggsave(
  "GSE251911_PCA.png",
  plot = p_pca,
  width = 9,
  height = 7,
  dpi = 300
)


# Export principal component coordinates
write.csv(
  pca_df,
  "GSE251911_PCA_coordinates.csv",
  row.names = FALSE
)


# Display analysis summary
cat("\nHIERARCHICAL CLUSTERING AND PCA COMPLETED\n")

cat("Samples analyzed:", ncol(dds_qc), "\n")
cat("Genes retained:", nrow(dds_qc), "\n")
cat("PC1 variance:", percent_var[1], "%\n")
cat("PC2 variance:", percent_var[2], "%\n")


# Define statistical and expression thresholds
fdr_cutoff <- 0.05
lfc_cutoff <- 0.58
min_samples <- min(table(metadata_clean$condition))

stopifnot(
  identical(colnames(counts), rownames(metadata_clean)),
  identical(colnames(dds), colnames(counts)),
  min_samples == 19,
  identical(
    levels(metadata_clean$condition),
    c("Untreated", "Cytokines")
  )
)


# Filter low-count genes before fitting DESeq2
keep <- rowSums(DESeq2::counts(dds) >= 10) >= min_samples
dds_de <- dds[keep, ]

cat("Genes before filtering:", nrow(dds), "\n")
cat("Genes after filtering:", nrow(dds_de), "\n")


# Fit a negative binomial model accounting for paired donors
dds_de <- DESeq(dds_de)

coef_name <- "condition_Cytokines_vs_Untreated"
stopifnot(coef_name %in% resultsNames(dds_de))


# Obtain standard DESeq2 results and shrink log2 fold changes
res <- results(
  dds_de,
  name = coef_name,
  alpha = fdr_cutoff
)

res_shrunk <- lfcShrink(
  dds_de,
  coef = coef_name,
  res = res,
  type = "apeglm"
)


# Add gene symbols and retain unshrunken fold changes
results_df <- as.data.frame(res_shrunk)

results_df$Ensembl_ID <- rownames(results_df)
results_df$log2FC_unshrunken <- res$log2FoldChange
results_df$Wald_statistic <- res$stat

ensembl_ids <- sub(
  "\\.[0-9]+$",
  "",
  results_df$Ensembl_ID
)

symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(ensembl_ids),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

results_df$gene_symbol <- unname(
  symbol_map[ensembl_ids]
)

results_df$gene_label <- ifelse(
  is.na(results_df$gene_symbol) |
    results_df$gene_symbol == "",
  results_df$Ensembl_ID,
  results_df$gene_symbol
)


# Classify DEGs using FDR and shrunken effect sizes
results_df$regulation <- "Not significant"

up <- !is.na(results_df$padj) &
  results_df$padj < fdr_cutoff &
  !is.na(results_df$log2FoldChange) &
  results_df$log2FoldChange >= lfc_cutoff

down <- !is.na(results_df$padj) &
  results_df$padj < fdr_cutoff &
  !is.na(results_df$log2FoldChange) &
  results_df$log2FoldChange <= -lfc_cutoff

results_df$regulation[up] <- "Upregulated"
results_df$regulation[down] <- "Downregulated"

results_df$regulation <- factor(
  results_df$regulation,
  levels = c(
    "Downregulated",
    "Not significant",
    "Upregulated"
  )
)

results_df <- results_df[
  order(results_df$padj, na.last = TRUE),
]

up_df <- subset(
  results_df,
  regulation == "Upregulated"
)

down_df <- subset(
  results_df,
  regulation == "Downregulated"
)

deg_df <- subset(
  results_df,
  regulation != "Not significant"
)


# Export standard and shrunken differential-expression results
write.csv(
  data.frame(
    Ensembl_ID = rownames(res),
    as.data.frame(res)
  ),
  "GSE251911_DESeq2_unshrunken_results.csv",
  row.names = FALSE
)

write.csv(
  results_df,
  "GSE251911_DESeq2_all_genes.csv",
  row.names = FALSE
)

write.csv(
  deg_df,
  "GSE251911_DESeq2_DEGs_FDR0.05_log2FC0.58.csv",
  row.names = FALSE
)

write.csv(
  up_df,
  "GSE251911_DESeq2_upregulated.csv",
  row.names = FALSE
)

write.csv(
  down_df,
  "GSE251911_DESeq2_downregulated.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    Ensembl_ID = rownames(dds_de),
    DESeq2::counts(dds_de, normalized = TRUE),
    check.names = FALSE
  ),
  "GSE251911_DESeq2_normalized_counts.csv",
  row.names = FALSE
)


# Display differential-expression summary
cat(
  "\nFDR-significant:",
  sum(results_df$padj < fdr_cutoff, na.rm = TRUE),
  "\n"
)

cat("Upregulated:", nrow(up_df), "\n")
cat("Downregulated:", nrow(down_df), "\n")
cat("Total selected DEGs:", nrow(deg_df), "\n")


# figure colors
reg_colors <- c(
  "Downregulated" = "#337EA1",
  "Not significant" = "#AAB2BD",
  "Upregulated" = "#D86551"
)


# Save figures in high-resolution PNG and vector PDF
save_figure <- function(p, stem, width = 11, height = 8) {
  ggsave(
    paste0(stem, ".png"),
    p,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )
  
  ggsave(
    paste0(stem, ".pdf"),
    p,
    width = width,
    height = height,
    bg = "white"
  )
}


# Prepare volcano plot with adjusted p-values and shrunken log2FC
volcano_df <- subset(
  results_df,
  is.finite(log2FoldChange) & !is.na(padj)
)

volcano_df$minus_log10_fdr <- -log10(
  pmax(volcano_df$padj, .Machine$double.xmin)
)

volcano_ceiling <- min(
  80,
  max(
    20,
    unname(quantile(
      volcano_df$minus_log10_fdr,
      0.995,
      na.rm = TRUE
    ))
  )
)

volcano_df$plot_y <- pmin(
  volcano_df$minus_log10_fdr,
  volcano_ceiling
)

label_ids <- unique(c(
  head(up_df$Ensembl_ID, 6),
  head(down_df$Ensembl_ID, 6)
))

label_df <- subset(
  volcano_df,
  Ensembl_ID %in% label_ids
)


# volcano plot
p_volcano <- ggplot(
  volcano_df,
  aes(x = log2FoldChange, y = plot_y)
) +
  geom_point(
    aes(color = regulation),
    size = 1,
    alpha = 0.6
  ) +
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = "dashed",
    linewidth = 0.4,
    color = "#606B79"
  ) +
  geom_hline(
    yintercept = -log10(fdr_cutoff),
    linetype = "dashed",
    linewidth = 0.4,
    color = "#606B79"
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = gene_label),
    color = "#17212F",
    size = 3.5,
    box.padding = 0.55,
    point.padding = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 2026
  ) +
  scale_color_manual(
    values = reg_colors,
    drop = FALSE
  ) +
  labs(
    title = "Cytokine-induced differential gene expression",
    subtitle = paste0(
      "Human pancreatic islets | 19 paired donors | Up: ",
      nrow(up_df),
      " | Down: ",
      nrow(down_df)
    ),
    x = "Shrunken log2 fold change (Cytokines / Untreated)",
    y = "-log10(adjusted p-value)",
    color = "Gene classification",
    caption = paste0(
      "FDR < 0.05; |shrunken log2FC| >= 0.58. ",
      "Values above ",
      round(volcano_ceiling, 1),
      " on y-axis are capped for display."
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(color = "#526070"),
    legend.position = "bottom",
    axis.title = element_text(face = "bold")
  )

print(p_volcano)

save_figure(
  p_volcano,
  "GSE251911_DESeq2_volcano"
)


# Prepare MA plot using mean normalized counts and shrunken log2FC
ma_df <- subset(
  results_df,
  is.finite(log2FoldChange) & baseMean > 0
)

ma_limit <- min(
  8,
  max(
    3,
    unname(quantile(
      abs(ma_df$log2FoldChange),
      0.995,
      na.rm = TRUE
    ))
  )
)


# Generate publication-style MA plot
p_ma <- ggplot(
  ma_df,
  aes(x = baseMean, y = log2FoldChange)
) +
  geom_point(
    aes(color = regulation),
    size = 1,
    alpha = 0.55
  ) +
  geom_hline(
    yintercept = 0,
    color = "#4B5563",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = c(-lfc_cutoff, lfc_cutoff),
    color = "#7B8794",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_x_log10(
    labels = scales::label_comma()
  ) +
  scale_color_manual(
    values = reg_colors,
    drop = FALSE
  ) +
  coord_cartesian(
    ylim = c(-ma_limit, ma_limit)
  ) +
  labs(
    title = "MA plot | Cytokines versus Untreated",
    subtitle = "DESeq2 mean normalized counts and shrunken log2FC",
    x = "Mean normalized count (log10 scale)",
    y = "Shrunken log2 fold change",
    color = "Gene classification",
    caption = paste0(
      "Display range +/-",
      round(ma_limit, 1),
      " log2FC; complete values are in the results table."
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "#526070"),
    legend.position = "bottom",
    axis.title = element_text(face = "bold")
  )

print(p_ma)

save_figure(
  p_ma,
  "GSE251911_DESeq2_MA_plot",
  height = 7.5
)


# Visualize top 30 DEGs using existing VST-transformed data
heatmap_ids <- head(
  deg_df$Ensembl_ID,
  30
)

if (length(heatmap_ids) >= 2) {
  
  heatmap_mat <- SummarizedExperiment::assay(vsd)[
    heatmap_ids,
    ,
    drop = FALSE
  ]
  
  rownames(heatmap_mat) <- make.unique(
    deg_df$gene_label[
      match(heatmap_ids, deg_df$Ensembl_ID)
    ]
  )
  
  colnames(heatmap_mat) <- paste0(
    metadata_clean$islet_id,
    "_",
    ifelse(
      metadata_clean$condition == "Cytokines",
      "CYT",
      "UNT"
    )
  )
  
  annotation_df <- data.frame(
    Condition = metadata_clean$condition,
    row.names = colnames(heatmap_mat)
  )
  
  # Export heatmap in PNG and PDF formats
  for (extension in c("png", "pdf")) {
    
    pheatmap(
      heatmap_mat,
      scale = "row",
      clustering_method = "complete",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      annotation_col = annotation_df,
      annotation_colors = list(
        Condition = c(
          Untreated = "#8B96A6",
          Cytokines = "#D86551"
        )
      ),
      color = colorRampPalette(
        c("#276B98", "#F7F7F7", "#CB5C4C")
      )(101),
      border_color = NA,
      show_colnames = FALSE,
      fontsize_row = 9,
      main = paste0(
        "Top ",
        length(heatmap_ids),
        " DEGs | VST row z-scores"
      ),
      filename = paste0(
        "GSE251911_DESeq2_top30_DEG_heatmap.",
        extension
      ),
      width = 11,
      height = 8
    )
  }
  
} else {
  message(
    "Fewer than two DEGs passed the cutoffs; heatmap skipped."
  )
}


# model and software versions
saveRDS(
  dds_de,
  "GSE251911_DESeq2_fitted_model.rds"
)

capture.output(
  sessionInfo(),
  file = "GSE251911_R_sessionInfo.txt"
)

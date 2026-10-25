# Colorectal Cancer Gene Expression Analysis
# Step 1: Load required packages

library(DESeq2)

# Confirm that DESeq2 is available
packageVersion("DESeq2")
# Step 2: Load the GSE251845 raw count matrix

counts_file <- "data/raw/GSE251845_htseq_raw_counts.csv.gz"

counts <- read.csv(
  counts_file,
  row.names = 1,
  check.names = FALSE
)

# Inspect the matrix dimensions
dim(counts)
# Step 3: Inspect the first row names

head(rownames(counts), 10)
# Step 4: Find HTSeq summary rows

grep("^__", rownames(counts), value = TRUE)
# Step 5: Remove HTSeq summary rows

summary_rows <- grep("^__", rownames(counts), value = TRUE)

counts <- counts[!rownames(counts) %in% summary_rows, ]

# Check the matrix dimensions after removal
dim(counts)
# Step 6: Validate raw gene counts

count_matrix <- as.matrix(counts)

cat("Are all counts non-negative?",
    all(count_matrix >= 0), "\n")

cat("Are all counts integers?",
    all(count_matrix == floor(count_matrix)), "\n")
    # Step 7: Load sample metadata

metadata <- read.csv(
  "data/gse251845_sample_metadata.csv",
  stringsAsFactors = FALSE
)

# Inspect the first six rows
head(metadata)
# Step 8: Check sample-name matching

cat(
  "All metadata samples match count columns:",
  all(metadata$sample %in% colnames(counts)),
  "\n"
)

cat(
  "Number of metadata samples:",
  nrow(metadata),
  "\n"
)

cat(
  "Number of count-matrix samples:",
  ncol(counts),
  "\n"
)
# Step 9: Verify matched patient pairs

metadata$condition <- factor(
  metadata$condition,
  levels = c("Normal", "Tumor")
)

pair_table <- table(
  metadata$patient_id,
  metadata$condition
)
metadata$patient_id <- factor(metadata$patient_id)
print(pair_table)

cat(
  "\nDoes every patient have one Normal and one Tumor sample?",
  all(pair_table == 1),
  "\n"
)
# Step 10: Prepare the DESeq2 dataset

# Ensure sample order matches the count matrix columns
metadata <- metadata[
  match(colnames(counts), metadata$sample),
]

rownames(metadata) <- metadata$sample

# Create the DESeq2 dataset with a paired design
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = metadata,
  design = ~ patient_id + condition
)

# Inspect the dataset
dds
# Step 11: Count genes with zero counts in every sample

zero_count_genes <- rowSums(counts) == 0

cat(
  "Genes with zero counts in all samples:",
  sum(zero_count_genes),
  "\n"
)

cat(
  "Genes with at least one count:",
  sum(!zero_count_genes),
  "\n"
)
# Step 12: Remove genes with zero counts in every sample

keep_genes <- rowSums(counts) > 0

counts_filtered <- counts[keep_genes, ]

cat("Genes before filtering:", nrow(counts), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")
# Step 13: Create DESeq2 dataset with filtered counts

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_filtered)),
  colData = metadata,
  design = ~ patient_id + condition
)

cat("Genes in DESeq2 dataset:", nrow(dds), "\n")
cat("Samples in DESeq2 dataset:", ncol(dds), "\n")
# Step 14: Run differential expression analysis

dds <- DESeq(dds)

# Extract Tumor versus Normal results
res <- results(
  dds,
  contrast = c("condition", "Tumor", "Normal")
)

# Inspect the results
summary(res)
# Step 15: Inspect the top differentially expressed genes

res_ordered <- res[order(res$padj), ]

head(res_ordered, 10)
# Step 16: Summarize adjusted p-values

cat(
  "Genes with an available adjusted p-value:",
  sum(!is.na(res$padj)),
  "\n"
)

cat(
  "Genes with adjusted p-value below 0.05:",
  sum(res$padj < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Genes with adjusted p-value below 0.1:",
  sum(res$padj < 0.1, na.rm = TRUE),
  "\n"
)
# Step 17: Export differential expression results

results_table <- as.data.frame(res)
results_table$gene_id <- rownames(results_table)

results_table <- results_table[
  order(results_table$padj, na.last = TRUE),
]

dir.create("results", showWarnings = FALSE)

write.csv(
  results_table,
  "results/GSE251845_deseq2_results.csv",
  row.names = FALSE
)

cat("Results saved to results/GSE251845_deseq2_results.csv\n")
# Step 18: Quality control - inspect size factors

size_factors <- sizeFactors(dds)

cat("Summary of DESeq2 size factors:\n")
print(summary(size_factors))
# Step 19: Transform counts for exploratory QC

vsd <- vst(dds, blind = FALSE)

# Calculate principal components
pca_data <- plotPCA(
  vsd,
  intgroup = c("condition", "patient_id"),
  returnData = TRUE
)

# Inspect the PCA coordinates
head(pca_data)
# Step 20: Plot PCA by tissue condition

p <- plotPCA(
  vsd,
  intgroup = "condition"
)

print(p)
# Save PCA plot
ggplot2::ggsave(
  filename = "results/GSE251845_PCA_condition.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

cat("PCA plot saved successfully.\n")
# Step 23: Check library sizes

library_sizes <- colSums(counts_filtered)

library_size_table <- data.frame(
  sample = names(library_sizes),
  library_size = as.numeric(library_sizes),
  condition = metadata[names(library_sizes), "condition"]
)

library_size_table <- library_size_table[
  order(library_size_table$library_size),
]

print(library_size_table)

write.csv(
  library_size_table,
  "results/GSE251845_library_sizes.csv",
  row.names = FALSE
)
# Step 24: Boxplot of normalized gene expression

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Convert to log2 scale
log_normalized_counts <- log2(normalized_counts + 1)

# Prepare data for boxplot
boxplot(
  log_normalized_counts,
  outline = FALSE,
  las = 2,
  main = "Normalized Gene Expression Across Samples",
  ylab = "log2 Normalized Counts",
  xlab = "Samples",
  cex.axis = 0.6
)

# Save boxplot
png(
  "results/GSE251845_expression_boxplot.png",
  width = 2400,
  height = 1400,
  res = 200
)

boxplot(
  log_normalized_counts,
  outline = FALSE,
  las = 2,
  main = "Normalized Gene Expression Across Samples",
  ylab = "log2 Normalized Counts",
  xlab = "Samples",
  cex.axis = 0.6
)

dev.off()

cat("Expression boxplot saved successfully.\n")
# Step 26: Paired PCA plot

library(ggplot2)

# Extract PCA coordinates
pca_data <- plotPCA(
  vsd,
  intgroup = c("condition", "patient_id"),
  returnData = TRUE
)

# Add patient ID and condition
pca_data$patient_id <- metadata[
  match(pca_data$name, metadata$sample),
  "patient_id"
]

pca_data$condition <- metadata[
  match(pca_data$name, metadata$sample),
  "condition"
]

# Calculate percentage variance
percentVar <- round(
  100 * attr(pca_data, "percentVar")
)

# Create paired PCA plot
paired_pca <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = condition)
) +
  geom_line(
    aes(group = patient_id),
    color = "grey60",
    linewidth = 0.6
  ) +
  geom_point(size = 3) +
  labs(
    title = "Paired PCA of Colorectal Cancer Samples",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Condition"
  ) +
  theme_bw()

# Save plot
ggsave(
  "results/GSE251845_paired_PCA.png",
  plot = paired_pca,
  width = 9,
  height = 7,
  dpi = 300
)

cat("Paired PCA plot saved successfully.\n")
# Step 28: Calculate paired PCA distances

patient_groups <- split(pca_data, pca_data$patient_id)

paired_distances <- do.call(
  rbind,
  lapply(patient_groups, function(x) {

    tumor <- x[x$condition == "Tumor", ]
    normal <- x[x$condition == "Normal", ]

    data.frame(
      patient_id = as.character(tumor$patient_id),
      delta_PC1 = tumor$PC1 - normal$PC1,
      delta_PC2 = tumor$PC2 - normal$PC2,
      distance = sqrt(
        (tumor$PC1 - normal$PC1)^2 +
        (tumor$PC2 - normal$PC2)^2
      )
    )
  })
)

# Sort by distance, largest first
paired_distances <- paired_distances[
  order(-paired_distances$distance),
]

print(paired_distances)

write.csv(
  paired_distances,
  "results/GSE251845_paired_PCA_distances.csv",
  row.names = FALSE
)

cat("Paired PCA distances calculated successfully.\n")
# Step 30: Sample-to-sample correlation

library(pheatmap)

# Calculate Pearson correlation
sample_cor <- cor(
  log_normalized_counts,
  method = "pearson"
)

# Check correlation values
cat("Correlation summary:\n")
print(summary(as.vector(sample_cor)))

# Generate and save heatmap
dir.create("results", showWarnings = FALSE, recursive = TRUE)
pheatmap(
  sample_cor,
  annotation_col = data.frame(
    Condition = metadata$condition,
    row.names = metadata$sample
  ),
  annotation_row = data.frame(
    Condition = metadata$condition,
    row.names = metadata$sample
  ),
  main = "Sample-to-Sample Correlation",
  border_color = NA,
  filename = "results/GSE251845_sample_correlation.png",
  width = 10,
  height = 9
)

cat("Correlation heatmap saved successfully.\n")
# Step 32: Investigate potential outlier samples

samples_to_check <- c(
  "50N_htseq.out",
  "47C_htseq.out"
)

cat("\nOutlier sample details:\n")

print(
  library_size_table[
    library_size_table$sample %in% samples_to_check,
  ]
)

cat("\nCorrelation with all other samples:\n")

print(
  sample_cor[samples_to_check, ]
)
# Step 34: Verify flagged patient pairs

patients_to_check <- c("47", "50")

cat("\nChecking patient metadata:\n")

print(
  metadata[
    metadata$patient_id %in% patients_to_check,
    c("sample", "patient_id", "condition")
  ]
)

cat("\nChecking size factors:\n")

print(
  sizeFactors(dds)[
    c("47C_htseq.out", "47N_htseq.out",
      "50C_htseq.out", "50N_htseq.out")
  ]
)


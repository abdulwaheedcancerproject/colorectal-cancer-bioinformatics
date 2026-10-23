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
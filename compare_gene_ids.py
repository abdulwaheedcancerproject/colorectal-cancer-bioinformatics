import pandas as pd

file_1 = "data/raw/GSE251845_htseq_raw_counts.csv.gz"
file_2 = "data/raw/GSE196006_raw_counts.csv.gz"

counts_1 = pd.read_csv(file_1, index_col=0)
counts_2 = pd.read_csv(file_2, index_col=0)

summary_rows = [
    "__no_feature",
    "__ambiguous",
    "__too_low_aQual",
    "__not_aligned",
    "__alignment_not_unique",
]

counts_1 = counts_1.drop(index=summary_rows, errors="ignore")
counts_2 = counts_2.drop(index=summary_rows, errors="ignore")

genes_1 = set(counts_1.index)
genes_2 = set(counts_2.index)

shared_genes = genes_1 & genes_2

print("GSE251845 gene IDs:", len(genes_1))
print("GSE196006 gene IDs:", len(genes_2))
print("Shared gene IDs:", len(shared_genes))
print(
    "Shared genes as percentage of GSE251845:",
    round(100 * len(shared_genes) / len(genes_1), 2),
)
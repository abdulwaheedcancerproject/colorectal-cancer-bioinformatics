import pandas as pd

file_path = "data/raw/GSE251845_htseq_raw_counts.csv.gz"
counts = pd.read_csv(file_path, index_col=0)

summary_rows = [
    "__no_feature",
    "__ambiguous",
    "__too_low_aQual",
    "__not_aligned",
    "__alignment_not_unique"
]

counts_clean = counts.drop(index=summary_rows)
library_sizes = counts_clean.sum(axis=0)

print("Five smallest samples:")
print(library_sizes.sort_values().head(5))

print("\nFive largest samples:")
print(library_sizes.sort_values(ascending=False).head(5))

cancer_samples = [
    sample for sample in library_sizes.index
    if sample.lower().endswith("c_htseq.out")
]

normal_samples = [
    sample for sample in library_sizes.index
    if sample.lower().endswith("n_htseq.out")
]

print("\nAverage cancer library size:")
print(library_sizes[cancer_samples].mean())

print("\nAverage normal library size:")
print(library_sizes[normal_samples].mean())
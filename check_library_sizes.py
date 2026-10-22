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

print("Total counts per sample:")
print(library_sizes)

print("\nSmallest library size:")
print(library_sizes.min())

print("\nLargest library size:")
print(library_sizes.max())

print("\nAverage library size:")
print(library_sizes.mean())
import pandas as pd

file_path = "data/raw/GSE251845_htseq_raw_counts.csv.gz"
counts = pd.read_csv(file_path, index_col=0)

print("Original dimensions:", counts.shape)

summary_rows = [
    "__no_feature",
    "__ambiguous",
    "__too_low_aQual",
    "__not_aligned",
    "__alignment_not_unique"
]

counts_clean = counts.drop(index=summary_rows, errors="raise")

print("Dimensions after removing summary rows:", counts_clean.shape)

print("\nAre all remaining row names gene IDs?")
print(counts_clean.index.str.startswith("ENSG").all())

print("\nAre gene IDs unique?")
print(counts_clean.index.is_unique)
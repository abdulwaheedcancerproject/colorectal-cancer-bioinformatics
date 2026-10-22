
import pandas as pd

file_path = "data/raw/GSE251845_htseq_raw_counts.csv.gz"

counts = pd.read_csv(file_path, index_col=0)

print("Matrix dimensions:", counts.shape)

print("\nFirst five genes:")
print(counts.head())

print("\nFirst ten sample names:")
print(counts.columns[:10].tolist())



import pandas as pd

file_path = "data/raw/GSE251845_htseq_raw_counts.csv.gz"

counts = pd.read_csv(file_path, index_col=0)

sample_names = counts.columns.tolist()

cancer_samples = [s for s in sample_names if s.endswith("C_htseq.out")]
normal_samples = [s for s in sample_names if s.endswith("N_htseq.out")]

print("Cancer samples:", len(cancer_samples))
print("Normal samples:", len(normal_samples))

print("\nCancer samples:")
print(cancer_samples)

print("\nNormal samples:")
print(normal_samples)

print("\nLast five rows of the count matrix:")
print(counts.tail())

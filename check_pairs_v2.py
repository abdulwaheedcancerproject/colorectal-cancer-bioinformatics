import pandas as pd

file_path = "data/raw/GSE251845_htseq_raw_counts.csv.gz"
counts = pd.read_csv(file_path, index_col=0)

sample_names = counts.columns.tolist()

cancer_ids = {
    name[:-len("_htseq.out")][:-1]
    for name in sample_names
    if name.lower().endswith("c_htseq.out")
}

normal_ids = {
    name[:-len("_htseq.out")][:-1]
    for name in sample_names
    if name.lower().endswith("n_htseq.out")
}

print("Number of cancer samples:", len(cancer_ids))
print("Number of normal samples:", len(normal_ids))

print("\nCancer samples without a normal match:")
print(sorted(cancer_ids - normal_ids))

print("\nNormal samples without a cancer match:")
print(sorted(normal_ids - cancer_ids))

print("\nAll sample names:")
print(sample_names)
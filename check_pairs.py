
import pandas as pd

file_path = "data/raw/GSE251845_htseq_raw_counts.csv.gz"
counts = pd.read_csv(file_path, index_col=0)

sample_names = counts.columns.tolist()

cancer_ids = {
    name.removesuffix("C_htseq.out")
    for name in sample_names
    if name.endswith("C_htseq.out")
}

normal_ids = {
    name.removesuffix("N_htseq.out")
    for name in sample_names
    if name.endswith("N_htseq.out")
}

print("Cancer samples without a normal match:")
print(sorted(cancer_ids - normal_ids))

print("\nNormal samples without a cancer match:")
print(sorted(normal_ids - cancer_ids))

print("\nSample names not matching either pattern:")
print([
    name for name in sample_names
    if not name.endswith("C_htseq.out")
    and not name.endswith("N_htseq.out")
])

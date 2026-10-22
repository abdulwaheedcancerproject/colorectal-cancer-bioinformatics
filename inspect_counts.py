
from pathlib import Path
import gzip

data_folder = Path("data/raw")

for file in data_folder.glob("*.csv.gz"):
    print("\nDataset:", file.name)

    with gzip.open(file, "rt", encoding="utf-8") as f:
        for i in range(5):
            line = f.readline()
            print(line[:200].rstrip())
import pandas as pd

file_path = "data/raw/GSE251845_htseq_raw_counts.csv.gz"
counts = pd.read_csv(file_path, index_col=0)

print("Matrix dimensions:", counts.shape)

print("\nFirst five gene IDs:")
print(counts.index[:5].tolist())

print("\nLast ten row IDs:")
print(counts.index[-10:].tolist())

print("\nAre all count values numeric?")
print(counts.apply(pd.to_numeric, errors="coerce").notna().all().all())

print("\nAre all count values non-negative?")
print((counts.select_dtypes(include="number") >= 0).all().all())            

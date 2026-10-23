import pandas as pd

files = {
    "GSE251845": "data/raw/GSE251845_htseq_raw_counts.csv.gz",
    "GSE196006": "data/raw/GSE196006_raw_counts.csv.gz",
}

summary_rows = [
    "__no_feature",
    "__ambiguous",
    "__too_low_aQual",
    "__not_aligned",
    "__alignment_not_unique",
]

for dataset, file_path in files.items():
    counts = pd.read_csv(file_path, index_col=0)
    counts = counts.drop(index=summary_rows, errors="ignore")

    print(f"\n{dataset}")
    print("Data types:")
    print(counts.dtypes.value_counts())

    print(
        "All values are non-negative:",
        (counts >= 0).all().all()
    )

    print(
        "All values are integers:",
        (counts % 1 == 0).all().all()
    )
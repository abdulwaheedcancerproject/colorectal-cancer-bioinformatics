import pandas as pd

file_path = "data/raw/GSE251845_htseq_raw_counts.csv.gz"
counts = pd.read_csv(file_path, index_col=0)

records = []

for sample in counts.columns:
    sample_id = sample.split("_")[0]
    patient_id = sample_id[:-1]
    condition_code = sample_id[-1].upper()

    if condition_code == "C":
        condition = "Tumor"
    elif condition_code == "N":
        condition = "Normal"
    else:
        condition = "Unknown"

    records.append({
        "sample": sample,
        "patient_id": patient_id,
        "condition": condition,
    })

metadata = pd.DataFrame(records)

print(metadata.to_string(index=False))
print("\nSamples per condition:")
print(metadata["condition"].value_counts())

print("\nUnique patient IDs:", metadata["patient_id"].nunique())
pair_counts = (
    metadata.groupby(["patient_id", "condition"])
    .size()
    .unstack(fill_value=0)
)

print("\nSamples per patient and condition:")
print(pair_counts)

valid_pairs = (
    (pair_counts["Tumor"] == 1)
    & (pair_counts["Normal"] == 1)
)

print("\nNumber of patients with exactly one tumor and one normal:")
print(valid_pairs.sum())

print("\nAre all patient pairs valid?")
print(valid_pairs.all())
metadata.to_csv("data/gse251845_sample_metadata.csv", index=False)
print("\nMetadata saved to data/gse251845_sample_metadata.csv")
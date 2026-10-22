
from pathlib import Path

data_folder = Path("data/raw")

for file in data_folder.glob("*"):
    print(file.name, file.stat().st_size, "bytes")

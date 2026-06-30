import os
import requests

# Настройки
UNIPROT_ID = "Q8WUB8" # PHF10
RAW_DIR = "data/01_raw"
os.makedirs(RAW_DIR, exist_ok=True)

print(f"Загрузка референса {UNIPROT_ID}...")

# 1. Скачиваем FASTA
fasta_url = f"https://rest.uniprot.org/uniprotkb/{UNIPROT_ID}.fasta"
fasta_path = os.path.join(RAW_DIR, f"{UNIPROT_ID}.fasta")
with open(fasta_path, "w") as f:
    f.write(requests.get(fasta_url).text)
print(f"FASTA сохранена в {fasta_path}")

# 2. Скачиваем PDB (AlphaFold v4)
af_url = f"https://alphafold.ebi.ac.uk/files/AF-{UNIPROT_ID}-F1-model_v4.pdb"
af_path = os.path.join(RAW_DIR, f"AF-{UNIPROT_ID}.pdb")
with open(af_path, "w") as f:
    f.write(requests.get(af_url).text)
print(f"Структура AlphaFold сохранена в {af_path}")

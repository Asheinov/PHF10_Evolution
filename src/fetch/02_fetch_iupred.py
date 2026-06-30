import os
import json
import requests
from pathlib import Path

UNIPROT_ID = "Q8WUB8"
# Используем длинный тип предсказания (long), так как IDR обычно протяженные
IUPRED_URL = f"https://iupred3.elte.hu/iupred3/long/{UNIPROT_ID}.json"
OUTPUT_PATH = Path(
    f"data/raw/iupred3_{UNIPROT_ID}.json"
)

print(f"Запрашиваем данные IUPred3 для {UNIPROT_ID}...")

response = requests.get(
    IUPRED_URL,
    timeout=60
)

response.raise_for_status()
if response.status_code == 200:
    data = response.json()

    print(data.keys())
    print(type(data["iupred2"]))
    print(type(data["exp_dis"]))

    print("iupred2 length:", len(data["iupred2"]))
    print("exp_dis length:", len(data["exp_dis"]))

    print("iupred2 first 5:", data["iupred2"][:5])
    print("exp_dis first 5:", data["exp_dis"][:5])
    print("sequence length:", len(data["sequence"]))    
    with open(OUTPUT_PATH, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Успешно сохранено в {OUTPUT_PATH}")
else:
    print(f"Ошибка API: {response.status_code}")

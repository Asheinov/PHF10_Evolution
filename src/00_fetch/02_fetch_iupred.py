import os
import json
import requests

UNIPROT_ID = "Q8WUB8"
# Используем длинный тип предсказания (long), так как IDR обычно протяженные
IUPRED_URL = f"https://iupred3.elte.hu/iupred3/long/{UNIPROT_ID}.json"
OUTPUT_PATH = f"data/01_raw/iupred3_{UNIPROT_ID}.json"

print(f"Запрашиваем данные IUPred3 для {UNIPROT_ID}...")

response = requests.get(IUPRED_URL)
if response.status_code == 200:
    data = response.json()
    with open(OUTPUT_PATH, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Успешно сохранено в {OUTPUT_PATH}")
else:
    print(f"Ошибка API: {response.status_code}")

"""
Ortholog FASTA QC module.

Назначение:
    Проверка качества набора ортологов PHF10 перед MSA и филогенетикой.

Что делает модуль:
    - читает FASTA файл ортологов
    - извлекает метаданные (protein_id, species, gene_id)
    - считает длины последовательностей
    - определяет выбросы по длине
    - сохраняет QC-таблицу
    - печатает базовую статистику

Важно:
    НЕ фильтрует данные.
    Только анализирует и помечает.
"""

from __future__ import annotations

from pathlib import Path
import csv
from typing import List, Dict, Tuple


def parse_fasta(fasta_path: str) -> List[Dict]:
    """
    Парсинг FASTA с ожидаемым форматом заголовка:

    >protein_id|species|gene_id
    """

    records = []

    with open(fasta_path, "r", encoding="utf-8") as f:

        header = None
        sequence_parts = []

        for line in f:
            line = line.strip()

            if line.startswith(">"):

                if header is not None:
                    records.append(_build_record(header, sequence_parts))

                header = line[1:]
                sequence_parts = []

            else:
                sequence_parts.append(line)

        if header is not None:
            records.append(_build_record(header, sequence_parts))

    return records


def _build_record(header: str, seq_parts: List[str]) -> Dict:
    """
    Convert FASTA entry to structured dict.
    """

    parts = header.split("|")

    protein_id = parts[0]
    species = parts[1] if len(parts) > 1 else None
    gene_id = parts[2] if len(parts) > 2 else None

    sequence = "".join(seq_parts)

    return {
        "protein_id": protein_id,
        "species": species,
        "gene_id": gene_id,
        "sequence": sequence,
        "length": len(sequence),
    }


def compute_outliers(lengths: List[int]) -> Tuple[int, int]:
    """
    Простая эвристика выбросов:
    mean ± 2*std
    """

    import statistics

    mean = statistics.mean(lengths)
    stdev = statistics.stdev(lengths) if len(lengths) > 1 else 0

    lower = mean - 2 * stdev
    upper = mean + 2 * stdev

    return lower, upper


def run_qc(
    fasta_path: str,
    output_tsv: str
) -> None:
    """
    Main QC pipeline.
    """

    records = parse_fasta(fasta_path)

    lengths = [r["length"] for r in records]

    lower, upper = compute_outliers(lengths)

    print(f"Sequences: {len(records)}")
    print(f"Min length: {min(lengths)}")
    print(f"Max length: {max(lengths)}")
    print(f"Mean length: {sum(lengths)/len(lengths):.2f}")
    print(f"Outlier range: [{lower:.2f}, {upper:.2f}]")

    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=[
                "protein_id",
                "species",
                "gene_id",
                "length",
                "is_outlier",
            ],
            delimiter="\t",
        )

        writer.writeheader()

        for r in records:

            is_outlier = not (lower <= r["length"] <= upper)

            writer.writerow({
                "protein_id": r["protein_id"],
                "species": r["species"],
                "gene_id": r["gene_id"],
                "length": r["length"],
                "is_outlier": is_outlier,
            })


if __name__ == "__main__":

    run_qc(
        fasta_path="data/raw/orthology/phf10_orthologs.fasta",
        output_tsv="data/interim/orthology/ortholog_qc.tsv",
    )

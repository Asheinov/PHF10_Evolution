"""
Sort PHF10 orthologs by sequence length.

Задачи:
    - извлечь последовательности из FASTA
    - вычислить длину
    - нормализовать название вида
    - отсортировать по длине
    - сохранить результат

Без фильтрации, без интерпретации.
"""

from pathlib import Path
import csv


def normalize_species(species: str) -> str:
    """
    Упрощённое имя вида (без биоинформатической перегрузки).
    """
    if species is None:
        return "unknown"

    return species.replace("_", " ").capitalize()


def parse_fasta(path: str):
    """
    FASTA format:
    >protein_id|species|gene_id
    SEQ
    """

    records = []

    with open(path, "r", encoding="utf-8") as f:

        header = None
        seq = []

        for line in f:
            line = line.strip()

            if line.startswith(">"):

                if header:
                    records.append(build(header, seq))

                header = line[1:]
                seq = []

            else:
                seq.append(line)

        if header:
            records.append(build(header, seq))

    return records


def build(header: str, seq_parts):

    parts = header.split("|")

    seq = "".join(seq_parts)

    return {
        "protein_id": parts[0],
        "species": parts[1] if len(parts) > 1 else None,
        "gene_id": parts[2] if len(parts) > 2 else None,
        "length": len(seq),
    }


def run(input_fasta: str, output_tsv: str):

    records = parse_fasta(input_fasta)

    for r in records:
        r["species_simple"] = normalize_species(r["species"])

    records.sort(key=lambda x: x["length"], reverse=True)

    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)

    with open(output_tsv, "w", newline="", encoding="utf-8") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=[
                "protein_id",
                "species",
                "species_simple",
                "gene_id",
                "length",
            ],
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(records)

    print(f"Saved: {output_tsv}")
    print(f"Total: {len(records)}")


if __name__ == "__main__":

    run(
        input_fasta="data/raw/orthology/phf10_orthologs.fasta",
        output_tsv="data/interim/orthology/sorted_by_length.tsv",
    )

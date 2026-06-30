"""
Download protein sequences for orthologs from Ensembl.

Pipeline:
    ensembl.tsv
        ->
    unique protein IDs
        ->
    Ensembl sequence API
        ->
    FASTA file

This module performs ONLY sequence export.

No:
    - orthology inference
    - filtering
    - MSA
    - trait analysis

Input:
    data/raw/orthology/ensembl.tsv

Output:
    data/raw/orthology/phf10_orthologs.fasta
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict

import requests


def load_ensembl_table(tsv_path: str) -> list[dict]:
    """
    Load Ensembl ortholog table.

    Parameters
    ----------
    tsv_path : str
        Path to ensembl.tsv

    Returns
    -------
    list[dict]
        Rows from TSV file.
    """

    rows = []

    with open(tsv_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            rows.append(row)

    return rows


def extract_unique_proteins(rows: list[dict]) -> Dict[str, dict]:
    """
    Extract unique protein IDs.

    Returns
    -------
    dict
        protein_id -> row
    """

    proteins = {}

    for row in rows:

        protein_id = row["target_protein_id"]

        if not protein_id:
            continue

        if protein_id not in proteins:
            proteins[protein_id] = row

    return proteins


def fetch_protein_sequence(protein_id: str) -> str:
    """
    Download protein sequence from Ensembl.

    Parameters
    ----------
    protein_id : str

    Returns
    -------
    str
        Amino acid sequence.
    """

    url = f"https://rest.ensembl.org/sequence/id/{protein_id}"

    headers = {
        "Content-Type": "text/plain"
    }

    response = requests.get(url, headers=headers)

    if not response.ok:
        raise RuntimeError(
            f"Failed to fetch sequence for {protein_id}"
        )

    return response.text.strip()


def write_fasta(
    proteins: Dict[str, dict],
    output_fasta: str
) -> None:
    """
    Download all sequences and write FASTA.
    """

    output_path = Path(output_fasta)

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(output_path, "w", encoding="utf-8") as fasta:

        total = len(proteins)

        for index, (protein_id, row) in enumerate(
            proteins.items(),
            start=1
        ):

            print(
                f"[{index}/{total}] "
                f"{protein_id}"
            )

            sequence = fetch_protein_sequence(
                protein_id
            )

            header = (
                f">{protein_id}"
                f"|{row['target_species']}"
                f"|{row['target_gene']}"
            )

            fasta.write(header + "\n")
            fasta.write(sequence + "\n")


def run_export_fasta(
    input_tsv: str,
    output_fasta: str
) -> None:
    """
    Main pipeline.
    """

    rows = load_ensembl_table(
        input_tsv
    )

    proteins = extract_unique_proteins(
        rows
    )

    print(
        f"Loaded rows: {len(rows)}"
    )

    print(
        f"Unique proteins: {len(proteins)}"
    )

    write_fasta(
        proteins=proteins,
        output_fasta=output_fasta
    )


if __name__ == "__main__":

    run_export_fasta(
        input_tsv="data/raw/orthology/ensembl.tsv",
        output_fasta="data/raw/orthology/phf10_orthologs.fasta"
    )

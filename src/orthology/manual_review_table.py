"""
Create a manual review table for PHF10 ortholog candidates before MSA.

Purpose:
    Join existing length QC, architecture split, and domain QC outputs
    into one TSV for manual inspection.

Input files:
    data/interim/orthology/length_clusters.tsv
    data/interim/orthology/architecture_split.tsv
    data/interim/orthology/domain_scan/core_domain_qc_summary.tsv
    data/interim/orthology/domain_scan/extended_domain_qc_summary.tsv

Output files:
    data/interim/orthology/manual_review_table.tsv

Pipeline position:
    length_distribution / architecture_split / domain_qc_summary
        ->
    manual_review_table
        ->
    manual inspection before MSA

Notes:
    This module does not filter proteins, define orthology, or perform
    biological interpretation. It only joins existing QC outputs into
    one review table.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


LENGTH_CLUSTERS = Path("data/interim/orthology/length_clusters.tsv")
ARCHITECTURE_SPLIT = Path("data/interim/orthology/architecture_split.tsv")
CORE_DOMAIN_QC = Path(
    "data/interim/orthology/domain_scan/core_domain_qc_summary.tsv"
)
EXTENDED_DOMAIN_QC = Path(
    "data/interim/orthology/domain_scan/extended_domain_qc_summary.tsv"
)
OUTPUT_TABLE = Path("data/interim/orthology/manual_review_table.tsv")

OUTPUT_COLUMNS = [
    "protein_id",
    "species",
    "gene_id",
    "group",
    "protein_length",
    "length_cluster",
    "cluster_probability",
    "n_domain_hits",
    "n_unique_pfam_ids",
    "unique_pfam_ids",
    "unique_domain_names",
    "has_any_domain",
    "has_phd_domain",
    "has_ini1_dna_bd",
    "min_domain_start",
    "max_domain_end",
    "manual_status",
    "manual_note",
]


def load_tsv(path: str) -> list[dict[str, str]]:
    """
    Load a TSV file into a list of dictionaries.
    """

    input_path = Path(path)

    if not input_path.exists():
        raise FileNotFoundError(f"Required input file is missing: {input_path}")

    with input_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def index_by_protein_id(
    records: list[dict[str, str]]
) -> dict[str, dict[str, str]]:
    """
    Index records by protein_id.
    """

    indexed = {}

    for record in records:
        protein_id = record.get("protein_id", "")

        if not protein_id:
            raise ValueError(f"Record is missing required protein_id: {record}")

        indexed[protein_id] = record

    return indexed


def build_review_records(
    architecture_records: list[dict[str, str]],
    length_records: list[dict[str, str]],
    core_domain_records: list[dict[str, str]],
    extended_domain_records: list[dict[str, str]],
) -> list[dict[str, Any]]:
    """
    Build manual review rows using architecture_split as the full set.
    """

    length_by_protein = index_by_protein_id(length_records)
    domain_by_protein = index_by_protein_id(
        core_domain_records + extended_domain_records
    )
    review_records = []

    for architecture_record in architecture_records:
        protein_id = architecture_record.get("protein_id", "")

        if not protein_id:
            raise ValueError(
                "architecture_split record is missing required protein_id: "
                f"{architecture_record}"
            )

        length_record = length_by_protein.get(protein_id, {})
        domain_record = domain_by_protein.get(protein_id, {})

        review_records.append(
            {
                "protein_id": protein_id,
                "species": (
                    architecture_record.get("species")
                    or domain_record.get("species")
                    or length_record.get("species")
                    or ""
                ),
                "gene_id": (
                    architecture_record.get("gene_id")
                    or domain_record.get("gene_id")
                    or length_record.get("gene_id")
                    or ""
                ),
                "group": (
                    architecture_record.get("group")
                    or domain_record.get("group")
                    or ""
                ),
                "protein_length": (
                    domain_record.get("protein_length")
                    or architecture_record.get("length")
                    or length_record.get("length")
                    or ""
                ),
                "length_cluster": (
                    length_record.get("cluster")
                    or architecture_record.get("cluster")
                    or ""
                ),
                "cluster_probability": (
                    length_record.get("cluster_probability")
                    or architecture_record.get("cluster_probability")
                    or ""
                ),
                "n_domain_hits": domain_record.get("n_domain_hits", ""),
                "n_unique_pfam_ids": domain_record.get(
                    "n_unique_pfam_ids",
                    "",
                ),
                "unique_pfam_ids": domain_record.get("unique_pfam_ids", ""),
                "unique_domain_names": domain_record.get(
                    "unique_domain_names",
                    "",
                ),
                "has_any_domain": domain_record.get("has_any_domain", ""),
                "has_phd_domain": domain_record.get("has_phd_domain", ""),
                "has_ini1_dna_bd": domain_record.get("has_ini1_dna_bd", ""),
                "min_domain_start": domain_record.get("min_domain_start", ""),
                "max_domain_end": domain_record.get("max_domain_end", ""),
                "manual_status": "",
                "manual_note": "",
            }
        )

    return review_records


def write_tsv(records: list[dict[str, Any]], output_path: str) -> None:
    """
    Write manual review records to a TSV file.
    """

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
        )
        writer.writeheader()

        for record in records:
            writer.writerow(
                {
                    column: record.get(column, "")
                    for column in OUTPUT_COLUMNS
                }
            )


def run() -> None:
    """
    Create the manual review table.
    """

    length_records = load_tsv(str(LENGTH_CLUSTERS))
    architecture_records = load_tsv(str(ARCHITECTURE_SPLIT))
    core_domain_records = load_tsv(str(CORE_DOMAIN_QC))
    extended_domain_records = load_tsv(str(EXTENDED_DOMAIN_QC))

    review_records = build_review_records(
        architecture_records=architecture_records,
        length_records=length_records,
        core_domain_records=core_domain_records,
        extended_domain_records=extended_domain_records,
    )

    write_tsv(review_records, str(OUTPUT_TABLE))

    groups = sorted(
        {
            record.get("group", "")
            for record in review_records
            if record.get("group", "")
        }
    )
    print(
        f"Saved {len(review_records)} manual review rows to {OUTPUT_TABLE}"
    )
    print(f"Groups present: {', '.join(groups)}")


def main() -> int:
    """
    Command-line entry point with concise error reporting.
    """

    try:
        run()
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

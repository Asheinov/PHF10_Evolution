"""
Prepare a manual decision template for duplicated-species PHF10 records.

Purpose:
    Copy all records from duplicate-species QC into a manual review table
    and append empty decision fields for later human review.

Input files:
    data/interim/orthology/msa/duplicate_species_qc.tsv

Output files:
    data/interim/orthology/msa/duplicate_species_decisions.tsv

Pipeline position:
    duplicate_species_qc
        ->
    prepare_duplicate_species_decisions

Notes:
    This module prepares a manual decision template only. It does not
    filter records, redefine orthology, interpret biology, or use domains,
    IDR, Ser-rich sequence, N-terminal extension, PTM motifs, or trait-level
    information.
"""

from __future__ import annotations

import csv
from pathlib import Path


INPUT_QC_TSV = Path("data/interim/orthology/msa/duplicate_species_qc.tsv")
OUTPUT_DECISIONS_TSV = Path(
    "data/interim/orthology/msa/duplicate_species_decisions.tsv"
)

DEFAULT_DUPLICATE_DECISION = "REVIEW_DUPLICATE"

REQUIRED_FIELDS = {
    "species",
    "species_record_count",
    "duplicate_class",
    "protein_id",
    "gene_id",
    "raw_length",
    "tree_readiness",
    "nterm_readiness",
    "ptm_readiness",
    "first_aligned_human_position",
    "nterm_window_coverage",
    "ptm_window_coverage",
    "sequence_qc_flags",
}
OUTPUT_FIELDNAMES = [
    "species",
    "species_record_count",
    "duplicate_class",
    "protein_id",
    "gene_id",
    "raw_length",
    "tree_readiness",
    "nterm_readiness",
    "ptm_readiness",
    "first_aligned_human_position",
    "nterm_window_coverage",
    "ptm_window_coverage",
    "sequence_qc_flags",
    "duplicate_decision",
    "duplicate_note",
]


def require_input_file(path: Path) -> None:
    """
    Fail clearly when an expected input file is missing.
    """

    if not path.exists():
        raise FileNotFoundError(f"Required input file is missing: {path}")


def validate_fields(
    path: Path,
    fieldnames: list[str] | None,
    required_fields: set[str],
) -> None:
    """
    Check that a TSV header contains all fields required by this module.
    """

    if fieldnames is None:
        raise ValueError(f"Input TSV has no header: {path}")

    missing_fields = sorted(required_fields.difference(fieldnames))
    if missing_fields:
        fields = ", ".join(missing_fields)
        raise ValueError(f"Input TSV is missing required fields: {path}: {fields}")


def load_duplicate_qc(path: Path) -> list[dict[str, str]]:
    """
    Load duplicate-species QC rows from a tab-separated file.
    """

    require_input_file(path)
    records = []

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_fields(path, reader.fieldnames, REQUIRED_FIELDS)

        for row in reader:
            protein_id = row["protein_id"].strip()
            if not protein_id:
                raise ValueError(f"Duplicate QC row has empty protein_id: {path}")
            records.append(row)

    if not records:
        raise ValueError(f"Input TSV has no duplicate QC records: {path}")

    return records


def build_decision_records(
    duplicate_qc_records: list[dict[str, str]],
) -> list[dict[str, str]]:
    """
    Add default manual decision fields to duplicate QC records.
    """

    decision_records = []
    for record in duplicate_qc_records:
        decision_record = {
            fieldname: record.get(fieldname, "").strip()
            for fieldname in OUTPUT_FIELDNAMES
            if fieldname not in {"duplicate_decision", "duplicate_note"}
        }
        decision_record["duplicate_decision"] = DEFAULT_DUPLICATE_DECISION
        decision_record["duplicate_note"] = ""
        decision_records.append(decision_record)

    return decision_records


def write_tsv(records: list[dict[str, str]], output_path: Path) -> None:
    """
    Write the manual decision template as a tab-separated file.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_FIELDNAMES,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    """
    Create the duplicated-species manual decision template.
    """

    duplicate_qc_records = load_duplicate_qc(INPUT_QC_TSV)
    decision_records = build_decision_records(duplicate_qc_records)
    write_tsv(decision_records, OUTPUT_DECISIONS_TSV)

    print(f"duplicate records read: {len(duplicate_qc_records)}")
    print(f"decision records written: {len(decision_records)}")
    print(f"saved: {OUTPUT_DECISIONS_TSV}")


if __name__ == "__main__":
    main()

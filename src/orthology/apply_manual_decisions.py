"""
Apply manual review decisions to the PHF10 ortholog review table.

Purpose:
    Merge manual_review_table.tsv with manual_decisions.tsv to create a
    final manual review table for downstream pipeline steps.

Input files:
    data/interim/orthology/manual_review_table.tsv
    data/interim/orthology/manual_decisions.tsv

Output files:
    data/interim/orthology/manual_review_final.tsv

Pipeline position:
    manual_review_table
        ->
    apply_manual_decisions
        ->
    later manual-QC-aware steps

Notes:
    This module preserves all rows from the review table. It does not
    filter proteins, define orthology, or perform biological
    interpretation.
"""

from __future__ import annotations

import csv
from pathlib import Path


REVIEW_TABLE = Path("data/interim/orthology/manual_review_table.tsv")
MANUAL_DECISIONS = Path("data/interim/orthology/manual_decisions.tsv")
OUTPUT_TABLE = Path("data/interim/orthology/manual_review_final.tsv")


def load_tsv(path: str) -> tuple[list[dict[str, str]], list[str]]:
    """
    Load a TSV file and return records plus field names.
    """

    input_path = Path(path)

    if not input_path.exists():
        raise FileNotFoundError(f"Required input file is missing: {input_path}")

    with input_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Input TSV has no header: {input_path}")

        return [dict(row) for row in reader], list(reader.fieldnames)


def index_decisions(
    decisions: list[dict[str, str]]
) -> dict[str, dict[str, str]]:
    """
    Index manual decisions by protein_id.
    """

    indexed = {}

    for decision in decisions:
        protein_id = decision.get("protein_id", "")

        if not protein_id:
            raise ValueError(
                f"Manual decision is missing required protein_id: {decision}"
            )

        indexed[protein_id] = decision

    return indexed


def apply_decisions(
    review_records: list[dict[str, str]],
    decisions_by_protein: dict[str, dict[str, str]],
) -> tuple[list[dict[str, str]], int]:
    """
    Replace manual fields where a matching decision exists.
    """

    merged_records = []
    applied_count = 0

    for record in review_records:
        protein_id = record.get("protein_id", "")

        if not protein_id:
            raise ValueError(
                f"Review table row is missing required protein_id: {record}"
            )

        merged = dict(record)
        decision = decisions_by_protein.get(protein_id)

        if decision is None:
            merged["manual_status"] = ""
            merged["manual_note"] = ""
        else:
            merged["manual_status"] = decision.get("manual_status", "")
            merged["manual_note"] = decision.get("manual_note", "")
            applied_count += 1

        merged_records.append(merged)

    return merged_records, applied_count


def write_tsv(
    records: list[dict[str, str]],
    output_path: str,
    fieldnames: list[str],
) -> None:
    """
    Write merged manual review records to a TSV file.
    """

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()

        for record in records:
            writer.writerow(
                {
                    field: record.get(field, "")
                    for field in fieldnames
                }
            )


def run() -> None:
    """
    Apply manual decisions and write the final review table.
    """

    review_records, review_fieldnames = load_tsv(str(REVIEW_TABLE))
    decisions, _ = load_tsv(str(MANUAL_DECISIONS))
    decisions_by_protein = index_decisions(decisions)

    if "manual_status" not in review_fieldnames:
        review_fieldnames.append("manual_status")

    if "manual_note" not in review_fieldnames:
        review_fieldnames.append("manual_note")

    merged_records, applied_count = apply_decisions(
        review_records,
        decisions_by_protein,
    )

    write_tsv(
        records=merged_records,
        output_path=str(OUTPUT_TABLE),
        fieldnames=review_fieldnames,
    )

    print(f"Manual decisions applied: {applied_count}")
    print(f"Saved: {OUTPUT_TABLE}")


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

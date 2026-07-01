"""
Build strict PHF10 MSA input from manually reviewed ortholog QC results.

Purpose:
    Create the first strict FASTA and metadata TSV for MSA preparation.

Input files:
    data/interim/orthology/manual_review_final.tsv
    data/raw/orthology/phf10_orthologs.fasta

Output files:
    data/interim/orthology/msa_input_strict.fasta
    data/interim/orthology/msa_input_strict.tsv

Pipeline position:
    apply_manual_decisions
        ->
    build_msa_input
        ->
    MSA

Notes:
    This module filters only by manual_status. It does not filter by
    domains, length, IDR, Ser-rich regions, or N-terminal extensions,
    and it does not perform biological interpretation.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


MANUAL_REVIEW_FINAL = Path("data/interim/orthology/manual_review_final.tsv")
SOURCE_FASTA = Path("data/raw/orthology/phf10_orthologs.fasta")
OUTPUT_FASTA = Path("data/interim/orthology/msa_input_strict.fasta")
OUTPUT_TSV = Path("data/interim/orthology/msa_input_strict.tsv")

TSV_COLUMNS = [
    "protein_id",
    "species",
    "gene_id",
    "group",
    "protein_length",
    "manual_status",
    "manual_note",
]


def load_review_table(path: str) -> list[dict[str, str]]:
    """
    Load the final manual review table.
    """

    input_path = Path(path)

    if not input_path.exists():
        raise FileNotFoundError(f"Required input file is missing: {input_path}")

    with input_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Input TSV has no header: {input_path}")

        return [dict(row) for row in reader]


def parse_fasta(path: str) -> dict[str, dict[str, str]]:
    """
    Load FASTA records keyed by protein ID while preserving headers.
    """

    input_path = Path(path)

    if not input_path.exists():
        raise FileNotFoundError(f"Required FASTA input is missing: {input_path}")

    records: dict[str, dict[str, str]] = {}
    header: str | None = None
    sequence_parts: list[str] = []

    def store_record(record_header: str, sequence: list[str]) -> None:
        protein_id = record_header.split("|", maxsplit=1)[0]

        if not protein_id:
            raise ValueError(f"FASTA header has empty protein_id: {record_header}")

        if protein_id in records:
            raise ValueError(f"Duplicate FASTA protein_id: {protein_id}")

        records[protein_id] = {
            "header": record_header,
            "sequence": "".join(sequence),
        }

    with input_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()

            if not stripped:
                continue

            if stripped.startswith(">"):
                if header is not None:
                    store_record(header, sequence_parts)

                header = stripped[1:]
                sequence_parts = []
            else:
                if header is None:
                    raise ValueError(
                        f"Sequence encountered before FASTA header in {input_path}"
                    )

                sequence_parts.append(stripped)

    if header is not None:
        store_record(header, sequence_parts)

    return records


def select_included_rows(
    rows: list[dict[str, str]]
) -> tuple[list[dict[str, str]], dict[str, int]]:
    """
    Select rows with empty manual_status and count excluded statuses.
    """

    included = []
    summary = {
        "total_rows": len(rows),
        "included_rows": 0,
        "excluded_review_rows": 0,
        "excluded_exclude_rows": 0,
    }

    for row in rows:
        status = row.get("manual_status", "").strip().upper()

        if status == "":
            included.append(row)
            summary["included_rows"] += 1
        elif status == "REVIEW":
            summary["excluded_review_rows"] += 1
        elif status == "EXCLUDE":
            summary["excluded_exclude_rows"] += 1

    return included, summary


def require_fasta_records(
    included_rows: list[dict[str, str]],
    fasta_records: dict[str, dict[str, str]],
) -> None:
    """
    Fail clearly if any included protein is missing from the FASTA.
    """

    missing = [
        row.get("protein_id", "")
        for row in included_rows
        if row.get("protein_id", "") not in fasta_records
    ]

    if missing:
        raise ValueError(
            "Included protein_id values missing from FASTA: "
            + ", ".join(sorted(missing))
        )


def write_fasta(
    included_rows: list[dict[str, str]],
    fasta_records: dict[str, dict[str, str]],
    output_path: str,
) -> int:
    """
    Write included sequences to FASTA with original headers.
    """

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    written = 0

    with path.open("w", encoding="utf-8") as handle:
        for row in included_rows:
            protein_id = row["protein_id"]
            record = fasta_records[protein_id]
            handle.write(f">{record['header']}\n")
            handle.write(f"{record['sequence']}\n")
            written += 1

    return written


def write_tsv(
    included_rows: list[dict[str, str]],
    output_path: str,
) -> None:
    """
    Write metadata for included MSA input records.
    """

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TSV_COLUMNS,
            delimiter="\t",
        )
        writer.writeheader()

        for row in included_rows:
            writer.writerow(
                {
                    column: row.get(column, "")
                    for column in TSV_COLUMNS
                }
            )


def run() -> None:
    """
    Build strict MSA input FASTA and metadata TSV.
    """

    review_rows = load_review_table(str(MANUAL_REVIEW_FINAL))
    fasta_records = parse_fasta(str(SOURCE_FASTA))
    included_rows, summary = select_included_rows(review_rows)

    require_fasta_records(included_rows, fasta_records)

    written_count = write_fasta(
        included_rows=included_rows,
        fasta_records=fasta_records,
        output_path=str(OUTPUT_FASTA),
    )
    write_tsv(included_rows, str(OUTPUT_TSV))

    print(f"total rows: {summary['total_rows']}")
    print(f"included rows: {summary['included_rows']}")
    print(f"excluded REVIEW rows: {summary['excluded_review_rows']}")
    print(f"excluded EXCLUDE rows: {summary['excluded_exclude_rows']}")
    print(f"written FASTA records: {written_count}")


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

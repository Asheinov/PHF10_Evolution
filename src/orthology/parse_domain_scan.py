"""
Parse HMMER Pfam domain scan output for PHF10 ortholog QC.

Purpose:
    Convert HMMER --domtblout files from domain_scan into clean TSV
    tables for downstream quality-control steps.

Input files:
    data/interim/orthology/domain_scan/core.domtblout
    data/interim/orthology/domain_scan/extended.domtblout

Output files:
    data/interim/orthology/domain_scan/core_domains.tsv
    data/interim/orthology/domain_scan/extended_domains.tsv

Pipeline position:
    domain_scan
        ->
    parse_domain_scan
        ->
    orthology_qc

Notes:
    This module parses scan output only. It does not filter hits,
    collapse overlapping domains, or interpret PHD, SAY, INI1, or any
    other domain biology. Domain scanning is QC only and must not define
    PHF10 orthology.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


CORE_DOMTBLOUT = Path("data/interim/orthology/domain_scan/core.domtblout")
EXTENDED_DOMTBLOUT = Path(
    "data/interim/orthology/domain_scan/extended.domtblout"
)
CORE_OUTPUT = Path("data/interim/orthology/domain_scan/core_domains.tsv")
EXTENDED_OUTPUT = Path(
    "data/interim/orthology/domain_scan/extended_domains.tsv"
)

OUTPUT_COLUMNS = [
    "protein_id",
    "species",
    "gene_id",
    "protein_length",
    "domain_name",
    "pfam_id",
    "domain_length",
    "full_sequence_evalue",
    "full_sequence_score",
    "domain_index",
    "domain_count",
    "domain_conditional_evalue",
    "domain_independent_evalue",
    "domain_score",
    "hmm_start",
    "hmm_end",
    "ali_start",
    "ali_end",
    "env_start",
    "env_end",
    "accuracy",
    "description",
]


def parse_query_name(query_name: str) -> tuple[str, str, str]:
    """
    Split a PHF10 FASTA query name into protein, species, and gene IDs.
    """

    parts = query_name.split("|")

    if len(parts) != 3:
        raise ValueError(
            "Expected query_name format protein_id|species|gene_id, "
            f"got: {query_name}"
        )

    protein_id, species, gene_id = parts
    return protein_id, species, gene_id


def parse_domtblout_line(line: str) -> dict[str, Any]:
    """
    Parse one non-comment HMMER domtblout line.
    """

    parts = line.strip().split(maxsplit=22)

    if len(parts) < 22:
        raise ValueError(
            "Malformed domtblout line with fewer than 22 fields: "
            f"{line.rstrip()}"
        )

    description = parts[22] if len(parts) > 22 else ""
    protein_id, species, gene_id = parse_query_name(parts[3])

    return {
        "protein_id": protein_id,
        "species": species,
        "gene_id": gene_id,
        "protein_length": int(parts[5]),
        "domain_name": parts[0],
        "pfam_id": parts[1],
        "domain_length": int(parts[2]),
        "full_sequence_evalue": float(parts[6]),
        "full_sequence_score": float(parts[7]),
        "domain_index": int(parts[9]),
        "domain_count": int(parts[10]),
        "domain_conditional_evalue": float(parts[11]),
        "domain_independent_evalue": float(parts[12]),
        "domain_score": float(parts[13]),
        "hmm_start": int(parts[15]),
        "hmm_end": int(parts[16]),
        "ali_start": int(parts[17]),
        "ali_end": int(parts[18]),
        "env_start": int(parts[19]),
        "env_end": int(parts[20]),
        "accuracy": float(parts[21]),
        "description": description,
    }


def parse_domtblout(path: str) -> list[dict[str, Any]]:
    """
    Parse all non-comment HMMER domtblout records from a file.
    """

    input_path = Path(path)

    if not input_path.exists():
        raise FileNotFoundError(
            f"Required domtblout input is missing: {input_path}"
        )

    records = []

    with input_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            try:
                records.append(parse_domtblout_line(line))
            except ValueError as exc:
                raise ValueError(
                    f"Failed to parse {input_path} line {line_number}: {exc}"
                ) from exc

    return records


def write_tsv(records: list[dict[str, Any]], output_path: str) -> None:
    """
    Write parsed domain records to a TSV file.
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
                    column: record.get(column)
                    for column in OUTPUT_COLUMNS
                }
            )


def print_summary(label: str, records: list[dict[str, Any]]) -> None:
    """
    Print a short parsing summary for one architecture group.
    """

    unique_proteins = {
        record["protein_id"]
        for record in records
    }

    print(
        f"{label}: {len(records)} parsed records; "
        f"{len(unique_proteins)} proteins with at least one domain hit"
    )


def run() -> None:
    """
    Parse core and extended domtblout files and write TSV outputs.
    """

    core_records = parse_domtblout(str(CORE_DOMTBLOUT))
    extended_records = parse_domtblout(str(EXTENDED_DOMTBLOUT))

    write_tsv(core_records, str(CORE_OUTPUT))
    write_tsv(extended_records, str(EXTENDED_OUTPUT))

    print_summary("core", core_records)
    print_summary("extended", extended_records)


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

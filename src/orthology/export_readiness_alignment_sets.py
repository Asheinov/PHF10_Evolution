"""
Export reference-aware PHF10 aligned FASTA subsets for manual inspection.

Purpose:
    Write aligned FASTA subsets from the reference-aware MAFFT alignment
    for visual review in alignment viewers such as AliView or Jalview.

Input files:
    data/interim/orthology/msa/
    msa_input_strict_with_reference.mafft.fasta
    data/interim/orthology/msa/
    msa_input_strict_with_reference.sequence_readiness_qc.tsv

Output files:
    data/interim/orthology/msa/readiness_sets/all_with_reference.aligned.fasta
    data/interim/orthology/msa/readiness_sets/tree_ok_with_reference.aligned.fasta
    data/interim/orthology/msa/readiness_sets/nterm_ok_with_reference.aligned.fasta
    data/interim/orthology/msa/readiness_sets/ptm_ok_with_reference.aligned.fasta
    data/interim/orthology/msa/readiness_sets/nterm_truncated_with_reference.aligned.fasta
    data/interim/orthology/msa/readiness_sets/ptm_not_ready_with_reference.aligned.fasta
    data/interim/orthology/msa/readiness_sets/x_window_review_with_reference.aligned.fasta
    data/interim/orthology/msa/readiness_sets/multiple_records_species_with_reference.aligned.fasta
    data/interim/orthology/msa/readiness_sets/readiness_alignment_sets_summary.tsv

Pipeline position:
    run_msa_with_reference
        ->
    sequence_readiness_qc
        ->
    export_readiness_alignment_sets

Notes:
    These subsets are for manual visual inspection only. This module
    preserves original aligned FASTA headers and aligned sequences,
    does not trim or modify alignment columns, does not filter the
    original alignment, does not redefine orthology, and does not use
    Ser-rich sequence, N-terminal extension, phosphosite motifs,
    domains, or trait-level information.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


ALIGNMENT_FASTA = Path(
    "data/interim/orthology/msa/"
    "msa_input_strict_with_reference.mafft.fasta"
)
READINESS_TSV = Path(
    "data/interim/orthology/msa/"
    "msa_input_strict_with_reference.sequence_readiness_qc.tsv"
)
OUTPUT_DIR = Path("data/interim/orthology/msa/readiness_sets")
SUMMARY_TSV = OUTPUT_DIR / "readiness_alignment_sets_summary.tsv"
REFERENCE_PROTEIN_ID = "Q8WUB8"


@dataclass(frozen=True)
class AlignedRecord:
    """
    One aligned FASTA record.
    """

    protein_id: str
    header: str
    sequence: str


@dataclass(frozen=True)
class ExportSet:
    """
    One readiness-based export definition.
    """

    name: str
    output_path: Path
    selector: Callable[[dict[str, str]], bool]


@dataclass(frozen=True)
class ExportSummary:
    """
    Summary for one exported FASTA subset.
    """

    set_name: str
    output_fasta: str
    records_written: int
    non_reference_records_written: int
    reference_included: bool


def require_file(path: Path) -> None:
    """
    Fail clearly if a required input file is missing.
    """

    if path.exists():
        return

    raise FileNotFoundError(f"Required input file is missing: {path}")


def parse_aligned_fasta(path: Path) -> list[AlignedRecord]:
    """
    Read aligned FASTA records manually.
    """

    require_file(path)

    records: list[AlignedRecord] = []
    seen_protein_ids: set[str] = set()
    header: str | None = None
    sequence_parts: list[str] = []

    def store_record(record_header: str, sequence: list[str]) -> None:
        protein_id = record_header.split("|", maxsplit=1)[0]

        if not protein_id:
            raise ValueError(f"FASTA header has empty protein_id: {record_header}")

        if protein_id in seen_protein_ids:
            raise ValueError(
                f"Duplicate protein_id found in alignment: {protein_id}"
            )

        seen_protein_ids.add(protein_id)
        records.append(
            AlignedRecord(
                protein_id=protein_id,
                header=record_header,
                sequence="".join(sequence),
            )
        )

    with path.open("r", encoding="utf-8") as handle:
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
                        f"Sequence encountered before FASTA header in {path}"
                    )

                sequence_parts.append(stripped)

    if header is not None:
        store_record(header, sequence_parts)

    if not records:
        raise ValueError(f"Aligned FASTA is empty: {path}")

    return records


def load_readiness_table(path: Path) -> list[dict[str, str]]:
    """
    Load sequence readiness QC rows from TSV.
    """

    require_file(path)

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Input TSV has no header: {path}")

        return [dict(row) for row in reader]


def index_alignment(
    records: list[AlignedRecord],
) -> dict[str, AlignedRecord]:
    """
    Index aligned records by protein_id.
    """

    return {
        record.protein_id: record
        for record in records
    }


def define_export_sets() -> list[ExportSet]:
    """
    Define all readiness-based aligned FASTA exports.
    """

    return [
        ExportSet(
            name="all_with_reference",
            output_path=OUTPUT_DIR / "all_with_reference.aligned.fasta",
            selector=lambda row: True,
        ),
        ExportSet(
            name="tree_ok_with_reference",
            output_path=OUTPUT_DIR / "tree_ok_with_reference.aligned.fasta",
            selector=lambda row: row.get("tree_readiness") == "TREE_OK",
        ),
        ExportSet(
            name="nterm_ok_with_reference",
            output_path=OUTPUT_DIR / "nterm_ok_with_reference.aligned.fasta",
            selector=lambda row: row.get("nterm_readiness") == "NTERM_OK",
        ),
        ExportSet(
            name="ptm_ok_with_reference",
            output_path=OUTPUT_DIR / "ptm_ok_with_reference.aligned.fasta",
            selector=lambda row: row.get("ptm_readiness") == "PTM_OK",
        ),
        ExportSet(
            name="nterm_truncated_with_reference",
            output_path=(
                OUTPUT_DIR / "nterm_truncated_with_reference.aligned.fasta"
            ),
            selector=lambda row: (
                row.get("nterm_readiness") == "NTERM_TRUNCATED"
            ),
        ),
        ExportSet(
            name="ptm_not_ready_with_reference",
            output_path=(
                OUTPUT_DIR / "ptm_not_ready_with_reference.aligned.fasta"
            ),
            selector=lambda row: row.get("ptm_readiness") == "PTM_NOT_READY",
        ),
        ExportSet(
            name="x_window_review_with_reference",
            output_path=OUTPUT_DIR / "x_window_review_with_reference.aligned.fasta",
            selector=lambda row: (
                int(row.get("x_count_nterm_window", "0")) > 0
                or int(row.get("x_count_ptm_window", "0")) > 0
            ),
        ),
        ExportSet(
            name="multiple_records_species_with_reference",
            output_path=(
                OUTPUT_DIR
                / "multiple_records_species_with_reference.aligned.fasta"
            ),
            selector=lambda row: row.get("is_multi_record_species") == "True",
        ),
    ]


def select_protein_ids(
    readiness_rows: list[dict[str, str]],
    export_set: ExportSet,
) -> list[str]:
    """
    Select protein IDs for one export set and always include reference.
    """

    selected_ids = [
        row["protein_id"]
        for row in readiness_rows
        if export_set.selector(row)
    ]

    if REFERENCE_PROTEIN_ID not in selected_ids:
        selected_ids.append(REFERENCE_PROTEIN_ID)

    seen = set()
    unique_ids = []

    for protein_id in selected_ids:
        if protein_id not in seen:
            seen.add(protein_id)
            unique_ids.append(protein_id)

    return unique_ids


def validate_selected_ids(
    selected_ids: list[str],
    alignment_by_protein: dict[str, AlignedRecord],
) -> None:
    """
    Fail clearly if selected records are missing from the alignment.
    """

    if REFERENCE_PROTEIN_ID not in alignment_by_protein:
        raise ValueError(
            f"Reference record is missing from alignment: {REFERENCE_PROTEIN_ID}"
        )

    missing = [
        protein_id
        for protein_id in selected_ids
        if protein_id not in alignment_by_protein
    ]

    if missing:
        raise ValueError(
            "Selected protein_id values missing from alignment: "
            + ", ".join(sorted(missing))
        )


def write_fasta(
    selected_ids: list[str],
    alignment_by_protein: dict[str, AlignedRecord],
    output_path: Path,
) -> None:
    """
    Write selected aligned records, preserving headers and columns.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        for protein_id in selected_ids:
            record = alignment_by_protein[protein_id]
            handle.write(f">{record.header}\n")
            handle.write(f"{record.sequence}\n")


def export_set(
    readiness_rows: list[dict[str, str]],
    alignment_by_protein: dict[str, AlignedRecord],
    export_definition: ExportSet,
) -> ExportSummary:
    """
    Export one aligned FASTA subset and return summary metadata.
    """

    selected_ids = select_protein_ids(readiness_rows, export_definition)
    validate_selected_ids(selected_ids, alignment_by_protein)
    write_fasta(
        selected_ids=selected_ids,
        alignment_by_protein=alignment_by_protein,
        output_path=export_definition.output_path,
    )

    reference_included = REFERENCE_PROTEIN_ID in selected_ids
    records_written = len(selected_ids)
    non_reference_records_written = records_written - int(reference_included)

    print(f"saved: {export_definition.output_path} records: {records_written}")

    return ExportSummary(
        set_name=export_definition.name,
        output_fasta=str(export_definition.output_path),
        records_written=records_written,
        non_reference_records_written=non_reference_records_written,
        reference_included=reference_included,
    )


def write_summary(summaries: list[ExportSummary], output_path: Path) -> None:
    """
    Write TSV summary for all exported sets.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "set_name",
                "output_fasta",
                "records_written",
                "non_reference_records_written",
                "reference_included",
            ],
            delimiter="\t",
        )
        writer.writeheader()

        for summary in summaries:
            writer.writerow(
                {
                    "set_name": summary.set_name,
                    "output_fasta": summary.output_fasta,
                    "records_written": summary.records_written,
                    "non_reference_records_written": (
                        summary.non_reference_records_written
                    ),
                    "reference_included": summary.reference_included,
                }
            )


def run() -> None:
    """
    Export all readiness alignment sets.
    """

    aligned_records = parse_aligned_fasta(ALIGNMENT_FASTA)
    readiness_rows = load_readiness_table(READINESS_TSV)
    alignment_by_protein = index_alignment(aligned_records)
    validate_selected_ids([REFERENCE_PROTEIN_ID], alignment_by_protein)

    summaries = [
        export_set(
            readiness_rows=readiness_rows,
            alignment_by_protein=alignment_by_protein,
            export_definition=export_definition,
        )
        for export_definition in define_export_sets()
    ]
    write_summary(summaries, SUMMARY_TSV)
    print(f"saved: {SUMMARY_TSV} records: {len(summaries)}")


def main() -> None:
    """
    Command-line entry point with concise error reporting.
    """

    try:
        run()
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}")
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()

"""
Export N-terminal form aligned FASTA subsets for manual inspection.

Purpose:
    Write aligned FASTA subsets from the reference-aware PHF10 MAFFT
    alignment using existing N-terminal form QC annotations. These subsets
    are intended for manual visual inspection in alignment viewers such as
    AliView or Jalview.

Input files:
    data/interim/orthology/msa/
    msa_input_strict_with_reference.mafft.fasta
    data/interim/orthology/msa/
    msa_input_strict_with_reference.nterm_form_qc.tsv

Output files:
    data/interim/orthology/msa/nterm_form_sets/
    all_with_reference.aligned.fasta
    full_nterm_extension_with_reference.aligned.fasta
    mgs_start_form_with_reference.aligned.fasta
    downstream_start_form_with_reference.aligned.fasta
    nterm_divergent_with_reference.aligned.fasta
    nterm_not_assessable_with_reference.aligned.fasta
    s50_substituted_with_reference.aligned.fasta
    s50_gap_with_reference.aligned.fasta
    s50_x_with_reference.aligned.fasta
    ptm_functional_ready_with_reference.aligned.fasta
    ptm_sequence_present_context_missing_with_reference.aligned.fasta
    ptm_context_review_with_reference.aligned.fasta
    ptm_not_ready_with_reference.aligned.fasta
    nterm_form_alignment_sets_summary.tsv

Pipeline position:
    run_msa_with_reference
        ->
    nterm_form_qc
        ->
    export_nterm_form_alignment_sets

Notes:
    This module exports inspection subsets only. It does not filter the
    main dataset, redefine orthology, infer evolutionary gain or loss,
    classify biological isoforms, trim alignment columns, or modify aligned
    sequences.
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
NTERM_FORM_QC_TSV = Path(
    "data/interim/orthology/msa/"
    "msa_input_strict_with_reference.nterm_form_qc.tsv"
)
OUTPUT_DIR = Path("data/interim/orthology/msa/nterm_form_sets")
SUMMARY_TSV = OUTPUT_DIR / "nterm_form_alignment_sets_summary.tsv"
REFERENCE_PROTEIN_ID = "Q8WUB8"
FASTA_WRAP_WIDTH = 80

REQUIRED_QC_FIELDS = {
    "protein_id",
    "nterm_architecture_class",
    "s50_equivalent_status",
    "ptm_context_class",
}
SUMMARY_FIELDNAMES = [
    "subset_name",
    "output_fasta",
    "non_reference_records",
    "reference_records",
    "total_records",
]


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
    One N-terminal form export definition.
    """

    name: str
    output_path: Path
    selector: Callable[[dict[str, str]], bool]


@dataclass(frozen=True)
class ExportSummary:
    """
    Summary metadata for one exported aligned FASTA subset.
    """

    subset_name: str
    output_fasta: str
    non_reference_records: int
    reference_records: int
    total_records: int


def require_file(path: Path) -> None:
    """
    Fail clearly if a required input file is missing.
    """

    if not path.exists():
        raise FileNotFoundError(f"Required input file is missing: {path}")


def parse_aligned_fasta(path: Path) -> list[AlignedRecord]:
    """
    Parse aligned FASTA records manually.
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


def load_nterm_form_qc(path: Path) -> list[dict[str, str]]:
    """
    Load N-terminal form QC rows from TSV.
    """

    require_file(path)

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_fields(path, reader.fieldnames, REQUIRED_QC_FIELDS)
        rows = [dict(row) for row in reader]

    if not rows:
        raise ValueError(f"Input TSV has no data records: {path}")

    protein_ids = [row["protein_id"].strip() for row in rows]
    if any(not protein_id for protein_id in protein_ids):
        raise ValueError(f"Input TSV contains an empty protein_id: {path}")

    duplicated_ids = sorted(
        {
            protein_id
            for protein_id in protein_ids
            if protein_ids.count(protein_id) > 1
        }
    )
    if duplicated_ids:
        raise ValueError(
            "Duplicate protein_id values found in N-terminal form QC: "
            + ", ".join(duplicated_ids)
        )

    reference_count = protein_ids.count(REFERENCE_PROTEIN_ID)
    if reference_count == 0:
        raise ValueError(
            f"Reference record is missing from N-terminal form QC: "
            f"{REFERENCE_PROTEIN_ID}"
        )
    if reference_count > 1:
        raise ValueError(
            f"Reference record appears more than once in N-terminal form QC: "
            f"{REFERENCE_PROTEIN_ID}"
        )

    return rows


def index_alignment(records: list[AlignedRecord]) -> dict[str, AlignedRecord]:
    """
    Index aligned records by protein_id.
    """

    alignment_by_protein = {
        record.protein_id: record
        for record in records
    }

    reference_count = sum(
        record.protein_id == REFERENCE_PROTEIN_ID
        for record in records
    )
    if reference_count == 0:
        raise ValueError(
            f"Reference record is missing from alignment: {REFERENCE_PROTEIN_ID}"
        )
    if reference_count > 1:
        raise ValueError(
            f"Reference record appears more than once in alignment: "
            f"{REFERENCE_PROTEIN_ID}"
        )

    return alignment_by_protein


def validate_qc_records_in_alignment(
    qc_rows: list[dict[str, str]],
    alignment_by_protein: dict[str, AlignedRecord],
) -> None:
    """
    Fail clearly if any QC protein_id is absent from the alignment.
    """

    missing = [
        row["protein_id"].strip()
        for row in qc_rows
        if row["protein_id"].strip() not in alignment_by_protein
    ]

    if missing:
        raise ValueError(
            "N-terminal form QC protein_id values missing from alignment: "
            + ", ".join(sorted(missing))
        )


def define_export_sets() -> list[ExportSet]:
    """
    Define all N-terminal form aligned FASTA exports.
    """

    return [
        ExportSet(
            name="all_with_reference",
            output_path=OUTPUT_DIR / "all_with_reference.aligned.fasta",
            selector=lambda row: True,
        ),
        ExportSet(
            name="full_nterm_extension_with_reference",
            output_path=(
                OUTPUT_DIR / "full_nterm_extension_with_reference.aligned.fasta"
            ),
            selector=lambda row: (
                row.get("nterm_architecture_class") == "FULL_NTERM_EXTENSION"
            ),
        ),
        ExportSet(
            name="mgs_start_form_with_reference",
            output_path=OUTPUT_DIR / "mgs_start_form_with_reference.aligned.fasta",
            selector=lambda row: (
                row.get("nterm_architecture_class") == "MGS_START_FORM"
            ),
        ),
        ExportSet(
            name="downstream_start_form_with_reference",
            output_path=(
                OUTPUT_DIR / "downstream_start_form_with_reference.aligned.fasta"
            ),
            selector=lambda row: (
                row.get("nterm_architecture_class") == "DOWNSTREAM_START_FORM"
            ),
        ),
        ExportSet(
            name="nterm_divergent_with_reference",
            output_path=(
                OUTPUT_DIR / "nterm_divergent_with_reference.aligned.fasta"
            ),
            selector=lambda row: (
                row.get("nterm_architecture_class") == "NTERM_DIVERGENT"
            ),
        ),
        ExportSet(
            name="nterm_not_assessable_with_reference",
            output_path=(
                OUTPUT_DIR / "nterm_not_assessable_with_reference.aligned.fasta"
            ),
            selector=lambda row: (
                row.get("nterm_architecture_class") == "NTERM_NOT_ASSESSABLE"
            ),
        ),
        ExportSet(
            name="s50_substituted_with_reference",
            output_path=OUTPUT_DIR / "s50_substituted_with_reference.aligned.fasta",
            selector=lambda row: (
                row.get("s50_equivalent_status") == "S50_SUBSTITUTED"
            ),
        ),
        ExportSet(
            name="s50_gap_with_reference",
            output_path=OUTPUT_DIR / "s50_gap_with_reference.aligned.fasta",
            selector=lambda row: row.get("s50_equivalent_status") == "S50_GAP",
        ),
        ExportSet(
            name="s50_x_with_reference",
            output_path=OUTPUT_DIR / "s50_x_with_reference.aligned.fasta",
            selector=lambda row: row.get("s50_equivalent_status") == "S50_X",
        ),
        ExportSet(
            name="ptm_functional_ready_with_reference",
            output_path=(
                OUTPUT_DIR / "ptm_functional_ready_with_reference.aligned.fasta"
            ),
            selector=lambda row: (
                row.get("ptm_context_class") == "PTM_FUNCTIONAL_READY"
            ),
        ),
        ExportSet(
            name="ptm_sequence_present_context_missing_with_reference",
            output_path=(
                OUTPUT_DIR
                / "ptm_sequence_present_context_missing_with_reference.aligned.fasta"
            ),
            selector=lambda row: (
                row.get("ptm_context_class")
                == "PTM_SEQUENCE_PRESENT_CONTEXT_MISSING"
            ),
        ),
        ExportSet(
            name="ptm_context_review_with_reference",
            output_path=(
                OUTPUT_DIR / "ptm_context_review_with_reference.aligned.fasta"
            ),
            selector=lambda row: (
                row.get("ptm_context_class") == "PTM_CONTEXT_REVIEW"
            ),
        ),
        ExportSet(
            name="ptm_not_ready_with_reference",
            output_path=OUTPUT_DIR / "ptm_not_ready_with_reference.aligned.fasta",
            selector=lambda row: row.get("ptm_context_class") == "PTM_NOT_READY",
        ),
    ]


def select_protein_ids(
    qc_rows: list[dict[str, str]],
    aligned_records: list[AlignedRecord],
    export_set: ExportSet,
) -> list[str]:
    """
    Select protein IDs for one export set with Q8WUB8 first.
    """

    selected_ids = {
        row["protein_id"].strip()
        for row in qc_rows
        if export_set.selector(row)
    }
    selected_ids.add(REFERENCE_PROTEIN_ID)

    ordered_ids = [REFERENCE_PROTEIN_ID]
    for record in aligned_records:
        if record.protein_id == REFERENCE_PROTEIN_ID:
            continue
        if record.protein_id in selected_ids:
            ordered_ids.append(record.protein_id)

    if len(ordered_ids) != len(set(ordered_ids)):
        raise ValueError(
            f"Output subset contains duplicate protein_id values: "
            f"{export_set.name}"
        )

    return ordered_ids


def write_wrapped_sequence(
    sequence: str,
    handle,
    width: int = FASTA_WRAP_WIDTH,
) -> None:
    """
    Write one aligned sequence wrapped to the requested width.
    """

    for start in range(0, len(sequence), width):
        handle.write(sequence[start:start + width] + "\n")


def write_fasta(
    selected_ids: list[str],
    alignment_by_protein: dict[str, AlignedRecord],
    output_path: Path,
) -> None:
    """
    Write selected aligned records, preserving headers and aligned columns.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        for protein_id in selected_ids:
            record = alignment_by_protein[protein_id]
            handle.write(f">{record.header}\n")
            write_wrapped_sequence(record.sequence, handle)


def export_alignment_set(
    qc_rows: list[dict[str, str]],
    aligned_records: list[AlignedRecord],
    alignment_by_protein: dict[str, AlignedRecord],
    export_set: ExportSet,
) -> ExportSummary:
    """
    Export one aligned FASTA subset and return summary metadata.
    """

    selected_ids = select_protein_ids(
        qc_rows=qc_rows,
        aligned_records=aligned_records,
        export_set=export_set,
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

    write_fasta(
        selected_ids=selected_ids,
        alignment_by_protein=alignment_by_protein,
        output_path=export_set.output_path,
    )

    reference_records = sum(
        protein_id == REFERENCE_PROTEIN_ID
        for protein_id in selected_ids
    )
    if reference_records != 1:
        raise ValueError(
            f"Expected one reference record in subset {export_set.name}; "
            f"found {reference_records}"
        )

    total_records = len(selected_ids)
    non_reference_records = total_records - reference_records

    print(f"saved: {export_set.output_path} records: {total_records}")

    return ExportSummary(
        subset_name=export_set.name,
        output_fasta=str(export_set.output_path),
        non_reference_records=non_reference_records,
        reference_records=reference_records,
        total_records=total_records,
    )


def write_summary(summaries: list[ExportSummary], output_path: Path) -> None:
    """
    Write TSV summary for all exported N-terminal form alignment sets.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=SUMMARY_FIELDNAMES,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for summary in summaries:
            writer.writerow(
                {
                    "subset_name": summary.subset_name,
                    "output_fasta": summary.output_fasta,
                    "non_reference_records": summary.non_reference_records,
                    "reference_records": summary.reference_records,
                    "total_records": summary.total_records,
                }
            )


def run() -> None:
    """
    Export all N-terminal form alignment subsets.
    """

    aligned_records = parse_aligned_fasta(ALIGNMENT_FASTA)
    qc_rows = load_nterm_form_qc(NTERM_FORM_QC_TSV)
    alignment_by_protein = index_alignment(aligned_records)
    validate_qc_records_in_alignment(qc_rows, alignment_by_protein)

    summaries = [
        export_alignment_set(
            qc_rows=qc_rows,
            aligned_records=aligned_records,
            alignment_by_protein=alignment_by_protein,
            export_set=export_set,
        )
        for export_set in define_export_sets()
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

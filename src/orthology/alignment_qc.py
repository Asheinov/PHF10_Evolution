"""
Summarize technical QC metrics for the strict PHF10 MAFFT alignment.

Purpose:
    Read the first strict PHF10 ortholog MAFFT alignment and summarize
    per-sequence gap metrics for technical QC.

Input files:
    data/interim/orthology/msa/msa_input_strict.mafft.fasta
    data/interim/orthology/msa_input_strict.tsv

Output files:
    data/interim/orthology/msa/msa_input_strict.alignment_qc.tsv
    data/interim/orthology/msa/msa_input_strict.alignment_qc_summary.txt

Pipeline position:
    build_msa_input
        ->
    run_msa
        ->
    alignment_qc

Notes:
    This module is technical QC only. It does not filter sequences,
    modify the alignment, make inclusion/exclusion decisions, use
    domains, IDR, Ser-rich regions, N-terminal extensions, or perform
    biological interpretation.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any


ALIGNMENT_FASTA = Path(
    "data/interim/orthology/msa/msa_input_strict.mafft.fasta"
)
METADATA_TSV = Path("data/interim/orthology/msa_input_strict.tsv")
QC_TSV = Path("data/interim/orthology/msa/msa_input_strict.alignment_qc.tsv")
SUMMARY_TXT = Path(
    "data/interim/orthology/msa/msa_input_strict.alignment_qc_summary.txt"
)

QC_COLUMNS = [
    "protein_id",
    "species",
    "gene_id",
    "group",
    "raw_length",
    "aligned_length",
    "non_gap_length",
    "gap_count",
    "gap_fraction",
    "leading_gap_count",
    "trailing_gap_count",
    "internal_gap_count",
    "internal_gap_fraction",
    "manual_status",
    "manual_note",
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
class AlignmentQcRecord:
    """
    Technical QC metrics for one aligned sequence.
    """

    protein_id: str
    species: str
    gene_id: str
    group: str
    raw_length: int
    aligned_length: int
    non_gap_length: int
    gap_count: int
    gap_fraction: float
    leading_gap_count: int
    trailing_gap_count: int
    internal_gap_count: int
    internal_gap_fraction: float
    manual_status: str
    manual_note: str


def require_file(path: Path) -> None:
    """
    Fail clearly if a required input file is missing.
    """

    if path.exists():
        return

    raise FileNotFoundError(f"Required input file is missing: {path}")


def parse_aligned_fasta(path: Path) -> list[AlignedRecord]:
    """
    Read an aligned FASTA file manually.
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


def load_metadata(path: Path) -> dict[str, dict[str, str]]:
    """
    Load strict MSA metadata keyed by protein_id.
    """

    require_file(path)

    records: dict[str, dict[str, str]] = {}

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Metadata TSV has no header: {path}")

        for row in reader:
            protein_id = row.get("protein_id", "")

            if not protein_id:
                raise ValueError(f"Metadata row is missing protein_id: {row}")

            if protein_id in records:
                raise ValueError(
                    f"Duplicate protein_id found in metadata: {protein_id}"
                )

            records[protein_id] = dict(row)

    return records


def validate_aligned_lengths(records: list[AlignedRecord]) -> int:
    """
    Ensure all aligned sequences have the same length.
    """

    lengths = {len(record.sequence) for record in records}

    if len(lengths) != 1:
        raise ValueError(
            "Aligned sequence lengths differ between records: "
            + ", ".join(str(length) for length in sorted(lengths))
        )

    return lengths.pop()


def count_terminal_gaps(sequence: str) -> tuple[int, int]:
    """
    Count leading and trailing gap characters.
    """

    first_non_gap = None
    last_non_gap = None

    for index, character in enumerate(sequence):
        if character != "-":
            first_non_gap = index
            break

    for index in range(len(sequence) - 1, -1, -1):
        if sequence[index] != "-":
            last_non_gap = index
            break

    if first_non_gap is None or last_non_gap is None:
        raise ValueError("Aligned sequence contains no non-gap characters")

    leading_gap_count = first_non_gap
    trailing_gap_count = len(sequence) - last_non_gap - 1

    return leading_gap_count, trailing_gap_count


def build_qc_record(
    aligned_record: AlignedRecord,
    metadata: dict[str, str],
) -> AlignmentQcRecord:
    """
    Calculate technical QC metrics for one aligned sequence.
    """

    raw_length_text = metadata.get("protein_length", "")

    if not raw_length_text:
        raise ValueError(
            f"Missing protein_length in metadata for {aligned_record.protein_id}"
        )

    raw_length = int(raw_length_text)
    aligned_length = len(aligned_record.sequence)
    non_gap_length = sum(1 for char in aligned_record.sequence if char != "-")

    if non_gap_length == 0:
        raise ValueError(
            f"Aligned sequence contains no non-gap characters: "
            f"{aligned_record.protein_id}"
        )

    gap_count = aligned_record.sequence.count("-")
    leading_gap_count, trailing_gap_count = count_terminal_gaps(
        aligned_record.sequence
    )
    internal_gap_count = gap_count - leading_gap_count - trailing_gap_count

    return AlignmentQcRecord(
        protein_id=aligned_record.protein_id,
        species=metadata.get("species", ""),
        gene_id=metadata.get("gene_id", ""),
        group=metadata.get("group", ""),
        raw_length=raw_length,
        aligned_length=aligned_length,
        non_gap_length=non_gap_length,
        gap_count=gap_count,
        gap_fraction=gap_count / aligned_length,
        leading_gap_count=leading_gap_count,
        trailing_gap_count=trailing_gap_count,
        internal_gap_count=internal_gap_count,
        internal_gap_fraction=internal_gap_count / aligned_length,
        manual_status=metadata.get("manual_status", ""),
        manual_note=metadata.get("manual_note", ""),
    )


def build_qc_records(
    aligned_records: list[AlignedRecord],
    metadata_by_protein: dict[str, dict[str, str]],
) -> list[AlignmentQcRecord]:
    """
    Build QC metrics for all aligned records.
    """

    qc_records = []

    for aligned_record in aligned_records:
        metadata = metadata_by_protein.get(aligned_record.protein_id)

        if metadata is None:
            raise ValueError(
                "Aligned protein_id is missing from metadata: "
                f"{aligned_record.protein_id}"
            )

        qc_records.append(build_qc_record(aligned_record, metadata))

    return qc_records


def format_fraction(value: float) -> str:
    """
    Format fraction values consistently for TSV and summary outputs.
    """

    return f"{value:.6f}"


def qc_record_to_row(record: AlignmentQcRecord) -> dict[str, Any]:
    """
    Convert a QC dataclass record to a TSV row.
    """

    return {
        "protein_id": record.protein_id,
        "species": record.species,
        "gene_id": record.gene_id,
        "group": record.group,
        "raw_length": record.raw_length,
        "aligned_length": record.aligned_length,
        "non_gap_length": record.non_gap_length,
        "gap_count": record.gap_count,
        "gap_fraction": format_fraction(record.gap_fraction),
        "leading_gap_count": record.leading_gap_count,
        "trailing_gap_count": record.trailing_gap_count,
        "internal_gap_count": record.internal_gap_count,
        "internal_gap_fraction": format_fraction(record.internal_gap_fraction),
        "manual_status": record.manual_status,
        "manual_note": record.manual_note,
    }


def write_qc_tsv(records: list[AlignmentQcRecord], output_path: Path) -> None:
    """
    Write one alignment QC row per aligned FASTA record.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=QC_COLUMNS,
            delimiter="\t",
        )
        writer.writeheader()

        for record in records:
            writer.writerow(qc_record_to_row(record))


def mean(values: list[float]) -> float:
    """
    Calculate the arithmetic mean.
    """

    if not values:
        raise ValueError("Cannot calculate mean for an empty list")

    return sum(values) / len(values)


def get_length_mismatches(
    records: list[AlignmentQcRecord],
) -> list[AlignmentQcRecord]:
    """
    Return records where non-gap alignment length differs from raw length.
    """

    return [
        record
        for record in records
        if record.non_gap_length != record.raw_length
    ]


def format_top_record(record: AlignmentQcRecord) -> str:
    """
    Format one top-list summary row.
    """

    return (
        f"{record.protein_id} {record.species} {record.raw_length} "
        f"{record.aligned_length} {format_fraction(record.gap_fraction)} "
        f"{format_fraction(record.internal_gap_fraction)}"
    )


def write_summary(
    records: list[AlignmentQcRecord],
    alignment_length: int,
    output_path: Path,
) -> None:
    """
    Write a text summary for alignment technical QC.
    """

    gap_fractions = [record.gap_fraction for record in records]
    internal_gap_fractions = [
        record.internal_gap_fraction
        for record in records
    ]
    mismatches = get_length_mismatches(records)
    top_by_gap = sorted(
        records,
        key=lambda record: (-record.gap_fraction, record.protein_id),
    )[:20]
    top_by_internal_gap = sorted(
        records,
        key=lambda record: (-record.internal_gap_fraction, record.protein_id),
    )[:20]

    lines = [
        "PHF10 strict MSA alignment QC",
        "QC only: no filtering or biological interpretation.",
        "",
        f"number of sequences: {len(records)}",
        f"alignment length: {alignment_length}",
        f"min gap_fraction: {format_fraction(min(gap_fractions))}",
        f"max gap_fraction: {format_fraction(max(gap_fractions))}",
        f"mean gap_fraction: {format_fraction(mean(gap_fractions))}",
        "min internal_gap_fraction: "
        f"{format_fraction(min(internal_gap_fractions))}",
        "max internal_gap_fraction: "
        f"{format_fraction(max(internal_gap_fractions))}",
        "mean internal_gap_fraction: "
        f"{format_fraction(mean(internal_gap_fractions))}",
        "records where non_gap_length != raw_length: "
        f"{len(mismatches)}",
        "",
        "Top 20 sequences by gap_fraction:",
        "protein_id species raw_length aligned_length "
        "gap_fraction internal_gap_fraction",
    ]

    lines.extend(format_top_record(record) for record in top_by_gap)
    lines.extend(
        [
            "",
            "Top 20 sequences by internal_gap_fraction:",
            "protein_id species raw_length aligned_length "
            "gap_fraction internal_gap_fraction",
        ]
    )
    lines.extend(format_top_record(record) for record in top_by_internal_gap)

    if mismatches:
        lines.extend(
            [
                "",
                "Records where non_gap_length != raw_length:",
                "protein_id species raw_length non_gap_length",
            ]
        )
        lines.extend(
            f"{record.protein_id} {record.species} "
            f"{record.raw_length} {record.non_gap_length}"
            for record in mismatches
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run() -> None:
    """
    Generate alignment QC TSV and summary outputs.
    """

    aligned_records = parse_aligned_fasta(ALIGNMENT_FASTA)
    metadata_by_protein = load_metadata(METADATA_TSV)
    alignment_length = validate_aligned_lengths(aligned_records)
    qc_records = build_qc_records(aligned_records, metadata_by_protein)

    write_qc_tsv(qc_records, QC_TSV)
    write_summary(qc_records, alignment_length, SUMMARY_TXT)

    mismatches = get_length_mismatches(qc_records)

    print(f"sequences: {len(qc_records)}")
    print(f"alignment_length: {alignment_length}")
    print(f"non_gap_length_mismatches: {len(mismatches)}")
    print(f"saved: {QC_TSV}")
    print(f"saved: {SUMMARY_TXT}")


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

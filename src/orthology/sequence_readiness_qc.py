"""
Classify strict PHF10 MSA sequences by downstream analysis readiness.

Purpose:
    Summarize technical readiness for tree analysis, N-terminal analysis,
    and PTM/phospho-cluster analysis using the strict MAFFT alignment.

Input files:
    data/interim/orthology/msa/msa_input_strict_with_reference.mafft.fasta
    data/interim/orthology/msa_input_strict_with_reference.tsv

Output files:
    data/interim/orthology/msa/
    msa_input_strict_with_reference.sequence_readiness_qc.tsv
    data/interim/orthology/msa/
    msa_input_strict_with_reference.sequence_readiness_qc_summary.txt

Pipeline position:
    build_msa_input
        ->
    run_msa
        ->
    alignment_qc
        ->
    sequence_readiness_qc

Notes:
    This module performs readiness QC only. It does not filter
    sequences, modify the alignment, define orthology, use Ser-rich
    sequence, N-terminal extension, phosphosite motifs, domains, or
    perform biological interpretation. TRUNCATED means not ready for a
    downstream analysis window, not junk and not non-orthology.
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any


ALIGNMENT_FASTA = Path(
    "data/interim/orthology/msa/"
    "msa_input_strict_with_reference.mafft.fasta"
)
METADATA_TSV = Path(
    "data/interim/orthology/msa_input_strict_with_reference.tsv"
)
QC_TSV = Path(
    "data/interim/orthology/msa/"
    "msa_input_strict_with_reference.sequence_readiness_qc.tsv"
)
SUMMARY_TXT = Path(
    "data/interim/orthology/msa/"
    "msa_input_strict_with_reference.sequence_readiness_qc_summary.txt"
)

NTERM_WINDOW_START = 1
NTERM_WINDOW_END = 80
PTM_WINDOW_START = 1
PTM_WINDOW_END = 80
MIN_NTERM_COVERAGE_FOR_OK = 0.80
MIN_NTERM_COVERAGE_FOR_REVIEW = 0.30
MAX_X_FRACTION_FOR_TREE_OK = 0.05

OUTPUT_COLUMNS = [
    "protein_id",
    "species",
    "gene_id",
    "group",
    "raw_length",
    "aligned_length",
    "non_gap_length",
    "species_record_count",
    "is_multi_record_species",
    "x_count_total",
    "x_fraction_total",
    "x_count_nterm_window",
    "x_count_ptm_window",
    "first_aligned_human_position",
    "last_aligned_human_position",
    "nterm_window_positions",
    "nterm_window_covered_positions",
    "nterm_window_coverage",
    "ptm_window_positions",
    "ptm_window_covered_positions",
    "ptm_window_coverage",
    "tree_readiness",
    "nterm_readiness",
    "ptm_readiness",
    "sequence_qc_flags",
    "sequence_qc_note",
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
class WindowMetrics:
    """
    Coverage and X-count metrics for a human-reference window.
    """

    positions: int
    covered_positions: int
    coverage: float
    x_count: int


@dataclass(frozen=True)
class ReadinessRecord:
    """
    Readiness QC metrics for one aligned sequence.
    """

    protein_id: str
    species: str
    gene_id: str
    group: str
    raw_length: int
    aligned_length: int
    non_gap_length: int
    species_record_count: int
    is_multi_record_species: bool
    x_count_total: int
    x_fraction_total: float
    x_count_nterm_window: int
    x_count_ptm_window: int
    first_aligned_human_position: int | None
    last_aligned_human_position: int | None
    nterm_window_positions: int
    nterm_window_covered_positions: int
    nterm_window_coverage: float
    ptm_window_positions: int
    ptm_window_covered_positions: int
    ptm_window_coverage: float
    tree_readiness: str
    nterm_readiness: str
    ptm_readiness: str
    sequence_qc_flags: list[str]
    sequence_qc_note: str


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


def validate_alignment_lengths(records: list[AlignedRecord]) -> int:
    """
    Ensure all aligned sequences have the same aligned length.
    """

    lengths = {len(record.sequence) for record in records}

    if len(lengths) != 1:
        raise ValueError(
            "Aligned sequence lengths differ between records: "
            + ", ".join(str(length) for length in sorted(lengths))
        )

    return lengths.pop()


def find_human_reference(
    metadata_by_protein: dict[str, dict[str, str]],
    aligned_by_protein: dict[str, AlignedRecord],
) -> AlignedRecord:
    """
    Find exactly one Homo sapiens aligned record for reference positions.
    """

    human_rows = [
        row
        for row in metadata_by_protein.values()
        if row.get("species") == "homo_sapiens"
    ]

    if len(human_rows) != 1:
        raise ValueError(
            "Expected exactly one Homo sapiens sequence with "
            'species == "homo_sapiens" in metadata, found '
            f"{len(human_rows)}"
        )

    protein_id = human_rows[0]["protein_id"]
    reference = aligned_by_protein.get(protein_id)

    if reference is None:
        raise ValueError(
            f"Homo sapiens protein_id is missing from alignment: {protein_id}"
        )

    return reference


def build_human_position_map(human_sequence: str) -> list[int | None]:
    """
    Map alignment columns to human ungapped residue positions.
    """

    position_map: list[int | None] = []
    human_position = 0

    for character in human_sequence:
        if character == "-":
            position_map.append(None)
        else:
            human_position += 1
            position_map.append(human_position)

    return position_map


def get_species_counts(
    metadata_by_protein: dict[str, dict[str, str]]
) -> Counter[str]:
    """
    Count metadata records per species.
    """

    return Counter(
        row.get("species", "")
        for row in metadata_by_protein.values()
    )


def get_aligned_human_span(
    sequence: str,
    human_position_map: list[int | None],
) -> tuple[int | None, int | None]:
    """
    Return first and last human positions where target has a residue.
    """

    covered_positions = [
        human_position
        for character, human_position in zip(sequence, human_position_map)
        if character != "-" and human_position is not None
    ]

    if not covered_positions:
        return None, None

    return min(covered_positions), max(covered_positions)


def calculate_window_metrics(
    sequence: str,
    human_position_map: list[int | None],
    start: int,
    end: int,
) -> WindowMetrics:
    """
    Calculate coverage and X count for one human-reference window.
    """

    window_columns = [
        index
        for index, human_position in enumerate(human_position_map)
        if human_position is not None and start <= human_position <= end
    ]
    positions = len(window_columns)
    covered_positions = sum(
        1
        for index in window_columns
        if sequence[index] != "-"
    )
    x_count = sum(
        1
        for index in window_columns
        if sequence[index] in ("X", "x")
    )
    coverage = covered_positions / positions if positions else 0.0

    return WindowMetrics(
        positions=positions,
        covered_positions=covered_positions,
        coverage=coverage,
        x_count=x_count,
    )


def classify_tree_readiness(x_fraction_total: float) -> str:
    """
    Classify tree readiness from total X fraction only.
    """

    if x_fraction_total > MAX_X_FRACTION_FOR_TREE_OK:
        return "TREE_REVIEW"

    return "TREE_OK"


def classify_nterm_readiness(metrics: WindowMetrics) -> str:
    """
    Classify N-terminal analysis readiness from coverage and X count.
    """

    if metrics.coverage < MIN_NTERM_COVERAGE_FOR_REVIEW:
        return "NTERM_TRUNCATED"

    if (
        metrics.coverage >= MIN_NTERM_COVERAGE_FOR_OK
        and metrics.x_count == 0
    ):
        return "NTERM_OK"

    return "NTERM_REVIEW"


def classify_ptm_readiness(metrics: WindowMetrics) -> str:
    """
    Classify PTM-window readiness from coverage and X count.
    """

    if metrics.coverage < MIN_NTERM_COVERAGE_FOR_REVIEW:
        return "PTM_NOT_READY"

    if (
        metrics.coverage >= MIN_NTERM_COVERAGE_FOR_OK
        and metrics.x_count == 0
    ):
        return "PTM_OK"

    return "PTM_REVIEW"


def build_flags(
    species_record_count: int,
    x_count_total: int,
    nterm_metrics: WindowMetrics,
    ptm_metrics: WindowMetrics,
) -> list[str]:
    """
    Build semicolon-separated QC flags for one sequence.
    """

    flags = []

    if species_record_count > 1:
        flags.append("MULTIPLE_RECORDS_PER_SPECIES")

    if x_count_total > 0:
        flags.append("HAS_X")

    if nterm_metrics.x_count > 0:
        flags.append("X_IN_NTERM_WINDOW")

    if ptm_metrics.x_count > 0:
        flags.append("X_IN_PTM_WINDOW")

    if nterm_metrics.coverage < MIN_NTERM_COVERAGE_FOR_REVIEW:
        flags.append("LOW_NTERM_COVERAGE")

    if ptm_metrics.coverage < MIN_NTERM_COVERAGE_FOR_REVIEW:
        flags.append("LOW_PTM_COVERAGE")

    return flags


def build_note(record: ReadinessRecord) -> str:
    """
    Create a concise non-interpretive QC note.
    """

    notes = []

    if record.is_multi_record_species:
        notes.append("multiple records for species")

    if record.x_count_total > 0:
        notes.append("contains X residues")

    if record.nterm_readiness == "NTERM_TRUNCATED":
        notes.append("low N-terminal window coverage")

    if record.ptm_readiness == "PTM_NOT_READY":
        notes.append("low PTM window coverage")

    return "; ".join(notes)


def build_readiness_record(
    aligned_record: AlignedRecord,
    metadata: dict[str, str],
    species_record_count: int,
    human_position_map: list[int | None],
) -> ReadinessRecord:
    """
    Calculate sequence readiness QC metrics for one aligned record.
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
            f"Aligned sequence contains no non-gap residues: "
            f"{aligned_record.protein_id}"
        )

    x_count_total = sum(
        1
        for char in aligned_record.sequence
        if char in ("X", "x")
    )
    x_fraction_total = x_count_total / non_gap_length
    nterm_metrics = calculate_window_metrics(
        aligned_record.sequence,
        human_position_map,
        NTERM_WINDOW_START,
        NTERM_WINDOW_END,
    )
    ptm_metrics = calculate_window_metrics(
        aligned_record.sequence,
        human_position_map,
        PTM_WINDOW_START,
        PTM_WINDOW_END,
    )
    first_position, last_position = get_aligned_human_span(
        aligned_record.sequence,
        human_position_map,
    )
    flags = build_flags(
        species_record_count=species_record_count,
        x_count_total=x_count_total,
        nterm_metrics=nterm_metrics,
        ptm_metrics=ptm_metrics,
    )
    record = ReadinessRecord(
        protein_id=aligned_record.protein_id,
        species=metadata.get("species", ""),
        gene_id=metadata.get("gene_id", ""),
        group=metadata.get("group", ""),
        raw_length=raw_length,
        aligned_length=aligned_length,
        non_gap_length=non_gap_length,
        species_record_count=species_record_count,
        is_multi_record_species=species_record_count > 1,
        x_count_total=x_count_total,
        x_fraction_total=x_fraction_total,
        x_count_nterm_window=nterm_metrics.x_count,
        x_count_ptm_window=ptm_metrics.x_count,
        first_aligned_human_position=first_position,
        last_aligned_human_position=last_position,
        nterm_window_positions=nterm_metrics.positions,
        nterm_window_covered_positions=nterm_metrics.covered_positions,
        nterm_window_coverage=nterm_metrics.coverage,
        ptm_window_positions=ptm_metrics.positions,
        ptm_window_covered_positions=ptm_metrics.covered_positions,
        ptm_window_coverage=ptm_metrics.coverage,
        tree_readiness=classify_tree_readiness(x_fraction_total),
        nterm_readiness=classify_nterm_readiness(nterm_metrics),
        ptm_readiness=classify_ptm_readiness(ptm_metrics),
        sequence_qc_flags=flags,
        sequence_qc_note="",
    )

    return ReadinessRecord(
        **{
            **record.__dict__,
            "sequence_qc_note": build_note(record),
        }
    )


def build_readiness_records(
    aligned_records: list[AlignedRecord],
    metadata_by_protein: dict[str, dict[str, str]],
    human_position_map: list[int | None],
) -> list[ReadinessRecord]:
    """
    Calculate readiness QC for all aligned records.
    """

    species_counts = get_species_counts(metadata_by_protein)
    records = []

    for aligned_record in aligned_records:
        metadata = metadata_by_protein.get(aligned_record.protein_id)

        if metadata is None:
            raise ValueError(
                "Aligned protein_id is missing from metadata: "
                f"{aligned_record.protein_id}"
            )

        species = metadata.get("species", "")
        records.append(
            build_readiness_record(
                aligned_record=aligned_record,
                metadata=metadata,
                species_record_count=species_counts[species],
                human_position_map=human_position_map,
            )
        )

    return records


def format_fraction(value: float) -> str:
    """
    Format fractions consistently.
    """

    return f"{value:.6f}"


def format_optional_int(value: int | None) -> str:
    """
    Format optional integer positions for TSV output.
    """

    if value is None:
        return ""

    return str(value)


def readiness_record_to_row(record: ReadinessRecord) -> dict[str, Any]:
    """
    Convert one readiness record to a TSV row.
    """

    return {
        "protein_id": record.protein_id,
        "species": record.species,
        "gene_id": record.gene_id,
        "group": record.group,
        "raw_length": record.raw_length,
        "aligned_length": record.aligned_length,
        "non_gap_length": record.non_gap_length,
        "species_record_count": record.species_record_count,
        "is_multi_record_species": record.is_multi_record_species,
        "x_count_total": record.x_count_total,
        "x_fraction_total": format_fraction(record.x_fraction_total),
        "x_count_nterm_window": record.x_count_nterm_window,
        "x_count_ptm_window": record.x_count_ptm_window,
        "first_aligned_human_position": format_optional_int(
            record.first_aligned_human_position
        ),
        "last_aligned_human_position": format_optional_int(
            record.last_aligned_human_position
        ),
        "nterm_window_positions": record.nterm_window_positions,
        "nterm_window_covered_positions": (
            record.nterm_window_covered_positions
        ),
        "nterm_window_coverage": format_fraction(
            record.nterm_window_coverage
        ),
        "ptm_window_positions": record.ptm_window_positions,
        "ptm_window_covered_positions": record.ptm_window_covered_positions,
        "ptm_window_coverage": format_fraction(record.ptm_window_coverage),
        "tree_readiness": record.tree_readiness,
        "nterm_readiness": record.nterm_readiness,
        "ptm_readiness": record.ptm_readiness,
        "sequence_qc_flags": ";".join(record.sequence_qc_flags),
        "sequence_qc_note": record.sequence_qc_note,
    }


def write_qc_tsv(records: list[ReadinessRecord], output_path: Path) -> None:
    """
    Write sequence readiness QC records to TSV.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
        )
        writer.writeheader()

        for record in records:
            writer.writerow(readiness_record_to_row(record))


def group_species_records(
    records: list[ReadinessRecord],
) -> dict[str, list[ReadinessRecord]]:
    """
    Group readiness records by species.
    """

    grouped: dict[str, list[ReadinessRecord]] = defaultdict(list)

    for record in records:
        grouped[record.species].append(record)

    return grouped


def count_by_field(records: list[ReadinessRecord], field: str) -> Counter[str]:
    """
    Count readiness states for a named dataclass field.
    """

    return Counter(str(getattr(record, field)) for record in records)


def format_counter(counter: Counter[str]) -> list[str]:
    """
    Format readiness count lines.
    """

    return [
        f"{key}: {counter[key]}"
        for key in sorted(counter)
    ]


def write_summary(records: list[ReadinessRecord], output_path: Path) -> None:
    """
    Write text summary for sequence readiness QC.
    """

    species_groups = group_species_records(records)
    reference_records = [
        record
        for record in records
        if record.group == "REFERENCE"
    ]
    non_reference_records = [
        record
        for record in records
        if record.group != "REFERENCE"
    ]
    multi_species = {
        species: species_records
        for species, species_records in species_groups.items()
        if len(species_records) > 1
    }
    records_with_x = [record for record in records if record.x_count_total > 0]
    records_with_x_nterm = [
        record for record in records if record.x_count_nterm_window > 0
    ]
    records_with_x_ptm = [
        record for record in records if record.x_count_ptm_window > 0
    ]
    nterm_truncated = [
        record
        for record in records
        if record.nterm_readiness == "NTERM_TRUNCATED"
    ]
    ptm_not_ready = [
        record
        for record in records
        if record.ptm_readiness == "PTM_NOT_READY"
    ]

    lines = [
        "PHF10 sequence readiness QC",
        "QC only: no orthology redefinition and no biological interpretation.",
        "Homo sapiens is used only as the coordinate reference.",
        "",
        f"total records: {len(records)}",
        f"reference records: {len(reference_records)}",
        f"non-reference records: {len(non_reference_records)}",
        f"number of sequences: {len(records)}",
        f"number of species: {len(species_groups)}",
        f"species with multiple records: {len(multi_species)}",
        "",
        "tree_readiness counts",
        *format_counter(count_by_field(records, "tree_readiness")),
        "",
        "nterm_readiness counts",
        *format_counter(count_by_field(records, "nterm_readiness")),
        "",
        "ptm_readiness counts",
        *format_counter(count_by_field(records, "ptm_readiness")),
        "",
        f"records with X: {len(records_with_x)}",
        "records with X in N-terminal window: "
        f"{len(records_with_x_nterm)}",
        f"records with X in PTM window: {len(records_with_x_ptm)}",
        f"records with NTERM_TRUNCATED: {len(nterm_truncated)}",
        f"records with PTM_NOT_READY: {len(ptm_not_ready)}",
        "",
        "Species with multiple records:",
        "species count protein_ids",
    ]

    if multi_species:
        for species in sorted(multi_species):
            protein_ids = sorted(
                record.protein_id
                for record in multi_species[species]
            )
            lines.append(
                f"{species} {len(protein_ids)} {';'.join(protein_ids)}"
            )
    else:
        lines.append("None")

    lines.extend(
        [
            "",
            "NTERM_TRUNCATED records:",
            "protein_id species raw_length nterm_window_coverage "
            "first_aligned_human_position",
        ]
    )
    if nterm_truncated:
        for record in nterm_truncated:
            lines.append(
                f"{record.protein_id} {record.species} {record.raw_length} "
                f"{format_fraction(record.nterm_window_coverage)} "
                f"{format_optional_int(record.first_aligned_human_position)}"
            )
    else:
        lines.append("None")

    lines.extend(
        [
            "",
            "PTM_NOT_READY records:",
            "protein_id species raw_length ptm_window_coverage "
            "first_aligned_human_position",
        ]
    )
    if ptm_not_ready:
        for record in ptm_not_ready:
            lines.append(
                f"{record.protein_id} {record.species} {record.raw_length} "
                f"{format_fraction(record.ptm_window_coverage)} "
                f"{format_optional_int(record.first_aligned_human_position)}"
            )
    else:
        lines.append("None")

    lines.extend(
        [
            "",
            "Records with X in N-terminal or PTM window:",
            "protein_id species x_count_nterm_window x_count_ptm_window",
        ]
    )
    x_window_records_by_id = {
        record.protein_id: record
        for record in records_with_x_nterm + records_with_x_ptm
    }
    x_window_records = sorted(
        x_window_records_by_id.values(),
        key=lambda record: record.protein_id,
    )
    if x_window_records:
        for record in x_window_records:
            lines.append(
                f"{record.protein_id} {record.species} "
                f"{record.x_count_nterm_window} {record.x_count_ptm_window}"
            )
    else:
        lines.append("None")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run() -> None:
    """
    Generate sequence readiness QC outputs.
    """

    aligned_records = parse_aligned_fasta(ALIGNMENT_FASTA)
    metadata_by_protein = load_metadata(METADATA_TSV)
    aligned_by_protein = {
        record.protein_id: record
        for record in aligned_records
    }
    validate_alignment_lengths(aligned_records)
    human_reference = find_human_reference(
        metadata_by_protein,
        aligned_by_protein,
    )
    human_position_map = build_human_position_map(human_reference.sequence)
    records = build_readiness_records(
        aligned_records,
        metadata_by_protein,
        human_position_map,
    )

    write_qc_tsv(records, QC_TSV)
    write_summary(records, SUMMARY_TXT)

    nterm_truncated = sum(
        1
        for record in records
        if record.nterm_readiness == "NTERM_TRUNCATED"
    )
    ptm_not_ready = sum(
        1
        for record in records
        if record.ptm_readiness == "PTM_NOT_READY"
    )
    species_count = len({record.species for record in records})

    print(f"sequences: {len(records)}")
    print(f"species: {species_count}")
    print(f"nterm_truncated: {nterm_truncated}")
    print(f"ptm_not_ready: {ptm_not_ready}")
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

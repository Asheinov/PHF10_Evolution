"""
Annotate PHF10 N-terminal architecture relative to the human reference.

Purpose:
    Use the reference-aware MAFFT alignment to calculate human-coordinate
    N-terminal window coverage, S50-equivalent status, phospho-cluster
    coverage, and PTM-context readiness annotations for downstream review.

Input files:
    data/interim/orthology/msa/msa_input_strict_with_reference.mafft.fasta
    data/interim/orthology/msa/
    msa_input_strict_with_reference.sequence_readiness_qc.tsv

Output files:
    data/interim/orthology/msa/
    msa_input_strict_with_reference.nterm_form_qc.tsv
    data/interim/orthology/msa/
    msa_input_strict_with_reference.nterm_form_qc_summary.txt

Pipeline position:
    run_msa_with_reference
        ->
    sequence_readiness_qc
        ->
    nterm_form_qc

Notes:
    This module is QC-only. It does not redefine orthology, filter records,
    infer evolutionary gain or loss, classify confirmed biological isoforms,
    or use domains, IDR, Ser-rich sequence, N-terminal extension, PTM motifs,
    or trait-level information as orthology criteria. The Q8WUB8 human PHF10
    reference is used only for coordinate mapping and is not treated as an
    ortholog.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


ALIGNMENT_FASTA = Path(
    "data/interim/orthology/msa/"
    "msa_input_strict_with_reference.mafft.fasta"
)
READINESS_QC_TSV = Path(
    "data/interim/orthology/msa/"
    "msa_input_strict_with_reference.sequence_readiness_qc.tsv"
)
OUTPUT_QC_TSV = Path(
    "data/interim/orthology/msa/"
    "msa_input_strict_with_reference.nterm_form_qc.tsv"
)
OUTPUT_SUMMARY_TXT = Path(
    "data/interim/orthology/msa/"
    "msa_input_strict_with_reference.nterm_form_qc_summary.txt"
)

REFERENCE_PROTEIN_ID = "Q8WUB8"

UPSTREAM_EXTENSION_START = 1
UPSTREAM_EXTENSION_END = 37
PHOSPHO_CLUSTER_START = 38
PHOSPHO_CLUSTER_END = 65
S50_POSITION = 50
DOWNSTREAM_ANCHOR_START = 80
DOWNSTREAM_ANCHOR_END = 120

REQUIRED_QC_FIELDS = {
    "protein_id",
    "species",
    "gene_id",
    "raw_length",
    "tree_readiness",
    "nterm_readiness",
    "ptm_readiness",
    "first_aligned_human_position",
    "nterm_window_coverage",
    "ptm_window_coverage",
}
OUTPUT_FIELDNAMES = [
    "protein_id",
    "species",
    "gene_id",
    "raw_length",
    "manual_status",
    "tree_readiness",
    "nterm_readiness",
    "ptm_readiness",
    "first_aligned_human_position",
    "nterm_window_coverage",
    "ptm_window_coverage",
    "upstream_extension_coverage",
    "phospho_cluster_coverage",
    "downstream_anchor_coverage",
    "upstream_extension_x_count",
    "phospho_cluster_x_count",
    "s50_equivalent_residue",
    "s50_equivalent_status",
    "nterm_architecture_class",
    "sequence_completeness_class",
    "ptm_context_class",
    "nterm_form_flags",
]


@dataclass(frozen=True)
class FastaRecord:
    """
    One aligned FASTA record.
    """

    protein_id: str
    header: str
    sequence: str


@dataclass(frozen=True)
class ReadinessRecord:
    """
    Fields from one sequence readiness QC row needed for N-terminal form QC.
    """

    protein_id: str
    species: str
    gene_id: str
    raw_length: str
    manual_status: str
    tree_readiness: str
    nterm_readiness: str
    ptm_readiness: str
    first_aligned_human_position: str
    nterm_window_coverage: str
    ptm_window_coverage: str


@dataclass(frozen=True)
class WindowMetrics:
    """
    Coverage and X-count for one human-coordinate window.
    """

    coverage: float
    x_count: int


@dataclass(frozen=True)
class NtermFormRecord:
    """
    N-terminal form QC annotation for one readiness input row.
    """

    protein_id: str
    species: str
    gene_id: str
    raw_length: str
    manual_status: str
    tree_readiness: str
    nterm_readiness: str
    ptm_readiness: str
    first_aligned_human_position: str
    nterm_window_coverage: str
    ptm_window_coverage: str
    upstream_extension_coverage: float
    phospho_cluster_coverage: float
    downstream_anchor_coverage: float
    upstream_extension_x_count: int
    phospho_cluster_x_count: int
    s50_equivalent_residue: str
    s50_equivalent_status: str
    nterm_architecture_class: str
    sequence_completeness_class: str
    ptm_context_class: str
    nterm_form_flags: list[str]

    def to_dict(self) -> dict[str, str]:
        """
        Convert the QC record to a DictWriter-compatible dictionary.
        """

        return {
            "protein_id": self.protein_id,
            "species": self.species,
            "gene_id": self.gene_id,
            "raw_length": self.raw_length,
            "manual_status": self.manual_status,
            "tree_readiness": self.tree_readiness,
            "nterm_readiness": self.nterm_readiness,
            "ptm_readiness": self.ptm_readiness,
            "first_aligned_human_position": self.first_aligned_human_position,
            "nterm_window_coverage": self.nterm_window_coverage,
            "ptm_window_coverage": self.ptm_window_coverage,
            "upstream_extension_coverage": (
                f"{self.upstream_extension_coverage:.6f}"
            ),
            "phospho_cluster_coverage": f"{self.phospho_cluster_coverage:.6f}",
            "downstream_anchor_coverage": (
                f"{self.downstream_anchor_coverage:.6f}"
            ),
            "upstream_extension_x_count": str(self.upstream_extension_x_count),
            "phospho_cluster_x_count": str(self.phospho_cluster_x_count),
            "s50_equivalent_residue": self.s50_equivalent_residue,
            "s50_equivalent_status": self.s50_equivalent_status,
            "nterm_architecture_class": self.nterm_architecture_class,
            "sequence_completeness_class": self.sequence_completeness_class,
            "ptm_context_class": self.ptm_context_class,
            "nterm_form_flags": ";".join(self.nterm_form_flags),
        }


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


def parse_fasta(path: Path) -> list[FastaRecord]:
    """
    Parse aligned FASTA records using only the standard library.
    """

    require_input_file(path)
    records = []
    header = ""
    chunks: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    records.append(
                        FastaRecord(
                            protein_id=header.split("|", 1)[0],
                            header=header,
                            sequence="".join(chunks),
                        )
                    )
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)

    if header:
        records.append(
            FastaRecord(
                protein_id=header.split("|", 1)[0],
                header=header,
                sequence="".join(chunks),
            )
        )

    if not records:
        raise ValueError(f"Aligned FASTA has no records: {path}")

    return records


def load_readiness_records(path: Path) -> list[ReadinessRecord]:
    """
    Load sequence readiness QC records from a TSV file.
    """

    require_input_file(path)
    records = []

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_fields(path, reader.fieldnames, REQUIRED_QC_FIELDS)

        for row in reader:
            protein_id = row["protein_id"].strip()
            if not protein_id:
                raise ValueError(f"Readiness QC row has empty protein_id: {path}")

            manual_status = row.get("manual_status", "").strip()
            if not manual_status and protein_id == REFERENCE_PROTEIN_ID:
                manual_status = "REFERENCE"

            records.append(
                ReadinessRecord(
                    protein_id=protein_id,
                    species=row["species"].strip(),
                    gene_id=row["gene_id"].strip(),
                    raw_length=row["raw_length"].strip(),
                    manual_status=manual_status,
                    tree_readiness=row["tree_readiness"].strip(),
                    nterm_readiness=row["nterm_readiness"].strip(),
                    ptm_readiness=row["ptm_readiness"].strip(),
                    first_aligned_human_position=row[
                        "first_aligned_human_position"
                    ].strip(),
                    nterm_window_coverage=row["nterm_window_coverage"].strip(),
                    ptm_window_coverage=row["ptm_window_coverage"].strip(),
                )
            )

    if not records:
        raise ValueError(f"Input TSV has no readiness QC records: {path}")

    return records


def index_alignment(records: list[FastaRecord]) -> dict[str, FastaRecord]:
    """
    Index aligned FASTA records by protein_id and fail on duplicates.
    """

    indexed = {}
    for record in records:
        if record.protein_id in indexed:
            raise ValueError(
                f"Duplicate protein_id in aligned FASTA: {record.protein_id}"
            )
        indexed[record.protein_id] = record
    return indexed


def get_reference_record(records: list[FastaRecord]) -> FastaRecord:
    """
    Return the unique Q8WUB8 reference alignment record.
    """

    reference_records = [
        record
        for record in records
        if record.protein_id == REFERENCE_PROTEIN_ID
    ]
    if not reference_records:
        raise ValueError(
            f"{REFERENCE_PROTEIN_ID} reference is absent from the alignment"
        )
    if len(reference_records) > 1:
        raise ValueError(
            f"Multiple {REFERENCE_PROTEIN_ID} reference records found"
        )
    return reference_records[0]


def map_alignment_columns_to_human_positions(
    reference_sequence: str,
) -> dict[int, int]:
    """
    Map 0-based alignment columns to 1-based ungapped human positions.
    """

    column_to_human_position = {}
    human_position = 0

    for column_index, residue in enumerate(reference_sequence):
        if residue != "-":
            human_position += 1
            column_to_human_position[column_index] = human_position

    return column_to_human_position


def map_human_positions_to_alignment_columns(
    column_to_human_position: dict[int, int],
) -> dict[int, int]:
    """
    Map 1-based ungapped human positions to 0-based alignment columns.
    """

    return {
        human_position: column_index
        for column_index, human_position in column_to_human_position.items()
    }


def window_columns(
    column_to_human_position: dict[int, int],
    start: int,
    end: int,
) -> list[int]:
    """
    Return alignment columns corresponding to a human-coordinate window.
    """

    return [
        column_index
        for column_index, human_position in column_to_human_position.items()
        if start <= human_position <= end
    ]


def calculate_window_metrics(
    sequence: str,
    columns: list[int],
    window_length: int,
) -> WindowMetrics:
    """
    Calculate coverage and X-count for one mapped human-coordinate window.
    """

    covered_positions = 0
    x_count = 0

    for column_index in columns:
        residue = sequence[column_index].upper()
        if residue != "-":
            covered_positions += 1
        if residue == "X":
            x_count += 1

    return WindowMetrics(
        coverage=covered_positions / window_length,
        x_count=x_count,
    )


def get_s50_residue(
    sequence: str,
    human_position_to_column: dict[int, int],
) -> str:
    """
    Extract the residue aligned to human S50, if that column is mapped.
    """

    column_index = human_position_to_column.get(S50_POSITION)
    if column_index is None:
        return ""
    return sequence[column_index].upper()


def classify_s50_status(residue: str) -> str:
    """
    Classify the residue aligned to human S50.
    """

    if not residue:
        return "S50_NOT_ASSESSABLE"
    if residue == "S":
        return "S50_PRESENT"
    if residue == "-":
        return "S50_GAP"
    if residue == "X":
        return "S50_X"
    return "S50_SUBSTITUTED"


def classify_nterm_architecture(
    protein_id: str,
    upstream_extension_coverage: float,
    phospho_cluster_coverage: float,
    downstream_anchor_coverage: float,
) -> str:
    """
    Classify N-terminal architecture from human-coordinate coverage metrics.
    """

    if protein_id == REFERENCE_PROTEIN_ID:
        return "REFERENCE"
    if upstream_extension_coverage >= 0.70 and phospho_cluster_coverage >= 0.70:
        return "FULL_NTERM_EXTENSION"
    if upstream_extension_coverage < 0.30 and phospho_cluster_coverage >= 0.70:
        return "MGS_START_FORM"
    if phospho_cluster_coverage < 0.30 and downstream_anchor_coverage >= 0.70:
        return "DOWNSTREAM_START_FORM"
    if upstream_extension_coverage >= 0.30 and phospho_cluster_coverage >= 0.30:
        return "NTERM_DIVERGENT"
    return "NTERM_NOT_ASSESSABLE"


def classify_sequence_completeness(
    protein_id: str,
    nterm_architecture_class: str,
    downstream_anchor_coverage: float,
) -> str:
    """
    Classify sequence completeness context from architecture and anchor metrics.
    """

    if protein_id == REFERENCE_PROTEIN_ID:
        return "REFERENCE"
    if (
        nterm_architecture_class == "FULL_NTERM_EXTENSION"
        and downstream_anchor_coverage >= 0.70
    ):
        return "COMPLETE_LIKE"
    if (
        nterm_architecture_class == "MGS_START_FORM"
        and downstream_anchor_coverage >= 0.70
    ):
        return "POSSIBLE_SHORT_ISOFORM"
    if nterm_architecture_class == "DOWNSTREAM_START_FORM":
        return "LIKELY_TRUNCATED_MODEL"
    if downstream_anchor_coverage < 0.70:
        return "FRAGMENT_REVIEW"
    return "COMPLETENESS_REVIEW"


def classify_ptm_context(
    protein_id: str,
    nterm_architecture_class: str,
    phospho_cluster_coverage: float,
    phospho_cluster_x_count: int,
    s50_equivalent_status: str,
) -> str:
    """
    Classify PTM-context readiness without making biological conclusions.
    """

    if protein_id == REFERENCE_PROTEIN_ID:
        return "REFERENCE"
    if (
        nterm_architecture_class == "FULL_NTERM_EXTENSION"
        and phospho_cluster_coverage >= 0.70
        and s50_equivalent_status == "S50_PRESENT"
        and phospho_cluster_x_count == 0
    ):
        return "PTM_FUNCTIONAL_READY"
    if (
        nterm_architecture_class == "MGS_START_FORM"
        and phospho_cluster_coverage >= 0.70
        and s50_equivalent_status == "S50_PRESENT"
    ):
        return "PTM_SEQUENCE_PRESENT_CONTEXT_MISSING"
    if (
        phospho_cluster_coverage >= 0.70
        and s50_equivalent_status in {"S50_SUBSTITUTED", "S50_X"}
    ):
        return "PTM_CONTEXT_REVIEW"
    if phospho_cluster_coverage < 0.30 or s50_equivalent_status == "S50_GAP":
        return "PTM_NOT_READY"
    return "PTM_NOT_ASSESSABLE"


def build_nterm_form_flags(
    upstream_extension_coverage: float,
    phospho_cluster_coverage: float,
    downstream_anchor_coverage: float,
    phospho_cluster_x_count: int,
    s50_equivalent_status: str,
    nterm_architecture_class: str,
) -> list[str]:
    """
    Build semicolon-separated technical QC flags for one record.
    """

    flags = []
    if upstream_extension_coverage < 0.30:
        flags.append("LOW_UPSTREAM_EXTENSION_COVERAGE")
    if phospho_cluster_coverage < 0.30:
        flags.append("LOW_PHOSPHO_CLUSTER_COVERAGE")
    if downstream_anchor_coverage < 0.70:
        flags.append("LOW_DOWNSTREAM_ANCHOR_COVERAGE")
    if phospho_cluster_x_count > 0:
        flags.append("PHOSPHO_CLUSTER_HAS_X")
    if s50_equivalent_status == "S50_SUBSTITUTED":
        flags.append("S50_SUBSTITUTED")
    if s50_equivalent_status == "S50_X":
        flags.append("S50_X")
    if s50_equivalent_status == "S50_GAP":
        flags.append("S50_GAP")
    if nterm_architecture_class == "MGS_START_FORM":
        flags.append("MGS_START_CONTEXT_REVIEW")
    return flags


def annotate_record(
    readiness_record: ReadinessRecord,
    aligned_record: FastaRecord,
    upstream_columns: list[int],
    phospho_cluster_columns: list[int],
    downstream_anchor_columns: list[int],
    human_position_to_column: dict[int, int],
) -> NtermFormRecord:
    """
    Annotate one readiness record with N-terminal form QC metrics.
    """

    upstream_metrics = calculate_window_metrics(
        aligned_record.sequence,
        upstream_columns,
        UPSTREAM_EXTENSION_END - UPSTREAM_EXTENSION_START + 1,
    )
    phospho_cluster_metrics = calculate_window_metrics(
        aligned_record.sequence,
        phospho_cluster_columns,
        PHOSPHO_CLUSTER_END - PHOSPHO_CLUSTER_START + 1,
    )
    downstream_anchor_metrics = calculate_window_metrics(
        aligned_record.sequence,
        downstream_anchor_columns,
        DOWNSTREAM_ANCHOR_END - DOWNSTREAM_ANCHOR_START + 1,
    )
    s50_residue = get_s50_residue(
        aligned_record.sequence,
        human_position_to_column,
    )
    s50_status = classify_s50_status(s50_residue)
    nterm_architecture_class = classify_nterm_architecture(
        protein_id=readiness_record.protein_id,
        upstream_extension_coverage=upstream_metrics.coverage,
        phospho_cluster_coverage=phospho_cluster_metrics.coverage,
        downstream_anchor_coverage=downstream_anchor_metrics.coverage,
    )
    sequence_completeness_class = classify_sequence_completeness(
        protein_id=readiness_record.protein_id,
        nterm_architecture_class=nterm_architecture_class,
        downstream_anchor_coverage=downstream_anchor_metrics.coverage,
    )
    ptm_context_class = classify_ptm_context(
        protein_id=readiness_record.protein_id,
        nterm_architecture_class=nterm_architecture_class,
        phospho_cluster_coverage=phospho_cluster_metrics.coverage,
        phospho_cluster_x_count=phospho_cluster_metrics.x_count,
        s50_equivalent_status=s50_status,
    )
    flags = build_nterm_form_flags(
        upstream_extension_coverage=upstream_metrics.coverage,
        phospho_cluster_coverage=phospho_cluster_metrics.coverage,
        downstream_anchor_coverage=downstream_anchor_metrics.coverage,
        phospho_cluster_x_count=phospho_cluster_metrics.x_count,
        s50_equivalent_status=s50_status,
        nterm_architecture_class=nterm_architecture_class,
    )

    return NtermFormRecord(
        protein_id=readiness_record.protein_id,
        species=readiness_record.species,
        gene_id=readiness_record.gene_id,
        raw_length=readiness_record.raw_length,
        manual_status=readiness_record.manual_status,
        tree_readiness=readiness_record.tree_readiness,
        nterm_readiness=readiness_record.nterm_readiness,
        ptm_readiness=readiness_record.ptm_readiness,
        first_aligned_human_position=(
            readiness_record.first_aligned_human_position
        ),
        nterm_window_coverage=readiness_record.nterm_window_coverage,
        ptm_window_coverage=readiness_record.ptm_window_coverage,
        upstream_extension_coverage=upstream_metrics.coverage,
        phospho_cluster_coverage=phospho_cluster_metrics.coverage,
        downstream_anchor_coverage=downstream_anchor_metrics.coverage,
        upstream_extension_x_count=upstream_metrics.x_count,
        phospho_cluster_x_count=phospho_cluster_metrics.x_count,
        s50_equivalent_residue=s50_residue,
        s50_equivalent_status=s50_status,
        nterm_architecture_class=nterm_architecture_class,
        sequence_completeness_class=sequence_completeness_class,
        ptm_context_class=ptm_context_class,
        nterm_form_flags=flags,
    )


def build_nterm_form_records(
    readiness_records: list[ReadinessRecord],
    alignment_by_protein_id: dict[str, FastaRecord],
    reference_record: FastaRecord,
) -> list[NtermFormRecord]:
    """
    Build one N-terminal form QC row per readiness QC row.
    """

    column_to_human_position = map_alignment_columns_to_human_positions(
        reference_record.sequence
    )
    human_position_to_column = map_human_positions_to_alignment_columns(
        column_to_human_position
    )
    upstream_columns = window_columns(
        column_to_human_position,
        UPSTREAM_EXTENSION_START,
        UPSTREAM_EXTENSION_END,
    )
    phospho_cluster_columns = window_columns(
        column_to_human_position,
        PHOSPHO_CLUSTER_START,
        PHOSPHO_CLUSTER_END,
    )
    downstream_anchor_columns = window_columns(
        column_to_human_position,
        DOWNSTREAM_ANCHOR_START,
        DOWNSTREAM_ANCHOR_END,
    )

    records = []
    for readiness_record in readiness_records:
        aligned_record = alignment_by_protein_id.get(readiness_record.protein_id)
        if aligned_record is None:
            raise ValueError(
                "QC protein_id is absent from the alignment: "
                f"{readiness_record.protein_id}"
            )

        records.append(
            annotate_record(
                readiness_record=readiness_record,
                aligned_record=aligned_record,
                upstream_columns=upstream_columns,
                phospho_cluster_columns=phospho_cluster_columns,
                downstream_anchor_columns=downstream_anchor_columns,
                human_position_to_column=human_position_to_column,
            )
        )

    if len(records) != len(readiness_records):
        raise ValueError(
            "Output row count differs from input QC row count: "
            f"{len(records)} != {len(readiness_records)}"
        )

    return records


def write_tsv(records: list[NtermFormRecord], output_path: Path) -> None:
    """
    Write N-terminal form QC rows to a tab-separated file.
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
        for record in records:
            writer.writerow(record.to_dict())


def write_summary(records: list[NtermFormRecord], output_path: Path) -> None:
    """
    Write a text summary of N-terminal form QC class counts.
    """

    reference_records = [
        record
        for record in records
        if record.protein_id == REFERENCE_PROTEIN_ID
    ]
    nterm_counts = Counter(
        record.nterm_architecture_class
        for record in records
    )
    completeness_counts = Counter(
        record.sequence_completeness_class
        for record in records
    )
    ptm_context_counts = Counter(
        record.ptm_context_class
        for record in records
    )
    s50_counts = Counter(
        record.s50_equivalent_status
        for record in records
    )

    lines = [
        "PHF10 N-terminal form QC",
        "QC only: no filtering, no orthology redefinition, no gain/loss inference.",
        "",
        f"total records: {len(records)}",
        f"reference records: {len(reference_records)}",
        f"non-reference records: {len(records) - len(reference_records)}",
        "",
        "nterm_architecture_class counts:",
    ]

    for label in sorted(nterm_counts):
        lines.append(f"{label}\t{nterm_counts[label]}")

    lines.extend(["", "sequence_completeness_class counts:"])
    for label in sorted(completeness_counts):
        lines.append(f"{label}\t{completeness_counts[label]}")

    lines.extend(["", "ptm_context_class counts:"])
    for label in sorted(ptm_context_counts):
        lines.append(f"{label}\t{ptm_context_counts[label]}")

    lines.extend(["", "s50_equivalent_status counts:"])
    for label in sorted(s50_counts):
        lines.append(f"{label}\t{s50_counts[label]}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """
    Run N-terminal form QC and write output files.
    """

    alignment_records = parse_fasta(ALIGNMENT_FASTA)
    readiness_records = load_readiness_records(READINESS_QC_TSV)
    alignment_by_protein_id = index_alignment(alignment_records)
    reference_record = get_reference_record(alignment_records)
    nterm_form_records = build_nterm_form_records(
        readiness_records=readiness_records,
        alignment_by_protein_id=alignment_by_protein_id,
        reference_record=reference_record,
    )

    write_tsv(nterm_form_records, OUTPUT_QC_TSV)
    write_summary(nterm_form_records, OUTPUT_SUMMARY_TXT)

    nterm_counts = Counter(
        record.nterm_architecture_class
        for record in nterm_form_records
    )
    print(f"total records: {len(nterm_form_records)}")
    print(f"reference records: {nterm_counts['REFERENCE']}")
    print(f"non-reference records: {len(nterm_form_records) - nterm_counts['REFERENCE']}")
    print(f"saved: {OUTPUT_QC_TSV}")
    print(f"saved: {OUTPUT_SUMMARY_TXT}")


if __name__ == "__main__":
    main()

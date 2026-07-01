"""
Summarize duplicate-species records in PHF10 sequence readiness QC.

Purpose:
    Identify species with multiple non-reference PHF10 candidate records
    in the reference-aware sequence readiness QC table and classify those
    duplicates by gene identifier structure for technical review.

Input files:
    data/interim/orthology/msa/
    msa_input_strict_with_reference.sequence_readiness_qc.tsv

Output files:
    data/interim/orthology/msa/duplicate_species_qc.tsv
    data/interim/orthology/msa/duplicate_species_qc_summary.txt

Pipeline position:
    sequence_readiness_qc
        ->
    duplicate_species_qc

Notes:
    This module performs QC only. It does not filter sequences, redefine
    orthology, or use Ser-rich sequence, N-terminal extension, phosphosite
    motifs, domains, IDR, or trait-level information. The Homo sapiens
    REFERENCE record is ignored for duplicate-species classification.
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


READINESS_QC_TSV = Path(
    "data/interim/orthology/msa/"
    "msa_input_strict_with_reference.sequence_readiness_qc.tsv"
)
OUTPUT_QC_TSV = Path("data/interim/orthology/msa/duplicate_species_qc.tsv")
OUTPUT_SUMMARY_TXT = Path(
    "data/interim/orthology/msa/duplicate_species_qc_summary.txt"
)

REQUIRED_FIELDS = {
    "protein_id",
    "species",
    "gene_id",
    "group",
    "raw_length",
    "tree_readiness",
    "nterm_readiness",
    "ptm_readiness",
    "first_aligned_human_position",
    "nterm_window_coverage",
    "ptm_window_coverage",
    "sequence_qc_flags",
    "sequence_qc_note",
}
OUTPUT_FIELDNAMES = [
    "species",
    "species_record_count",
    "duplicate_class",
    "protein_id",
    "gene_id",
    "group",
    "raw_length",
    "tree_readiness",
    "nterm_readiness",
    "ptm_readiness",
    "first_aligned_human_position",
    "nterm_window_coverage",
    "ptm_window_coverage",
    "sequence_qc_flags",
    "sequence_qc_note",
]

SAME_GENE_MULTIPLE_PROTEINS = "SAME_GENE_MULTIPLE_PROTEINS"
DIFFERENT_GENES_SAME_SPECIES = "DIFFERENT_GENES_SAME_SPECIES"
MISSING_GENE_ID_REVIEW = "MISSING_GENE_ID_REVIEW"


@dataclass(frozen=True)
class ReadinessRecord:
    """
    Fields from one sequence readiness QC row needed for duplicate QC.
    """

    protein_id: str
    species: str
    gene_id: str
    group: str
    raw_length: str
    tree_readiness: str
    nterm_readiness: str
    ptm_readiness: str
    first_aligned_human_position: str
    nterm_window_coverage: str
    ptm_window_coverage: str
    sequence_qc_flags: str
    sequence_qc_note: str


@dataclass(frozen=True)
class DuplicateSpeciesSummary:
    """
    Summary of one species with multiple non-reference records.
    """

    species: str
    species_record_count: int
    duplicate_class: str
    gene_ids: list[str]
    protein_ids: list[str]


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


def load_readiness_records(path: Path) -> list[ReadinessRecord]:
    """
    Load sequence readiness QC records from a TSV file.
    """

    require_input_file(path)
    records = []

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_fields(path, reader.fieldnames, REQUIRED_FIELDS)

        for row in reader:
            protein_id = row["protein_id"].strip()
            species = row["species"].strip()
            if not protein_id:
                raise ValueError(f"Readiness row has empty protein_id: {path}")
            if not species:
                raise ValueError(f"Readiness row has empty species: {path}")

            records.append(
                ReadinessRecord(
                    protein_id=protein_id,
                    species=species,
                    gene_id=row["gene_id"].strip(),
                    group=row["group"].strip(),
                    raw_length=row["raw_length"].strip(),
                    tree_readiness=row["tree_readiness"].strip(),
                    nterm_readiness=row["nterm_readiness"].strip(),
                    ptm_readiness=row["ptm_readiness"].strip(),
                    first_aligned_human_position=row[
                        "first_aligned_human_position"
                    ].strip(),
                    nterm_window_coverage=row["nterm_window_coverage"].strip(),
                    ptm_window_coverage=row["ptm_window_coverage"].strip(),
                    sequence_qc_flags=row["sequence_qc_flags"].strip(),
                    sequence_qc_note=row["sequence_qc_note"].strip(),
                )
            )

    if not records:
        raise ValueError(f"Input TSV has no data records: {path}")

    return records


def is_reference_record(record: ReadinessRecord) -> bool:
    """
    Return whether a record is the coordinate reference sequence.
    """

    return record.group == "REFERENCE"


def group_non_reference_records_by_species(
    records: list[ReadinessRecord],
) -> dict[str, list[ReadinessRecord]]:
    """
    Group non-reference records by species.
    """

    records_by_species: dict[str, list[ReadinessRecord]] = defaultdict(list)
    for record in records:
        if not is_reference_record(record):
            records_by_species[record.species].append(record)
    return dict(records_by_species)


def classify_duplicate_species(records: list[ReadinessRecord]) -> str:
    """
    Classify duplicate records for one species by gene identifier pattern.
    """

    gene_ids = [record.gene_id for record in records]
    if any(not gene_id for gene_id in gene_ids):
        return MISSING_GENE_ID_REVIEW

    unique_gene_ids = set(gene_ids)
    if len(unique_gene_ids) == 1:
        return SAME_GENE_MULTIPLE_PROTEINS

    return DIFFERENT_GENES_SAME_SPECIES


def build_duplicate_summaries(
    records_by_species: dict[str, list[ReadinessRecord]],
) -> list[DuplicateSpeciesSummary]:
    """
    Build one summary object per species with multiple non-reference records.
    """

    summaries = []
    for species in sorted(records_by_species):
        species_records = records_by_species[species]
        if len(species_records) <= 1:
            continue

        gene_ids = sorted(
            {record.gene_id if record.gene_id else "<empty>" for record in species_records}
        )
        protein_ids = sorted(record.protein_id for record in species_records)
        summaries.append(
            DuplicateSpeciesSummary(
                species=species,
                species_record_count=len(species_records),
                duplicate_class=classify_duplicate_species(species_records),
                gene_ids=gene_ids,
                protein_ids=protein_ids,
            )
        )

    return summaries


def build_output_records(
    records_by_species: dict[str, list[ReadinessRecord]],
    summaries: list[DuplicateSpeciesSummary],
) -> list[dict[str, str]]:
    """
    Build one output TSV row per record in duplicated species.
    """

    duplicate_class_by_species = {
        summary.species: summary.duplicate_class
        for summary in summaries
    }
    output_records = []

    for species in sorted(duplicate_class_by_species):
        species_records = sorted(
            records_by_species[species],
            key=lambda record: record.protein_id,
        )
        species_record_count = str(len(species_records))
        duplicate_class = duplicate_class_by_species[species]

        for record in species_records:
            output_records.append(
                {
                    "species": record.species,
                    "species_record_count": species_record_count,
                    "duplicate_class": duplicate_class,
                    "protein_id": record.protein_id,
                    "gene_id": record.gene_id,
                    "group": record.group,
                    "raw_length": record.raw_length,
                    "tree_readiness": record.tree_readiness,
                    "nterm_readiness": record.nterm_readiness,
                    "ptm_readiness": record.ptm_readiness,
                    "first_aligned_human_position": (
                        record.first_aligned_human_position
                    ),
                    "nterm_window_coverage": record.nterm_window_coverage,
                    "ptm_window_coverage": record.ptm_window_coverage,
                    "sequence_qc_flags": record.sequence_qc_flags,
                    "sequence_qc_note": record.sequence_qc_note,
                }
            )

    return output_records


def write_tsv(records: list[dict[str, str]], output_path: Path) -> None:
    """
    Write duplicate-species QC rows to a tab-separated file.
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


def write_summary(
    non_reference_record_count: int,
    species_count: int,
    summaries: list[DuplicateSpeciesSummary],
    output_path: Path,
) -> None:
    """
    Write a text summary of duplicated species and duplicate classes.
    """

    class_counts = Counter(summary.duplicate_class for summary in summaries)
    lines = [
        "PHF10 duplicate species QC",
        "QC only: no orthology redefinition and no biological interpretation.",
        "",
        f"number of non-reference records: {non_reference_record_count}",
        f"number of species: {species_count}",
        f"number of species with multiple records: {len(summaries)}",
        (
            f"{SAME_GENE_MULTIPLE_PROTEINS} species count: "
            f"{class_counts[SAME_GENE_MULTIPLE_PROTEINS]}"
        ),
        (
            f"{DIFFERENT_GENES_SAME_SPECIES} species count: "
            f"{class_counts[DIFFERENT_GENES_SAME_SPECIES]}"
        ),
        (
            f"{MISSING_GENE_ID_REVIEW} species count: "
            f"{class_counts[MISSING_GENE_ID_REVIEW]}"
        ),
        "",
        "Duplicated species:",
        "species\tspecies_record_count\tduplicate_class\tgene_ids\tprotein_ids",
    ]

    for summary in summaries:
        lines.append(
            "\t".join(
                [
                    summary.species,
                    str(summary.species_record_count),
                    summary.duplicate_class,
                    ";".join(summary.gene_ids),
                    ";".join(summary.protein_ids),
                ]
            )
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """
    Run duplicate-species QC and write output files.
    """

    records = load_readiness_records(READINESS_QC_TSV)
    records_by_species = group_non_reference_records_by_species(records)
    non_reference_record_count = sum(
        len(species_records)
        for species_records in records_by_species.values()
    )
    summaries = build_duplicate_summaries(records_by_species)
    output_records = build_output_records(records_by_species, summaries)

    write_tsv(output_records, OUTPUT_QC_TSV)
    write_summary(
        non_reference_record_count=non_reference_record_count,
        species_count=len(records_by_species),
        summaries=summaries,
        output_path=OUTPUT_SUMMARY_TXT,
    )

    print(f"non_reference_records: {non_reference_record_count}")
    print(f"duplicated_species: {len(summaries)}")
    print(f"saved: {OUTPUT_QC_TSV}")
    print(f"saved: {OUTPUT_SUMMARY_TXT}")


if __name__ == "__main__":
    main()

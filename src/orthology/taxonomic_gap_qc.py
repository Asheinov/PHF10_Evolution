"""
Summarize taxonomic gaps in the reference-aware PHF10 readiness dataset.

Purpose:
    Compare the current sequence readiness QC table against a manually
    curated target taxonomic backbone. The output reports expected clades
    and species that are present, missing, or represented by records that
    are weak for existing tree, N-terminal, or PTM readiness fields.

Input files:
    data/interim/orthology/msa/
    msa_input_strict_with_reference.sequence_readiness_qc.tsv

    data/interim/orthology/taxonomy/target_taxonomic_backbone.tsv

Output files:
    data/interim/orthology/taxonomy/taxonomic_gap_qc.tsv
    data/interim/orthology/taxonomy/taxonomic_gap_summary.tsv

Pipeline position:
    sequence_readiness_qc
        ->
    taxonomic_gap_qc

Biological limitation:
    This module is QC-only. It does not infer motif origin, search for
    phosphomodule motifs, redefine orthology, filter proteins, or classify
    sequences by Ser-rich sequence, N-terminal extension, domains, IDR, or
    trait-level information. It only summarizes taxonomic representation
    using existing readiness fields.
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


READINESS_TSV = Path(
    "data/interim/orthology/msa/"
    "msa_input_strict_with_reference.sequence_readiness_qc.tsv"
)
BACKBONE_TSV = Path(
    "data/interim/orthology/taxonomy/target_taxonomic_backbone.tsv"
)
OUTPUT_QC_TSV = Path(
    "data/interim/orthology/taxonomy/taxonomic_gap_qc.tsv"
)
OUTPUT_SUMMARY_TSV = Path(
    "data/interim/orthology/taxonomy/taxonomic_gap_summary.tsv"
)

READINESS_REQUIRED_FIELDS = {
    "protein_id",
    "species",
    "tree_readiness",
    "nterm_readiness",
    "ptm_readiness",
    "sequence_qc_flags",
}
BACKBONE_REQUIRED_FIELDS = {
    "clade_order",
    "clade_id",
    "clade_name",
    "priority",
    "expected_species",
    "why_needed",
}
QC_FIELDNAMES = [
    "clade_order",
    "clade_id",
    "clade_name",
    "priority",
    "expected_species",
    "observed_species",
    "missing_species",
    "expected_species_count",
    "observed_species_count",
    "missing_species_count",
    "total_records",
    "tree_ok_records",
    "nterm_ok_records",
    "ptm_ok_records",
    "nterm_truncated_records",
    "ptm_not_ready_records",
    "gap_status",
    "why_needed",
]
SUMMARY_FIELDNAMES = [
    "gap_status",
    "clade_count",
]


@dataclass(frozen=True)
class ReadinessRecord:
    """
    Existing readiness QC fields needed for taxonomic gap summaries.
    """

    protein_id: str
    species: str
    tree_readiness: str
    nterm_readiness: str
    ptm_readiness: str
    sequence_qc_flags: str


@dataclass(frozen=True)
class BackboneClade:
    """
    One expected taxonomic backbone row.
    """

    clade_order: str
    clade_id: str
    clade_name: str
    priority: str
    expected_species: list[str]
    why_needed: str


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
    Load the reference-aware sequence readiness QC table.
    """

    require_input_file(path)
    records = []

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_fields(path, reader.fieldnames, READINESS_REQUIRED_FIELDS)

        for row in reader:
            species = row["species"].strip()
            protein_id = row["protein_id"].strip()
            if not species:
                raise ValueError(f"Readiness row has empty species: {path}")
            if not protein_id:
                raise ValueError(f"Readiness row has empty protein_id: {path}")

            records.append(
                ReadinessRecord(
                    protein_id=protein_id,
                    species=species,
                    tree_readiness=row["tree_readiness"].strip(),
                    nterm_readiness=row["nterm_readiness"].strip(),
                    ptm_readiness=row["ptm_readiness"].strip(),
                    sequence_qc_flags=row["sequence_qc_flags"].strip(),
                )
            )

    if not records:
        raise ValueError(f"Input TSV has no data records: {path}")

    return records


def load_backbone(path: Path) -> list[BackboneClade]:
    """
    Load the target taxonomic backbone table.
    """

    require_input_file(path)
    clades = []

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_fields(path, reader.fieldnames, BACKBONE_REQUIRED_FIELDS)

        for row in reader:
            expected_species = [
                species.strip()
                for species in row["expected_species"].split(";")
                if species.strip()
            ]
            if not expected_species:
                raise ValueError(
                    f"Backbone row has no expected species: {path}: "
                    f"{row.get('clade_id', '')}"
                )

            clades.append(
                BackboneClade(
                    clade_order=row["clade_order"].strip(),
                    clade_id=row["clade_id"].strip(),
                    clade_name=row["clade_name"].strip(),
                    priority=row["priority"].strip(),
                    expected_species=expected_species,
                    why_needed=row["why_needed"].strip(),
                )
            )

    if not clades:
        raise ValueError(f"Input TSV has no data records: {path}")

    return clades


def index_records_by_species(
    records: list[ReadinessRecord],
) -> dict[str, list[ReadinessRecord]]:
    """
    Group readiness records by species.
    """

    records_by_species: dict[str, list[ReadinessRecord]] = defaultdict(list)
    for record in records:
        records_by_species[record.species].append(record)
    return dict(records_by_species)


def assign_gap_status(
    priority: str,
    observed_species_count: int,
    tree_ok_records: int,
    nterm_ok_records: int,
    ptm_ok_records: int,
    nterm_truncated_records: int,
    ptm_not_ready_records: int,
) -> str:
    """
    Assign the taxonomic representation status for one backbone clade.
    """

    if priority == "CRITICAL" and observed_species_count == 0:
        return "MISSING_CRITICAL"
    if observed_species_count == 0:
        return "MISSING"
    if (
        tree_ok_records > 0
        and nterm_ok_records > 0
        and ptm_ok_records > 0
    ):
        return "PRESENT_STRONG"
    if tree_ok_records > 0 and ptm_ok_records == 0:
        return "PRESENT_TREE_ONLY"
    if nterm_ok_records == 0 and nterm_truncated_records > 0:
        return "PRESENT_BUT_NTERM_WEAK"
    if ptm_ok_records == 0 and ptm_not_ready_records > 0:
        return "PRESENT_BUT_PTM_WEAK"
    return "PRESENT_REVIEW"


def summarize_clade(
    clade: BackboneClade,
    records_by_species: dict[str, list[ReadinessRecord]],
) -> dict[str, str]:
    """
    Summarize one expected backbone clade against the readiness records.
    """

    observed_species = [
        species
        for species in clade.expected_species
        if species in records_by_species
    ]
    missing_species = [
        species
        for species in clade.expected_species
        if species not in records_by_species
    ]
    clade_records = [
        record
        for species in observed_species
        for record in records_by_species[species]
    ]

    tree_ok_records = sum(
        record.tree_readiness == "TREE_OK"
        for record in clade_records
    )
    nterm_ok_records = sum(
        record.nterm_readiness == "NTERM_OK"
        for record in clade_records
    )
    ptm_ok_records = sum(
        record.ptm_readiness == "PTM_OK"
        for record in clade_records
    )
    nterm_truncated_records = sum(
        record.nterm_readiness == "NTERM_TRUNCATED"
        for record in clade_records
    )
    ptm_not_ready_records = sum(
        record.ptm_readiness == "PTM_NOT_READY"
        for record in clade_records
    )

    gap_status = assign_gap_status(
        priority=clade.priority,
        observed_species_count=len(observed_species),
        tree_ok_records=tree_ok_records,
        nterm_ok_records=nterm_ok_records,
        ptm_ok_records=ptm_ok_records,
        nterm_truncated_records=nterm_truncated_records,
        ptm_not_ready_records=ptm_not_ready_records,
    )

    return {
        "clade_order": clade.clade_order,
        "clade_id": clade.clade_id,
        "clade_name": clade.clade_name,
        "priority": clade.priority,
        "expected_species": ";".join(clade.expected_species),
        "observed_species": ";".join(observed_species),
        "missing_species": ";".join(missing_species),
        "expected_species_count": str(len(clade.expected_species)),
        "observed_species_count": str(len(observed_species)),
        "missing_species_count": str(len(missing_species)),
        "total_records": str(len(clade_records)),
        "tree_ok_records": str(tree_ok_records),
        "nterm_ok_records": str(nterm_ok_records),
        "ptm_ok_records": str(ptm_ok_records),
        "nterm_truncated_records": str(nterm_truncated_records),
        "ptm_not_ready_records": str(ptm_not_ready_records),
        "gap_status": gap_status,
        "why_needed": clade.why_needed,
    }


def build_gap_qc_records(
    clades: list[BackboneClade],
    readiness_records: list[ReadinessRecord],
) -> list[dict[str, str]]:
    """
    Build one taxonomic gap QC output row per backbone clade.
    """

    records_by_species = index_records_by_species(readiness_records)
    return [
        summarize_clade(clade, records_by_species)
        for clade in clades
    ]


def build_gap_summary(records: list[dict[str, str]]) -> list[dict[str, str]]:
    """
    Count clades by gap status.
    """

    counts = Counter(record["gap_status"] for record in records)
    return [
        {
            "gap_status": gap_status,
            "clade_count": str(counts[gap_status]),
        }
        for gap_status in sorted(counts)
    ]


def write_tsv(
    records: list[dict[str, str]],
    output_path: Path,
    fieldnames: list[str],
) -> None:
    """
    Write dictionaries to a tab-separated output file.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    """
    Run taxonomic gap QC and write output tables.
    """

    readiness_records = load_readiness_records(READINESS_TSV)
    clades = load_backbone(BACKBONE_TSV)
    gap_records = build_gap_qc_records(clades, readiness_records)
    summary_records = build_gap_summary(gap_records)

    write_tsv(gap_records, OUTPUT_QC_TSV, QC_FIELDNAMES)
    write_tsv(summary_records, OUTPUT_SUMMARY_TSV, SUMMARY_FIELDNAMES)

    status_counts = Counter(record["gap_status"] for record in gap_records)
    missing_critical_count = status_counts["MISSING_CRITICAL"]
    missing_count = status_counts["MISSING"]

    print(f"clades checked: {len(gap_records)}")
    print(f"MISSING_CRITICAL clades: {missing_critical_count}")
    print(f"MISSING clades: {missing_count}")
    print(f"saved: {OUTPUT_QC_TSV}")
    print(f"saved: {OUTPUT_SUMMARY_TSV}")


if __name__ == "__main__":
    main()

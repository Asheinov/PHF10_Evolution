"""
Create a reviewable fetch plan for missing lamprey PHF10 candidates.

Purpose:
    Build an explicit target list for future lamprey PHF10 candidate
    searches from the current taxonomic gap QC result and a manually
    curated target-species table. This module creates a plan only; it does
    not query Ensembl, NCBI, UniProt, OrthoDB, or any other external
    service.

Input files:
    data/interim/orthology/taxonomy/taxonomic_gap_qc.tsv
    data/interim/orthology/taxonomy/missing_lamprey_targets.tsv

Output files:
    data/interim/orthology/taxonomy/missing_lamprey_fetch_plan.tsv

Pipeline position:
    taxonomic_gap_qc
        ->
    missing_lamprey_fetch_plan

Biological limitation:
    This module does not infer orthology, fetch sequences, classify
    proteins, or search by the phosphomodule. It only records which
    lamprey species should be searched later and which sources should be
    considered.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


TAXONOMIC_GAP_QC_TSV = Path(
    "data/interim/orthology/taxonomy/taxonomic_gap_qc.tsv"
)
MISSING_LAMPREY_TARGETS_TSV = Path(
    "data/interim/orthology/taxonomy/missing_lamprey_targets.tsv"
)
OUTPUT_FETCH_PLAN_TSV = Path(
    "data/interim/orthology/taxonomy/missing_lamprey_fetch_plan.tsv"
)

LAMPREY_CLADE_ID = "cyclostomata_lamprey"
REQUIRED_LAMPREY_GAP_STATUS = "MISSING_CRITICAL"
SEARCH_STAGE = "TARGET_DEFINED"
FETCH_STATUS = "NOT_FETCHED"

GAP_QC_REQUIRED_FIELDS = {
    "clade_id",
    "gap_status",
}
TARGET_REQUIRED_FIELDS = {
    "species",
    "taxon_group",
    "priority",
    "preferred_sources",
    "why_needed",
    "manual_note",
}
OUTPUT_FIELDNAMES = [
    "species",
    "taxon_group",
    "priority",
    "gap_status_from_taxonomic_qc",
    "preferred_sources",
    "search_stage",
    "fetch_status",
    "why_needed",
    "manual_note",
]


@dataclass(frozen=True)
class LampreyTarget:
    """
    One manually defined lamprey search target.
    """

    species: str
    taxon_group: str
    priority: str
    preferred_sources: str
    why_needed: str
    manual_note: str


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


def get_lamprey_gap_status(path: Path) -> str:
    """
    Return the upstream taxonomic gap status for the lamprey clade.
    """

    require_input_file(path)

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_fields(path, reader.fieldnames, GAP_QC_REQUIRED_FIELDS)

        matching_rows = [
            row
            for row in reader
            if row["clade_id"].strip() == LAMPREY_CLADE_ID
        ]

    if not matching_rows:
        raise ValueError(
            f"{LAMPREY_CLADE_ID} is absent from taxonomic gap QC: {path}"
        )
    if len(matching_rows) > 1:
        raise ValueError(
            f"{LAMPREY_CLADE_ID} appears more than once in taxonomic gap QC: "
            f"{path}"
        )

    gap_status = matching_rows[0]["gap_status"].strip()
    if gap_status != REQUIRED_LAMPREY_GAP_STATUS:
        raise ValueError(
            f"{LAMPREY_CLADE_ID} gap_status is {gap_status!r}; expected "
            f"{REQUIRED_LAMPREY_GAP_STATUS!r}"
        )

    return gap_status


def load_lamprey_targets(path: Path) -> list[LampreyTarget]:
    """
    Load manually defined lamprey target species.
    """

    require_input_file(path)
    targets = []

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_fields(path, reader.fieldnames, TARGET_REQUIRED_FIELDS)

        for row in reader:
            species = row["species"].strip()
            if not species:
                raise ValueError(f"Lamprey target row has empty species: {path}")

            targets.append(
                LampreyTarget(
                    species=species,
                    taxon_group=row["taxon_group"].strip(),
                    priority=row["priority"].strip(),
                    preferred_sources=row["preferred_sources"].strip(),
                    why_needed=row["why_needed"].strip(),
                    manual_note=row["manual_note"].strip(),
                )
            )

    if not targets:
        raise ValueError(f"Input TSV has no lamprey target rows: {path}")

    return targets


def build_fetch_plan_records(
    targets: list[LampreyTarget],
    gap_status: str,
) -> list[dict[str, str]]:
    """
    Build one fetch-plan row per lamprey target species.
    """

    return [
        {
            "species": target.species,
            "taxon_group": target.taxon_group,
            "priority": target.priority,
            "gap_status_from_taxonomic_qc": gap_status,
            "preferred_sources": target.preferred_sources,
            "search_stage": SEARCH_STAGE,
            "fetch_status": FETCH_STATUS,
            "why_needed": target.why_needed,
            "manual_note": target.manual_note,
        }
        for target in targets
    ]


def write_tsv(records: list[dict[str, str]], output_path: Path) -> None:
    """
    Write the fetch plan as a tab-separated file.
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


def main() -> None:
    """
    Create the missing-lamprey fetch plan.
    """

    gap_status = get_lamprey_gap_status(TAXONOMIC_GAP_QC_TSV)
    targets = load_lamprey_targets(MISSING_LAMPREY_TARGETS_TSV)
    records = build_fetch_plan_records(targets, gap_status)
    write_tsv(records, OUTPUT_FETCH_PLAN_TSV)

    print(f"lamprey targets written: {len(records)}")
    print(f"output: {OUTPUT_FETCH_PLAN_TSV}")


if __name__ == "__main__":
    main()

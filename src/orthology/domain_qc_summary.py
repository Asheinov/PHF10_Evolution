"""
Summarize parsed Pfam/HMMER domain scan results for PHF10 ortholog QC.

Purpose:
    Build per-protein QC summary tables from parsed domain_scan output.
    This module summarizes domain hits for later manual inspection only.

Input files:
    data/interim/orthology/core_orthologs.fasta
    data/interim/orthology/extended_orthologs.fasta
    data/interim/orthology/domain_scan/core_domains.tsv
    data/interim/orthology/domain_scan/extended_domains.tsv

Output files:
    data/interim/orthology/domain_scan/core_domain_qc_summary.tsv
    data/interim/orthology/domain_scan/extended_domain_qc_summary.tsv
    data/interim/orthology/domain_scan/domain_qc_report.txt

Pipeline position:
    parse_domain_scan
        ->
    domain_qc_summary
        ->
    manual QC / later orthology QC

Notes:
    This module does not filter proteins, define orthology, collapse
    domains, or perform biological interpretation. Domain scanning is
    QC only and must not define PHF10 orthology.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


CORE_FASTA = Path("data/interim/orthology/core_orthologs.fasta")
EXTENDED_FASTA = Path("data/interim/orthology/extended_orthologs.fasta")
CORE_DOMAINS = Path("data/interim/orthology/domain_scan/core_domains.tsv")
EXTENDED_DOMAINS = Path(
    "data/interim/orthology/domain_scan/extended_domains.tsv"
)
CORE_OUTPUT = Path(
    "data/interim/orthology/domain_scan/core_domain_qc_summary.tsv"
)
EXTENDED_OUTPUT = Path(
    "data/interim/orthology/domain_scan/extended_domain_qc_summary.tsv"
)
REPORT_OUTPUT = Path(
    "data/interim/orthology/domain_scan/domain_qc_report.txt"
)

OUTPUT_COLUMNS = [
    "protein_id",
    "species",
    "gene_id",
    "group",
    "protein_length",
    "n_domain_hits",
    "n_unique_pfam_ids",
    "unique_pfam_ids",
    "unique_domain_names",
    "has_any_domain",
    "has_phd_domain",
    "has_ini1_dna_bd",
    "min_domain_start",
    "max_domain_end",
]


def parse_fasta(fasta_path: str) -> dict[str, dict[str, Any]]:
    """
    Load FASTA records keyed by protein ID and calculate sequence length.
    """

    path = Path(fasta_path)

    if not path.exists():
        raise FileNotFoundError(f"Required FASTA input is missing: {path}")

    records: dict[str, dict[str, Any]] = {}
    header: str | None = None
    sequence_parts: list[str] = []

    def store_record(record_header: str, sequence: list[str]) -> None:
        parts = record_header.split("|")

        if len(parts) != 3:
            raise ValueError(
                "Expected FASTA header format "
                f"protein_id|species|gene_id, got: {record_header}"
            )

        protein_id, species, gene_id = parts

        if protein_id in records:
            raise ValueError(f"Duplicate FASTA protein ID: {protein_id}")

        records[protein_id] = {
            "protein_id": protein_id,
            "species": species,
            "gene_id": gene_id,
            "protein_length": len("".join(sequence)),
        }

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

    return records


def load_domain_table(tsv_path: str) -> list[dict[str, str]]:
    """
    Load parsed domain hit rows from a TSV file.
    """

    path = Path(tsv_path)

    if not path.exists():
        raise FileNotFoundError(
            f"Required parsed domain table is missing: {path}"
        )

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def summarize_group(
    group: str,
    fasta_path: str,
    domains_tsv: str,
) -> list[dict[str, Any]]:
    """
    Summarize parsed domain hits for every protein in one FASTA group.
    """

    fasta_records = parse_fasta(fasta_path)
    domain_records = load_domain_table(domains_tsv)
    domains_by_protein: dict[str, list[dict[str, str]]] = {}

    for domain in domain_records:
        protein_id = domain.get("protein_id", "")
        domains_by_protein.setdefault(protein_id, []).append(domain)

    summaries = []

    for protein_id in sorted(fasta_records):
        fasta_record = fasta_records[protein_id]
        domains = domains_by_protein.get(protein_id, [])

        pfam_ids = sorted(
            {
                domain["pfam_id"]
                for domain in domains
                if domain.get("pfam_id")
            }
        )
        domain_names = sorted(
            {
                domain["domain_name"]
                for domain in domains
                if domain.get("domain_name")
            }
        )
        env_starts = [
            int(domain["env_start"])
            for domain in domains
            if domain.get("env_start")
        ]
        env_ends = [
            int(domain["env_end"])
            for domain in domains
            if domain.get("env_end")
        ]

        has_phd_domain = any(
            "phd" in domain_name.lower()
            for domain_name in domain_names
        )
        has_ini1_dna_bd = any(
            "ini1_dna-bd" in domain_name.lower()
            for domain_name in domain_names
        )

        summaries.append(
            {
                "protein_id": protein_id,
                "species": fasta_record["species"],
                "gene_id": fasta_record["gene_id"],
                "group": group,
                "protein_length": fasta_record["protein_length"],
                "n_domain_hits": len(domains),
                "n_unique_pfam_ids": len(pfam_ids),
                "unique_pfam_ids": ";".join(pfam_ids),
                "unique_domain_names": ";".join(domain_names),
                "has_any_domain": len(domains) > 0,
                "has_phd_domain": has_phd_domain,
                "has_ini1_dna_bd": has_ini1_dna_bd,
                "min_domain_start": min(env_starts) if env_starts else "",
                "max_domain_end": max(env_ends) if env_ends else "",
            }
        )

    return summaries


def write_tsv(records: list[dict[str, Any]], output_path: str) -> None:
    """
    Write per-protein QC summary records to a TSV file.
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


def count_true(records: list[dict[str, Any]], field: str) -> int:
    """
    Count records where a boolean summary field is true.
    """

    return sum(1 for record in records if record[field])


def proteins_without_hits(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """
    Return summary records for proteins without parsed domain hits.
    """

    return [
        record
        for record in records
        if not record["has_any_domain"]
    ]


def write_report(
    core_records: list[dict[str, Any]],
    extended_records: list[dict[str, Any]],
    output_path: str,
) -> None:
    """
    Write a plain-text QC report for domain summary tables.
    """

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    core_without_hits = proteins_without_hits(core_records)
    extended_without_hits = proteins_without_hits(extended_records)

    lines = [
        "PHF10 domain QC summary",
        "QC only: domain hits are not orthology criteria.",
        "",
        f"CORE FASTA proteins: {len(core_records)}",
        f"EXTENDED FASTA proteins: {len(extended_records)}",
        "",
        "Proteins with at least one domain hit:",
        f"CORE: {count_true(core_records, 'has_any_domain')}",
        f"EXTENDED: {count_true(extended_records, 'has_any_domain')}",
        "",
        "Proteins without any domain hit:",
        f"CORE: {len(core_without_hits)}",
        f"EXTENDED: {len(extended_without_hits)}",
        "",
        "CORE proteins without domain hits:",
    ]

    if core_without_hits:
        for record in core_without_hits:
            lines.append(
                "\t".join(
                    [
                        str(record["protein_id"]),
                        str(record["species"]),
                        str(record["gene_id"]),
                    ]
                )
            )
    else:
        lines.append("None")

    lines.extend(
        [
            "",
            "EXTENDED proteins without domain hits:",
        ]
    )

    if extended_without_hits:
        for record in extended_without_hits:
            lines.append(
                "\t".join(
                    [
                        str(record["protein_id"]),
                        str(record["species"]),
                        str(record["gene_id"]),
                    ]
                )
            )
    else:
        lines.append("None")

    lines.extend(
        [
            "",
            "Proteins with PHD-like hits:",
            f"CORE: {count_true(core_records, 'has_phd_domain')}",
            f"EXTENDED: {count_true(extended_records, 'has_phd_domain')}",
            "",
            "Proteins with INI1_DNA-bd-like hits:",
            f"CORE: {count_true(core_records, 'has_ini1_dna_bd')}",
            f"EXTENDED: {count_true(extended_records, 'has_ini1_dna_bd')}",
            "",
        ]
    )

    path.write_text("\n".join(lines), encoding="utf-8")


def run() -> None:
    """
    Create domain QC summary tables and a text report.
    """

    core_records = summarize_group(
        group="CORE",
        fasta_path=str(CORE_FASTA),
        domains_tsv=str(CORE_DOMAINS),
    )
    extended_records = summarize_group(
        group="EXTENDED",
        fasta_path=str(EXTENDED_FASTA),
        domains_tsv=str(EXTENDED_DOMAINS),
    )

    write_tsv(core_records, str(CORE_OUTPUT))
    write_tsv(extended_records, str(EXTENDED_OUTPUT))
    write_report(core_records, extended_records, str(REPORT_OUTPUT))

    print(
        f"CORE: {len(core_records)} proteins; "
        f"{count_true(core_records, 'has_any_domain')} with domain hits"
    )
    print(
        f"EXTENDED: {len(extended_records)} proteins; "
        f"{count_true(extended_records, 'has_any_domain')} with domain hits"
    )
    print(f"Saved report: {REPORT_OUTPUT}")


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

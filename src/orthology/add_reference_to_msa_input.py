"""
Add Homo sapiens PHF10 reference sequence to strict MSA input.

Purpose
-------
This module creates a separate MSA input dataset that contains the strict
ortholog set plus the human PHF10 reference sequence. The reference is needed
only for coordinate mapping in downstream sequence/readiness QC.

Inputs
------
- data/interim/orthology/msa_input_strict.fasta
- data/interim/orthology/msa_input_strict.tsv
- data/raw/Q8WUB8.fasta

Outputs
-------
- data/interim/orthology/msa_input_strict_with_reference.fasta
- data/interim/orthology/msa_input_strict_with_reference.tsv

Pipeline position
-----------------
Runs after:
- src/orthology/build_msa_input.py

Runs before:
- MAFFT alignment with reference
- src/orthology/sequence_readiness_qc.py

Notes
-----
This module does not modify the original strict MSA input files.
It does not filter sequences, define orthology, or interpret biology.
The Homo sapiens sequence is added as a REFERENCE record, not as an ortholog.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


STRICT_FASTA = Path("data/interim/orthology/msa_input_strict.fasta")
STRICT_TSV = Path("data/interim/orthology/msa_input_strict.tsv")
REFERENCE_FASTA = Path("data/raw/Q8WUB8.fasta")

OUTPUT_FASTA = Path("data/interim/orthology/msa_input_strict_with_reference.fasta")
OUTPUT_TSV = Path("data/interim/orthology/msa_input_strict_with_reference.tsv")

REFERENCE_PROTEIN_ID = "Q8WUB8"
REFERENCE_SPECIES = "homo_sapiens"
REFERENCE_GENE_ID = "PHF10"
REFERENCE_GROUP = "REFERENCE"
REFERENCE_HEADER = "Q8WUB8|homo_sapiens|PHF10|REFERENCE"
REFERENCE_NOTE = "Human PHF10 reference sequence for coordinate mapping"


@dataclass(frozen=True)
class FastaRecord:
    """Simple FASTA record."""

    header: str
    sequence: str

    @property
    def protein_id(self) -> str:
        """Return protein_id as the first pipe-delimited header field."""
        return self.header.split("|", 1)[0]


def require_file(path: Path) -> None:
    """Fail clearly if an input file does not exist."""
    if not path.exists():
        raise FileNotFoundError(f"Required input file not found: {path}")


def read_fasta(path: Path) -> list[FastaRecord]:
    """Read FASTA records without external libraries."""
    require_file(path)

    records: list[FastaRecord] = []
    current_header: str | None = None
    current_sequence: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    records.append(
                        FastaRecord(
                            header=current_header,
                            sequence="".join(current_sequence),
                        )
                    )
                current_header = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line)

    if current_header is not None:
        records.append(
            FastaRecord(
                header=current_header,
                sequence="".join(current_sequence),
            )
        )

    if not records:
        raise ValueError(f"No FASTA records found in: {path}")

    return records


def write_fasta(records: list[FastaRecord], path: Path) -> None:
    """Write FASTA records with wrapped sequence lines."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.header}\n")
            for index in range(0, len(record.sequence), 80):
                handle.write(record.sequence[index : index + 80] + "\n")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read TSV rows and preserve fieldnames."""
    require_file(path)

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")

        rows = list(reader)

    return reader.fieldnames, rows


def write_tsv(fieldnames: list[str], rows: list[dict[str, str]], path: Path) -> None:
    """Write TSV rows with preserved field order."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def get_single_reference_sequence(path: Path) -> str:
    """Read Q8WUB8 reference FASTA and return exactly one sequence."""
    records = read_fasta(path)

    if len(records) != 1:
        raise ValueError(
            f"Expected exactly one reference FASTA record in {path}, "
            f"found {len(records)}"
        )

    sequence = records[0].sequence.strip()

    if not sequence:
        raise ValueError(f"Reference sequence is empty: {path}")

    return sequence


def validate_no_reference_already_present(
    fasta_records: list[FastaRecord],
    metadata_rows: list[dict[str, str]],
) -> None:
    """Prevent accidental duplicate Homo sapiens reference records."""
    fasta_ids = {record.protein_id for record in fasta_records}
    metadata_ids = {row["protein_id"] for row in metadata_rows}

    if REFERENCE_PROTEIN_ID in fasta_ids:
        raise ValueError(
            f"Reference protein_id already present in FASTA: {REFERENCE_PROTEIN_ID}"
        )

    if REFERENCE_PROTEIN_ID in metadata_ids:
        raise ValueError(
            f"Reference protein_id already present in TSV: {REFERENCE_PROTEIN_ID}"
        )

    human_rows = [
        row for row in metadata_rows if row.get("species") == REFERENCE_SPECIES
    ]
    if human_rows:
        raise ValueError(
            f"Input TSV already contains {REFERENCE_SPECIES} rows: "
            f"{len(human_rows)}"
        )


def build_reference_metadata_row(
    fieldnames: list[str],
    reference_length: int,
) -> dict[str, str]:
    """Create metadata row for the human PHF10 reference."""
    row = {fieldname: "" for fieldname in fieldnames}

    required_fields = [
        "protein_id",
        "species",
        "gene_id",
        "group",
        "protein_length",
        "manual_status",
        "manual_note",
    ]

    missing_fields = [field for field in required_fields if field not in fieldnames]
    if missing_fields:
        raise ValueError(
            "Input TSV is missing required fields: "
            + ", ".join(missing_fields)
        )

    row["protein_id"] = REFERENCE_PROTEIN_ID
    row["species"] = REFERENCE_SPECIES
    row["gene_id"] = REFERENCE_GENE_ID
    row["group"] = REFERENCE_GROUP
    row["protein_length"] = str(reference_length)
    row["manual_status"] = ""
    row["manual_note"] = REFERENCE_NOTE

    return row


def main() -> None:
    """Create strict MSA input with Homo sapiens PHF10 reference."""
    strict_records = read_fasta(STRICT_FASTA)
    fieldnames, metadata_rows = read_tsv(STRICT_TSV)
    reference_sequence = get_single_reference_sequence(REFERENCE_FASTA)

    validate_no_reference_already_present(strict_records, metadata_rows)

    reference_record = FastaRecord(
        header=REFERENCE_HEADER,
        sequence=reference_sequence,
    )
    reference_metadata = build_reference_metadata_row(
        fieldnames=fieldnames,
        reference_length=len(reference_sequence),
    )

    output_records = strict_records + [reference_record]
    output_rows = metadata_rows + [reference_metadata]

    write_fasta(output_records, OUTPUT_FASTA)
    write_tsv(fieldnames, output_rows, OUTPUT_TSV)

    print(f"ortholog records: {len(strict_records)}")
    print("reference records: 1")
    print(f"total FASTA records: {len(output_records)}")
    print(f"total TSV rows: {len(output_rows)}")
    print(f"reference length: {len(reference_sequence)}")
    print(f"saved: {OUTPUT_FASTA}")
    print(f"saved: {OUTPUT_TSV}")


if __name__ == "__main__":
    main()

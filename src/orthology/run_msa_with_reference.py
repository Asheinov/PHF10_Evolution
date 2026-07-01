"""
Run MAFFT on strict PHF10 MSA input plus Homo sapiens reference.

Purpose:
    Produce a reproducible MAFFT alignment from the strict PHF10
    ortholog FASTA after adding the Homo sapiens reference sequence.

Input files:
    data/interim/orthology/msa_input_strict_with_reference.fasta

Output files:
    data/interim/orthology/msa/msa_input_strict_with_reference.mafft.fasta
    data/interim/orthology/msa/msa_input_strict_with_reference.mafft.log
    data/interim/orthology/msa/msa_input_strict_with_reference.mafft.run.tsv

Pipeline position:
    add_reference_to_msa_input
        ->
    run_msa_with_reference
        ->
    reference-aware MSA QC

Notes:
    This module runs MAFFT only. It does not filter sequences, use
    domains, IDR, Ser-rich regions, N-terminal extensions, PTM motifs,
    or perform biological interpretation.
"""

from __future__ import annotations

import csv
from pathlib import Path
import shutil
import subprocess


INPUT_FASTA = Path(
    "data/interim/orthology/msa_input_strict_with_reference.fasta"
)
OUTPUT_DIR = Path("data/interim/orthology/msa")
OUTPUT_ALIGNMENT = OUTPUT_DIR / "msa_input_strict_with_reference.mafft.fasta"
LOG_FILE = OUTPUT_DIR / "msa_input_strict_with_reference.mafft.log"
RUN_METADATA = OUTPUT_DIR / "msa_input_strict_with_reference.mafft.run.tsv"
MAFFT_COMMAND = [
    "mafft",
    "--auto",
    str(INPUT_FASTA),
]


def require_input_file(path: Path) -> None:
    """
    Fail clearly if the required input FASTA is missing.
    """

    if path.exists():
        return

    raise FileNotFoundError(f"Required input FASTA is missing: {path}")


def require_mafft() -> None:
    """
    Fail clearly if MAFFT is not available on PATH.
    """

    if shutil.which("mafft"):
        return

    raise RuntimeError(
        "MAFFT executable was not found on PATH. Install MAFFT or add it "
        "to PATH before running MSA."
    )


def count_fasta_records(path: Path) -> int:
    """
    Count FASTA records by header lines starting with '>'.
    """

    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.startswith(">"))


def get_mafft_version() -> str:
    """
    Return the MAFFT version string reported by the executable.
    """

    result = subprocess.run(
        ["mafft", "--version"],
        capture_output=True,
        text=True,
    )
    version_text = (result.stdout + result.stderr).strip()

    if not version_text:
        return "unknown"

    return version_text.splitlines()[0]


def run_mafft() -> int:
    """
    Run MAFFT and write captured stdout/stderr to output files.
    """

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        MAFFT_COMMAND,
        capture_output=True,
        text=True,
    )

    OUTPUT_ALIGNMENT.write_text(result.stdout, encoding="utf-8")
    LOG_FILE.write_text(result.stderr, encoding="utf-8")

    return result.returncode


def write_run_metadata(
    mafft_version: str,
    input_records: int,
    output_records: int,
    exit_code: int,
) -> None:
    """
    Write reproducibility metadata for the MAFFT run.
    """

    RUN_METADATA.parent.mkdir(parents=True, exist_ok=True)

    metadata = [
        ("input_fasta", str(INPUT_FASTA)),
        ("output_alignment", str(OUTPUT_ALIGNMENT)),
        ("log_file", str(LOG_FILE)),
        ("command", " ".join(MAFFT_COMMAND)),
        ("mafft_version", mafft_version),
        ("input_records", str(input_records)),
        ("output_records", str(output_records)),
        ("exit_code", str(exit_code)),
    ]

    with RUN_METADATA.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["key", "value"])
        writer.writerows(metadata)


def run() -> None:
    """
    Run MAFFT and save alignment, log, and metadata outputs.
    """

    require_input_file(INPUT_FASTA)
    require_mafft()

    input_records = count_fasta_records(INPUT_FASTA)
    mafft_version = get_mafft_version()
    exit_code = run_mafft()
    output_records = count_fasta_records(OUTPUT_ALIGNMENT)

    write_run_metadata(
        mafft_version=mafft_version,
        input_records=input_records,
        output_records=output_records,
        exit_code=exit_code,
    )

    if exit_code != 0:
        raise RuntimeError(
            f"MAFFT failed with exit code {exit_code}. See log: {LOG_FILE}"
        )

    print(f"input_records: {input_records}")
    print(f"output_records: {output_records}")
    print(f"alignment: {OUTPUT_ALIGNMENT}")
    print(f"log: {LOG_FILE}")
    print(f"metadata: {RUN_METADATA}")


def main() -> int:
    """
    Command-line entry point with concise error reporting.
    """

    try:
        run()
    except (FileNotFoundError, RuntimeError) as exc:
        print(f"ERROR: {exc}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

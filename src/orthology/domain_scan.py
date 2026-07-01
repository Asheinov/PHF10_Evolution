"""
Run HMMER hmmscan against Pfam for PHF10 ortholog domain QC.

Purpose:
    Run hmmscan on the core and extended PHF10 ortholog FASTA files.
    This module only creates HMMER domtblout files; it does not parse
    hits or perform biological interpretation.

Input files:
    data/interim/orthology/core_orthologs.fasta
    data/interim/orthology/extended_orthologs.fasta
    data/external/pfam/Pfam-A.hmm

Output files:
    data/interim/orthology/domain_scan/core.domtblout
    data/interim/orthology/domain_scan/extended.domtblout

Pipeline position:
    architecture_split
        ->
    domain_scan
        ->
    orthology_qc

Notes:
    Domain scanning is QC only. Domain presence must not be used to
    define PHF10 orthology.
"""

from __future__ import annotations

from pathlib import Path
import shutil
import subprocess


CORE_FASTA = Path("data/interim/orthology/core_orthologs.fasta")
EXTENDED_FASTA = Path("data/interim/orthology/extended_orthologs.fasta")
PFAM_HMM = Path("data/external/pfam/Pfam-A.hmm")
OUTPUT_DIR = Path("data/interim/orthology/domain_scan")
CORE_DOMTBLOUT = OUTPUT_DIR / "core.domtblout"
EXTENDED_DOMTBLOUT = OUTPUT_DIR / "extended.domtblout"


def require_hmmscan() -> str:
    """
    Return the hmmscan executable path or raise a clear setup error.
    """

    executable = shutil.which("hmmscan")

    if executable:
        return executable

    raise RuntimeError(
        "hmmscan was not found on PATH. Install HMMER before running "
        "domain scanning."
    )


def require_file(path: Path, label: str) -> None:
    """
    Raise a clear setup error if a required input file is missing.
    """

    if path.exists():
        return

    raise RuntimeError(
        f"Required {label} is missing: {path}"
    )


def run_hmmscan(
    input_fasta: str | Path,
    pfam_hmm: str | Path,
    output_domtblout: str | Path,
) -> None:
    """
    Run hmmscan for one FASTA file and save HMMER domtblout output.
    """

    hmmscan = require_hmmscan()
    input_path = Path(input_fasta)
    pfam_path = Path(pfam_hmm)
    output_path = Path(output_domtblout)

    require_file(input_path, "input FASTA")
    require_file(pfam_path, "Pfam HMM database")

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    command = [
        hmmscan,
        "--noali",
        "--domtblout",
        str(output_path),
        str(pfam_path),
        str(input_path),
    ]

    subprocess.run(
        command,
        check=True,
    )


def scan_core() -> None:
    """
    Run hmmscan for core PHF10 ortholog candidates.
    """

    run_hmmscan(
        input_fasta=CORE_FASTA,
        pfam_hmm=PFAM_HMM,
        output_domtblout=CORE_DOMTBLOUT,
    )


def scan_extended() -> None:
    """
    Run hmmscan for extended PHF10 ortholog candidates.
    """

    run_hmmscan(
        input_fasta=EXTENDED_FASTA,
        pfam_hmm=PFAM_HMM,
        output_domtblout=EXTENDED_DOMTBLOUT,
    )


def main() -> int:
    """
    Run hmmscan for both architecture groups.
    """

    try:
        scan_core()
        scan_extended()
    except RuntimeError as exc:
        print(f"ERROR: {exc}")
        return 1
    except subprocess.CalledProcessError as exc:
        command = " ".join(str(part) for part in exc.cmd)
        print(
            "ERROR: hmmscan failed with "
            f"exit code {exc.returncode}: {command}"
        )
        return exc.returncode

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

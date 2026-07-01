"""
Check external domain scanning resources for PHF10 ortholog QC.

Purpose:
    Report whether the local HMMER/Pfam resources needed by
    domain_scan are available.

Input files:
    data/external/pfam/Pfam-A.hmm
    data/external/pfam/Pfam-A.hmm.h3f
    data/external/pfam/Pfam-A.hmm.h3i
    data/external/pfam/Pfam-A.hmm.h3m
    data/external/pfam/Pfam-A.hmm.h3p

Output files:
    None. This module prints a console report only.

Pipeline position:
    Resource setup check before domain_scan.

Notes:
    Domain scanning is QC only. Domain presence must not be used to
    define PHF10 orthology.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil


PFAM_HMM = Path("data/external/pfam/Pfam-A.hmm")
PRESSED_PFAM_FILES = (
    Path("data/external/pfam/Pfam-A.hmm.h3f"),
    Path("data/external/pfam/Pfam-A.hmm.h3i"),
    Path("data/external/pfam/Pfam-A.hmm.h3m"),
    Path("data/external/pfam/Pfam-A.hmm.h3p"),
)
HMMER_EXECUTABLES = (
    "hmmscan",
    "hmmpress",
    "hmmsearch",
)


@dataclass(frozen=True)
class ResourceStatus:
    """
    Status for one required executable or file.
    """

    name: str
    ok: bool
    detail: str


def check_executable(executable: str) -> ResourceStatus:
    """
    Check whether an executable is available on PATH.
    """

    resolved_path = shutil.which(executable)

    if resolved_path:
        return ResourceStatus(
            name=executable,
            ok=True,
            detail=resolved_path,
        )

    return ResourceStatus(
        name=executable,
        ok=False,
        detail="not found on PATH",
    )


def check_file(path: Path) -> ResourceStatus:
    """
    Check whether a required local file exists.
    """

    if path.exists():
        return ResourceStatus(
            name=str(path),
            ok=True,
            detail="found",
        )

    return ResourceStatus(
        name=str(path),
        ok=False,
        detail="missing",
    )


def check_domain_resources() -> list[ResourceStatus]:
    """
    Check HMMER executables, Pfam-A.hmm, and pressed Pfam files.
    """

    statuses = [
        check_executable(executable)
        for executable in HMMER_EXECUTABLES
    ]

    statuses.append(check_file(PFAM_HMM))
    statuses.extend(
        check_file(path)
        for path in PRESSED_PFAM_FILES
    )

    return statuses


def print_status_report(statuses: list[ResourceStatus]) -> None:
    """
    Print a clear console report for domain scanning prerequisites.
    """

    print("Domain scanning resource check")
    print("QC only: domain presence is not an orthology criterion.")
    print("")

    for status in statuses:
        marker = "OK" if status.ok else "MISSING"
        print(f"[{marker}] {status.name}: {status.detail}")

    print("")

    if all(status.ok for status in statuses):
        print("All required domain scanning resources are available.")
    else:
        print("Missing required resources. Domain scanning is not ready.")


def main() -> int:
    """
    Run the resource check and return a shell exit code.
    """

    statuses = check_domain_resources()
    print_status_report(statuses)

    if all(status.ok for status in statuses):
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())

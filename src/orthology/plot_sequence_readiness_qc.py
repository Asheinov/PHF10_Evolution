"""
Create sequence readiness QC figures for the reference-aware PHF10 MSA.

Purpose:
    Plot technical sequence/readiness QC metrics from the reference-aware
    sequence readiness QC table.

Input files:
    data/interim/orthology/msa/
    msa_input_strict_with_reference.sequence_readiness_qc.tsv

Output files:
    results/orthology/figures/nterm_coverage_all_records.png
    results/orthology/figures/ptm_readiness_counts.png
    results/orthology/figures/species_duplication_counts.png
    results/orthology/figures/x_in_nterm_ptm_window.png

Pipeline position:
    sequence_readiness_qc
        ->
    plot_sequence_readiness_qc

Notes:
    This module plots sequence/readiness QC metrics already present in
    the TSV. It does not filter records, redefine orthology, use
    domains, IDR, Ser-rich regions, N-terminal extensions, PTM motifs,
    or perform biological interpretation.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt


INPUT_QC_TSV = Path(
    "data/interim/orthology/msa/"
    "msa_input_strict_with_reference.sequence_readiness_qc.tsv"
)
OUTPUT_DIR = Path("results/orthology/figures")
NTERM_COVERAGE_PLOT = OUTPUT_DIR / "nterm_coverage_all_records.png"
PTM_READINESS_COUNTS = OUTPUT_DIR / "ptm_readiness_counts.png"
SPECIES_DUPLICATION_COUNTS = OUTPUT_DIR / "species_duplication_counts.png"
X_WINDOW_COUNTS = OUTPUT_DIR / "x_in_nterm_ptm_window.png"


@dataclass(frozen=True)
class ReadinessQcRecord:
    """
    Plotting fields from one sequence readiness QC row.
    """

    protein_id: str
    species: str
    nterm_window_coverage: float
    ptm_readiness: str
    species_record_count: int
    x_count_nterm_window: int
    x_count_ptm_window: int

    @property
    def label(self) -> str:
        """
        Return a compact record label for plots.
        """

        return f"{self.species} {self.protein_id}"


def load_qc_records(path: Path) -> list[ReadinessQcRecord]:
    """
    Load readiness QC records from a TSV file.
    """

    if not path.exists():
        raise FileNotFoundError(f"Required input file is missing: {path}")

    records = []

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Input TSV has no header: {path}")

        for row in reader:
            records.append(
                ReadinessQcRecord(
                    protein_id=row.get("protein_id", ""),
                    species=row.get("species", ""),
                    nterm_window_coverage=float(
                        row.get("nterm_window_coverage", "0")
                    ),
                    ptm_readiness=row.get("ptm_readiness", ""),
                    species_record_count=int(
                        row.get("species_record_count", "0")
                    ),
                    x_count_nterm_window=int(
                        row.get("x_count_nterm_window", "0")
                    ),
                    x_count_ptm_window=int(
                        row.get("x_count_ptm_window", "0")
                    ),
                )
            )

    if not records:
        raise ValueError(f"Input TSV has no data records: {path}")

    return records


def save_nterm_coverage_plot(
    records: list[ReadinessQcRecord],
    output_path: Path,
) -> None:
    """
    Plot N-terminal window coverage across all records.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)
    values = [
        record.nterm_window_coverage
        for record in records
    ]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(
        range(1, len(values) + 1),
        values,
        marker="o",
        linewidth=1,
        markersize=3,
        color="#4C78A8",
    )
    ax.set_xlabel("record index")
    ax.set_ylabel("nterm_window_coverage")
    ax.set_title("N-terminal coverage across all records")
    ax.set_ylim(0, 1.05)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def save_ptm_readiness_counts(
    records: list[ReadinessQcRecord],
    output_path: Path,
) -> None:
    """
    Plot PTM readiness status counts.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter(record.ptm_readiness for record in records)
    labels = sorted(counts)
    values = [counts[label] for label in labels]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(
        labels,
        values,
        color="#59A14F",
        edgecolor="black",
    )
    ax.set_xlabel("ptm_readiness")
    ax.set_ylabel("number of records")
    ax.set_title("PTM readiness status counts")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def save_species_duplication_counts(
    records: list[ReadinessQcRecord],
    output_path: Path,
) -> None:
    """
    Plot species with more than one sequence record.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)
    species_counts = {
        record.species: record.species_record_count
        for record in records
        if record.species_record_count > 1
    }
    sorted_items = sorted(
        species_counts.items(),
        key=lambda item: (-item[1], item[0]),
    )
    labels = [
        item[0]
        for item in reversed(sorted_items)
    ]
    values = [
        item[1]
        for item in reversed(sorted_items)
    ]
    height = max(5, 0.35 * max(1, len(labels)))

    fig, ax = plt.subplots(figsize=(9, height))
    ax.barh(
        labels,
        values,
        color="#F28E2B",
        edgecolor="black",
    )
    ax.set_xlabel("sequence records per species")
    ax.set_ylabel("species")
    ax.set_title("Species with more than one sequence record")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def save_x_window_plot(
    records: list[ReadinessQcRecord],
    output_path: Path,
) -> None:
    """
    Plot records with X symbols in N-terminal and/or PTM windows.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)
    x_records = [
        record
        for record in records
        if (
            record.x_count_nterm_window > 0
            or record.x_count_ptm_window > 0
        )
    ]
    x_records = sorted(
        x_records,
        key=lambda record: (
            -(record.x_count_nterm_window + record.x_count_ptm_window),
            record.protein_id,
        ),
    )
    plot_records = list(reversed(x_records))
    labels = [record.label for record in plot_records]
    nterm_values = [
        record.x_count_nterm_window
        for record in plot_records
    ]
    ptm_values = [
        record.x_count_ptm_window
        for record in plot_records
    ]
    y_positions = list(range(len(plot_records)))
    height = max(5, 0.35 * max(1, len(labels)))

    fig, ax = plt.subplots(figsize=(10, height))
    ax.barh(
        y_positions,
        nterm_values,
        color="#4C78A8",
        edgecolor="black",
        label="N-terminal window",
    )
    ax.barh(
        y_positions,
        ptm_values,
        left=nterm_values,
        color="#E15759",
        edgecolor="black",
        label="PTM window",
    )
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels)
    ax.set_xlabel("X count")
    ax.set_title("Records with X symbols in N-terminal/PTM windows")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def run() -> None:
    """
    Create all sequence readiness QC figures.
    """

    records = load_qc_records(INPUT_QC_TSV)

    save_nterm_coverage_plot(records, NTERM_COVERAGE_PLOT)
    save_ptm_readiness_counts(records, PTM_READINESS_COUNTS)
    save_species_duplication_counts(records, SPECIES_DUPLICATION_COUNTS)
    save_x_window_plot(records, X_WINDOW_COUNTS)

    print(f"records read: {len(records)}")
    print(f"saved: {NTERM_COVERAGE_PLOT}")
    print(f"saved: {PTM_READINESS_COUNTS}")
    print(f"saved: {SPECIES_DUPLICATION_COUNTS}")
    print(f"saved: {X_WINDOW_COUNTS}")


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

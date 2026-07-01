"""
Create technical QC figures for the strict PHF10 MAFFT alignment.

Purpose:
    Plot gap-fraction distributions and top-gap sequences from the
    strict MAFFT alignment QC table.

Input files:
    data/interim/orthology/msa/msa_input_strict.alignment_qc.tsv

Output files:
    results/orthology/figures/strict_msa_gap_fraction_hist.png
    results/orthology/figures/strict_msa_internal_gap_fraction_hist.png
    results/orthology/figures/strict_msa_top_gap_fraction.png
    results/orthology/figures/strict_msa_top_internal_gap_fraction.png

Pipeline position:
    alignment_qc
        ->
    plot_alignment_qc

Notes:
    This module plots technical QC metrics only. It does not filter
    sequences, modify QC data, make inclusion/exclusion decisions, use
    domains, IDR, Ser-rich regions, N-terminal extensions, or perform
    biological interpretation.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt


INPUT_QC_TSV = Path(
    "data/interim/orthology/msa/msa_input_strict.alignment_qc.tsv"
)
OUTPUT_DIR = Path("results/orthology/figures")
GAP_FRACTION_HIST = OUTPUT_DIR / "strict_msa_gap_fraction_hist.png"
INTERNAL_GAP_FRACTION_HIST = (
    OUTPUT_DIR / "strict_msa_internal_gap_fraction_hist.png"
)
TOP_GAP_FRACTION = OUTPUT_DIR / "strict_msa_top_gap_fraction.png"
TOP_INTERNAL_GAP_FRACTION = (
    OUTPUT_DIR / "strict_msa_top_internal_gap_fraction.png"
)


@dataclass(frozen=True)
class AlignmentQcRecord:
    """
    Plotting fields from one alignment QC row.
    """

    protein_id: str
    species: str
    gap_fraction: float
    internal_gap_fraction: float

    @property
    def label(self) -> str:
        """
        Return the y-axis label for top-20 sequence plots.
        """

        return f"{self.species} {self.protein_id}"


def load_qc_records(path: Path) -> list[AlignmentQcRecord]:
    """
    Load alignment QC records from a TSV file.
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
                AlignmentQcRecord(
                    protein_id=row.get("protein_id", ""),
                    species=row.get("species", ""),
                    gap_fraction=float(row.get("gap_fraction", "0")),
                    internal_gap_fraction=float(
                        row.get("internal_gap_fraction", "0")
                    ),
                )
            )

    if not records:
        raise ValueError(f"Input TSV has no data records: {path}")

    return records


def save_histogram(
    values: list[float],
    output_path: Path,
    x_label: str,
    title: str,
) -> None:
    """
    Save a histogram for one QC metric.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(
        values,
        bins=20,
        color="#4C78A8",
        edgecolor="black",
    )
    ax.set_xlabel(x_label)
    ax.set_ylabel("number of sequences")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def save_top_bar_plot(
    records: list[AlignmentQcRecord],
    output_path: Path,
    metric_name: str,
    title: str,
) -> None:
    """
    Save a horizontal bar plot for the top 20 records by one metric.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)

    top_records = sorted(
        records,
        key=lambda record: (
            -getattr(record, metric_name),
            record.protein_id,
        ),
    )[:20]
    plot_records = list(reversed(top_records))
    labels = [record.label for record in plot_records]
    values = [getattr(record, metric_name) for record in plot_records]

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.barh(
        labels,
        values,
        color="#59A14F",
        edgecolor="black",
    )
    ax.set_xlabel(metric_name)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def run() -> None:
    """
    Create all alignment QC figures.
    """

    records = load_qc_records(INPUT_QC_TSV)
    gap_fractions = [record.gap_fraction for record in records]
    internal_gap_fractions = [
        record.internal_gap_fraction
        for record in records
    ]

    save_histogram(
        values=gap_fractions,
        output_path=GAP_FRACTION_HIST,
        x_label="gap_fraction",
        title="Strict MSA gap fraction distribution",
    )
    save_histogram(
        values=internal_gap_fractions,
        output_path=INTERNAL_GAP_FRACTION_HIST,
        x_label="internal_gap_fraction",
        title="Strict MSA internal gap fraction distribution",
    )
    save_top_bar_plot(
        records=records,
        output_path=TOP_GAP_FRACTION,
        metric_name="gap_fraction",
        title="Top 20 strict MSA sequences by gap fraction",
    )
    save_top_bar_plot(
        records=records,
        output_path=TOP_INTERNAL_GAP_FRACTION,
        metric_name="internal_gap_fraction",
        title="Top 20 strict MSA sequences by internal gap fraction",
    )

    print(f"records read: {len(records)}")
    print(f"saved: {GAP_FRACTION_HIST}")
    print(f"saved: {INTERNAL_GAP_FRACTION_HIST}")
    print(f"saved: {TOP_GAP_FRACTION}")
    print(f"saved: {TOP_INTERNAL_GAP_FRACTION}")


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

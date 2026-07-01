"""
Split PHF10 ortholog candidates into
core and extended sets based on
length-cluster assignment.

Pipeline position:

length_distribution
    ->
architecture_split
    ->
domain_scan
    ->
orthology_qc

This module performs no biological
interpretation.

It only partitions data.
"""

from pathlib import Path

import pandas as pd


INPUT_FASTA = "data/raw/orthology/phf10_orthologs.fasta"

INPUT_CLUSTERS = (
    "data/interim/orthology/length_clusters.tsv"
)

OUTPUT_TABLE = (
    "data/interim/orthology/architecture_split.tsv"
)

OUTPUT_CORE_FASTA = (
    "data/interim/orthology/core_orthologs.fasta"
)

OUTPUT_EXTENDED_FASTA = (
    "data/interim/orthology/extended_orthologs.fasta"
)


def load_cluster_map(cluster_file: str):

    df = pd.read_csv(cluster_file, sep="\t")

    cluster_map = {}

    for _, row in df.iterrows():

        protein_id = row["protein_id"]
        cluster = int(row["cluster"])

        if cluster in (0, 2):
            group = "CORE"
        else:
            group = "EXTENDED"

        cluster_map[protein_id] = group

    return df, cluster_map


def split_fasta(
    fasta_path: str,
    cluster_map: dict,
    core_out: str,
    extended_out: str,
):

    core_records = []
    extended_records = []

    with open(fasta_path, "r") as f:

        header = None
        sequence = []

        for line in f:

            if line.startswith(">"):

                if header is not None:

                    protein_id = (
                        header[1:].split("|")[0]
                    )

                    record = (
                        header
                        + "".join(sequence)
                    )

                    if (
                        cluster_map.get(protein_id)
                        == "CORE"
                    ):
                        core_records.append(record)
                    else:
                        extended_records.append(record)

                header = line
                sequence = []

            else:
                sequence.append(line)

        if header is not None:

            protein_id = (
                header[1:].split("|")[0]
            )

            record = (
                header
                + "".join(sequence)
            )

            if (
                cluster_map.get(protein_id)
                == "CORE"
            ):
                core_records.append(record)
            else:
                extended_records.append(record)

    Path(core_out).parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with open(core_out, "w") as f:
        f.writelines(core_records)

    with open(extended_out, "w") as f:
        f.writelines(extended_records)

    print(
        f"CORE records: {len(core_records)}"
    )

    print(
        f"EXTENDED records: "
        f"{len(extended_records)}"
    )


def save_table(
    df: pd.DataFrame,
    output_file: str,
):

    df = df.copy()

    df["group"] = df["cluster"].apply(
        lambda x:
        "CORE"
        if x in (0, 2)
        else "EXTENDED"
    )

    df.to_csv(
        output_file,
        sep="\t",
        index=False,
    )

    print(
        f"Saved: {output_file}"
    )


def main():

    df, cluster_map = load_cluster_map(
        INPUT_CLUSTERS
    )

    save_table(
        df,
        OUTPUT_TABLE,
    )

    split_fasta(
        INPUT_FASTA,
        cluster_map,
        OUTPUT_CORE_FASTA,
        OUTPUT_EXTENDED_FASTA,
    )


if __name__ == "__main__":
    main()

"""
Join GMM length clusters with taxonomy information.

Цель:
    - связать длину, кластер и вид
    - выявить таксономические паттерны в длине
    - отделить биологические сигналы от аннотационных артефактов

Без интерпретации — только структура данных.
"""

from pathlib import Path
import pandas as pd


def load(tsv_path: str) -> pd.DataFrame:
    return pd.read_csv(tsv_path, sep="\t")


def run(cluster_file: str, sorted_file: str, output_file: str):

    clusters = load(cluster_file)
    sorted_df = load(sorted_file)

    df = clusters.merge(
        sorted_df[["protein_id", "species_simple"]],
        on="protein_id",
        how="left"
    )

    # -----------------------------
    # 1. Cluster summary by species
    # -----------------------------

    print("\nSpecies distribution per cluster:\n")

    species_counts = df.groupby(["cluster", "species_simple"]).size()

    print(species_counts.sort_values(ascending=False).head(30))

    # -----------------------------
    # 2. Outlier detection per species
    # -----------------------------

    print("\nExtreme lengths per species:\n")

    outliers = df.groupby("species_simple")["length"].agg(["min", "max", "count"])
    outliers["range"] = outliers["max"] - outliers["min"]

    print(outliers.sort_values("range", ascending=False).head(20))

    # -----------------------------
    # 3. Singleton clusters
    # -----------------------------

    singletons = df[df["cluster"].map(df["cluster"].value_counts()) == 1]

    print("\nSingleton cluster entries:\n")
    print(singletons[["protein_id", "species_simple", "length", "cluster"]])

    # -----------------------------
    # Save
    # -----------------------------

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(output_file, sep="\t", index=False)

    print(f"\nSaved: {output_file}")
    print(f"Total rows: {len(df)}")


if __name__ == "__main__":

    run(
        cluster_file="data/interim/orthology/length_clusters.tsv",
        sorted_file="data/interim/orthology/sorted_by_length.tsv",
        output_file="data/interim/orthology/length_cluster_taxonomy.tsv",
    )

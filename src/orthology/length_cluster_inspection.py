"""
Inspection of GMM-based length clustering for PHF10 orthologs.

Цель:
    - проверить, что кластеры реально отражают структуру данных
    - исключить случайную сегментацию распределения
    - оценить связь кластеров с длинами и таксонами

Важно:
    Это диагностический модуль.
    Никаких решений о фильтрации не принимает.
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def load_data(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def summarize_clusters(df: pd.DataFrame):

    print("\nCluster summary:\n")

    summary = df.groupby("cluster").agg(
        count=("protein_id", "count"),
        mean_length=("length", "mean"),
        min_length=("length", "min"),
        max_length=("length", "max")
    )

    print(summary)


def plot_distribution(df: pd.DataFrame, output_path: str):

    plt.figure()

    for c in sorted(df["cluster"].unique()):
        subset = df[df["cluster"] == c]
        plt.scatter(subset["length"], [c] * len(subset), s=10, alpha=0.6)

    plt.xlabel("Protein length")
    plt.ylabel("Cluster")
    plt.title("PHF10 ortholog length clustering (GMM)")
    plt.tight_layout()

    plt.savefig(output_path, dpi=200)
    print(f"\nPlot saved to: {output_path}")


def run(input_tsv: str, output_plot: str):

    df = load_data(input_tsv)

    summarize_clusters(df)
    plot_distribution(df, output_plot)


if __name__ == "__main__":

    run(
        input_tsv="data/interim/orthology/length_clusters.tsv",
        output_plot="data/interim/orthology/length_clusters.png",
    )

"""
Length distribution analysis for PHF10 orthologs.

Цель:
    - изучить распределение длин белков
    - выявить скрытые кластеры (isoforms / truncations / full-length / fusions)
    - без ручных порогов
    - без биологических интерпретаций на этом этапе

Метод:
    Gaussian Mixture Model (GMM) + выбор числа компонент по BIC
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import csv
from typing import List, Dict, Tuple

from sklearn.mixture import GaussianMixture


# -----------------------------
# FASTA PARSING
# -----------------------------

def parse_fasta(fasta_path: str) -> List[Dict]:
    """
    FASTA format:
    >protein_id|species|gene_id
    SEQUENCE
    """

    records = []

    with open(fasta_path, "r", encoding="utf-8") as f:

        header = None
        seq_parts = []

        for line in f:
            line = line.strip()

            if line.startswith(">"):

                if header is not None:
                    records.append(_build(header, seq_parts))

                header = line[1:]
                seq_parts = []

            else:
                seq_parts.append(line)

        if header is not None:
            records.append(_build(header, seq_parts))

    return records


def _build(header: str, seq_parts: List[str]) -> Dict:
    parts = header.split("|")

    return {
        "protein_id": parts[0],
        "species": parts[1] if len(parts) > 1 else None,
        "gene_id": parts[2] if len(parts) > 2 else None,
        "sequence": "".join(seq_parts),
        "length": len("".join(seq_parts)),
    }


# -----------------------------
# CLUSTERING
# -----------------------------

def fit_gmm(lengths: np.ndarray, max_k: int = 5) -> Tuple[GaussianMixture, int]:

    """
    Выбор числа кластеров через BIC.
    """

    lengths = lengths.reshape(-1, 1)

    best_model = None
    best_bic = np.inf
    best_k = 1

    for k in range(1, max_k + 1):

        model = GaussianMixture(
            n_components=k,
            covariance_type="full",
            random_state=42
        )

        model.fit(lengths)
        bic = model.bic(lengths)

        if bic < best_bic:
            best_bic = bic
            best_model = model
            best_k = k

    return best_model, best_k


def assign_clusters(model: GaussianMixture, lengths: np.ndarray):

    X = lengths.reshape(-1, 1)

    clusters = model.predict(X)
    probs = model.predict_proba(X).max(axis=1)

    return clusters, probs


# -----------------------------
# MAIN PIPELINE
# -----------------------------

def run(
    fasta_path: str,
    output_tsv: str
) -> None:

    records = parse_fasta(fasta_path)

    lengths = np.array([r["length"] for r in records])

    model, best_k = fit_gmm(lengths, max_k=6)

    clusters, probs = assign_clusters(model, lengths)

    print(f"Sequences: {len(records)}")
    print(f"Selected clusters (GMM): {best_k}")
    print(f"BIC: {model.bic(lengths.reshape(-1,1)):.2f}")

    # cluster sizes
    unique, counts = np.unique(clusters, return_counts=True)
    for c, n in zip(unique, counts):
        print(f"Cluster {c}: {n} proteins")

    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=[
                "protein_id",
                "species",
                "gene_id",
                "length",
                "cluster",
                "cluster_probability"
            ],
            delimiter="\t"
        )

        writer.writeheader()

        for r, c, p in zip(records, clusters, probs):

            writer.writerow({
                "protein_id": r["protein_id"],
                "species": r["species"],
                "gene_id": r["gene_id"],
                "length": r["length"],
                "cluster": int(c),
                "cluster_probability": float(p)
            })


if __name__ == "__main__":

    run(
        fasta_path="data/raw/orthology/phf10_orthologs.fasta",
        output_tsv="data/interim/orthology/length_clusters.tsv",
    )

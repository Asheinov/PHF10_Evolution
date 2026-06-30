# src/idr/run_iupred.py

import json
import pandas as pd
from pathlib import Path

from src.utils.fasta import read_fasta


def load_iupred(json_file: str):

    with open(json_file) as f:
        data = json.load(f)

    return data


def build_dataframe(json_file: str, fasta_file: str):

    data = load_iupred(json_file)

    sequence = read_fasta(fasta_file)

    iupred_scores = data["iupred2"]

    if len(sequence) != len(iupred_scores):
        raise ValueError(
            f"Length mismatch: "
            f"sequence={len(sequence)}, "
            f"iupred={len(iupred_scores)}"
        )

    df = pd.DataFrame({
        "position": range(1, len(sequence) + 1),
        "aa": list(sequence),
        "iupred": iupred_scores
    })

    return df


def main():

    accession = "Q8WUB8"

    json_file = f"data/raw/iupred3_{accession}.json"
    fasta_file = f"data/raw/{accession}.fasta"

    df = build_dataframe(
        json_file,
        fasta_file
    )

    print(df.head())

    outdir = Path("data/02_interim/iupred")
    outdir.mkdir(parents=True, exist_ok=True)

    outfile = outdir / f"{accession}_iupred.csv"

    df.to_csv(outfile, index=False)

    print(f"Saved: {outfile}")


if __name__ == "__main__":
    main()

import requests
import pandas as pd
from pathlib import Path
from src.utils.fasta import read_fasta



def get_aiupred(accession: str):

    url = "https://aiupred.elte.hu/rest_api"

    r = requests.get(
        url,
        params={
            "accession": accession,
            "analysis_type": "binding"
        },
        timeout=60
    )

    r.raise_for_status()

    return r.json()


def build_dataframe(accession: str, fasta_file: str):


    data = get_aiupred(accession)

    disorder = data["AIUPred"]
    binding = data["AIUPred-binding"]

    
    sequence = read_fasta(fasta_file)
    if len(sequence) != len(disorder):
        raise ValueError(
        f"Length mismatch: sequence={len(sequence)}, aiupred={len(disorder)}"
    )

    df = pd.DataFrame({
        "position": range(1, len(disorder) + 1),
        "aa": list(sequence),
        "disorder": disorder,
        "binding": binding
    })

    return df


def main():

    accession = "Q8WUB8"
    fasta_file = "data/01_raw/Q8WUB8.fasta"

    df = build_dataframe(
        accession,
        fasta_file
    )

    print(df.head())

    outdir = Path("data/02_interim/aiupred")
    outdir.mkdir(parents=True, exist_ok=True)

    outfile = outdir / f"{accession}_aiupred.csv"

    df.to_csv(outfile, index=False)

    print(f"Saved: {outfile}")

    


if __name__ == "__main__":
    main()

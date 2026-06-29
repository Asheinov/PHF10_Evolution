from pathlib import Path
import json

import numpy as np
import pandas as pd
import metapredict as meta


# ==========================================================
# CONFIG
# ==========================================================

DATA_DIR = Path("data/01_raw")

FASTA_FILE = DATA_DIR / "Q8WUB8.fasta"
PDB_FILE = DATA_DIR / "AF-Q8WUB8.pdb"
IUPRED_FILE = DATA_DIR / "iupred3_Q8WUB8.json"
FLDPNN_FILE = DATA_DIR / "fldpnn_Q8WUB8.csv"

DISORDER_THRESHOLD = 0.5
PLDDT_THRESHOLD = 50.0

MIN_IDR_LENGTH = 20


# ==========================================================
# LOADERS
# ==========================================================

def read_fasta(path):
    seq = []

    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())

    return "".join(seq)


def extract_plddt(path):

    plddt = []

    with open(path) as f:
        for line in f:

            if not line.startswith("ATOM"):
                continue

            atom_name = line[12:16].strip()

            if atom_name != "CA":
                continue

            plddt.append(float(line[60:66]))

    if len(plddt) == 0:
        raise ValueError(f"No pLDDT values found in {path}")

    return np.array(plddt)


def load_iupred(path):

    with open(path) as f:
        data = json.load(f)

    scores = np.array(data["iupred2"])

    return scores


def load_fldpnn(path):

    with open(path) as f:

        for line in f:

            clean = line.replace(",", "").strip()

            if (
                len(clean) > 0
                and all(c in "01" for c in clean)
            ):

                return np.array(
                    [int(x) for x in clean],
                    dtype=float
                )

    raise ValueError(
        f"Cannot parse flDPnn output: {path}"
    )


# ==========================================================
# HELPERS
# ==========================================================

def validate_lengths(**arrays):

    lengths = {
        k: len(v)
        for k, v in arrays.items()
    }

    unique_lengths = set(lengths.values())

    if len(unique_lengths) != 1:
        raise ValueError(
            f"Length mismatch:\n{lengths}"
        )


def find_regions(mask, min_len=20):

    regions = []

    start = None

    for i, val in enumerate(mask):

        if val == 1 and start is None:
            start = i

        elif val == 0 and start is not None:

            end = i

            if end - start >= min_len:
                regions.append(
                    (start + 1, end)
                )

            start = None

    if start is not None:

        end = len(mask)

        if end - start >= min_len:
            regions.append(
                (start + 1, end)
            )

    return regions


# ==========================================================
# MAIN
# ==========================================================

def main():

    print("Loading PHF10...")

    sequence = read_fasta(FASTA_FILE)

    meta_scores = np.array(
        meta.predict_disorder(sequence)
    )

    iupred_scores = load_iupred(IUPRED_FILE)

    fldpnn_scores = load_fldpnn(FLDPNN_FILE)

    plddt_scores = extract_plddt(PDB_FILE)

    validate_lengths(
        sequence=sequence,
        metapredict=meta_scores,
        iupred=iupred_scores,
        fldpnn=fldpnn_scores,
        plddt=plddt_scores,
    )

    # ---------------------------------------
    # Binary predictions
    # ---------------------------------------

    meta_bin = (
        meta_scores > DISORDER_THRESHOLD
    ).astype(int)

    iupred_bin = (
        iupred_scores > DISORDER_THRESHOLD
    ).astype(int)

    fldpnn_bin = (
        fldpnn_scores > 0.5
    ).astype(int)

    plddt_bin = (
        plddt_scores < PLDDT_THRESHOLD
    ).astype(int)

    # ---------------------------------------
    # Consensus
    # ---------------------------------------

    consensus_votes = (
        meta_bin
        + iupred_bin
        + fldpnn_bin
        + plddt_bin
    )

    consensus_mask = (
        consensus_votes >= 3
    ).astype(int)

    confidence = (
        consensus_votes / 4.0
    )

    # ---------------------------------------
    # Results table
    # ---------------------------------------

    df = pd.DataFrame({

        "position":
            np.arange(
                1,
                len(sequence) + 1
            ),

        "aa":
            list(sequence),

        "metapredict":
            meta_scores,

        "iupred":
            iupred_scores,

        "fldpnn":
            fldpnn_scores,

        "plddt":
            plddt_scores,

        "votes":
            consensus_votes,

        "confidence":
            confidence,

        "consensus_idr":
            consensus_mask
    })

    # ---------------------------------------
    # Regions
    # ---------------------------------------

    idr_regions = find_regions(
        consensus_mask,
        min_len=MIN_IDR_LENGTH
    )

    print("\nConsensus IDRs:")

    for start, end in idr_regions:

        print(
            f"IDR {start}-{end} "
            f"(length={end-start+1})"
        )

    print(
        f"\nTotal IDR residues: "
        f"{consensus_mask.sum()}"
    )

    # ---------------------------------------
    # Save
    # ---------------------------------------

    df.to_csv(
        "phf10_idr_consensus.csv",
        index=False
    )

    print(
        "\nSaved: phf10_idr_consensus.csv"
    )


if __name__ == "__main__":
    main()
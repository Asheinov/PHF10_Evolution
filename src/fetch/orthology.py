"""
Orthology fetch pipeline for PHF10 Evolution project.

Pipeline:
1. Map UniProt ID → Ensembl Gene ID
2. Fetch orthologs from:
   - Ensembl Compara
   - OMA
   - OrthoDB
3. Normalize all results into unified TSV schema
4. Save per-source files into data/raw/orthology/

This module performs ONLY:
- data retrieval (fetch)
- normalization
- saving raw standardized outputs

No aggregation or confidence scoring is performed here.
"""

from dataclasses import asdict
from typing import List, Dict, Any
import os
import csv
import requests



from src.orthology.models import OrthologHit


# ------------------------------------------------------------
# 1. ID MAPPING
# ------------------------------------------------------------

def map_uniprot_to_ensembl_gene(uniprot_id: str) -> str:
    """
    UniProt → Ensembl Gene ID via UniProt cross-reference
    """

    import requests

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"

    r = requests.get(url)

    if not r.ok:
        raise RuntimeError(f"Failed UniProt fetch: {uniprot_id}")

    data = r.json()

    for ref in data.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "Ensembl":

            for prop in ref.get("properties", []):
                if prop.get("key") == "GeneId":
                    return prop.get("value")

    raise RuntimeError(f"No Ensembl Gene ID found for {uniprot_id}")


# ------------------------------------------------------------
# 2. FETCH LAYERS (stubs for now)
# ------------------------------------------------------------

def fetch_ensembl_orthologs(ensembl_gene_id: str):

    import requests

    headers = {"Content-Type": "application/json"}

    clean_id = ensembl_gene_id.split(".")[0]

    species = "homo_sapiens"

    url = (
        f"https://rest.ensembl.org/homology/id/"
        f"{species}/{clean_id}?type=orthologues"
    )

    # print("REQUEST URL:", url)

    r = requests.get(url, headers=headers)

    # print("STATUS:", r.status_code)
    # print("RESPONSE SAMPLE:", r.text[:300])

    if not r.ok:
        raise RuntimeError(f"Ensembl fetch failed for {clean_id}")

    data = r.json()

    results = []

    for entry in data.get("data", []):

        query_id = entry.get("id")

        for homology in entry.get("homologies", []):

            # иногда Ensembl кладёт данные в target
            target = homology.get("target") or homology.get("source")

            if not target:
                continue

            target_gene = target.get("id") or target.get("gene_id")

            # 1. строгая проверка валидности
            if target_gene is None:
                continue

            # 2. фильтр self-hit (человек → человек тот же ген)
            if target_gene == query_id:
                continue

            results.append({
                "target_gene": target_gene,
                "target_species": target.get("species"),
                "target_protein_id": target.get("protein_id"),

                "orthology_type": homology.get("type"),
                "raw_score": homology.get("confidence", 0.0),
        })

    return results
    


def fetch_oma_orthologs(ensembl_gene_id: str) -> List[Dict[str, Any]]:
    """
    Fetch orthologs from OMA.
    """
    raise NotImplementedError("Implement OMA fetch")


def fetch_orthodb_orthologs(ensembl_gene_id: str) -> List[Dict[str, Any]]:
    """
    Fetch orthologs from OrthoDB.
    """
    raise NotImplementedError("Implement OrthoDB fetch")


# ------------------------------------------------------------
# 3. NORMALIZATION
# ------------------------------------------------------------

def normalize_to_standard_schema(
    raw_records: List[Dict[str, Any]],
    source: str,
    query_gene: str,
    query_species: str,
    query_protein_id: str
) -> List[OrthologHit]:
    """
    Convert raw API/dump records into OrthologHit objects.
    """

    normalized = []

    for r in raw_records:

        hit = OrthologHit(
            source=source,

            query_gene=query_gene,
            query_species=query_species,
            query_protein_id=query_protein_id,

            target_gene=r.get("target_gene"),
            target_species=r.get("target_species"),
            target_protein_id=r.get("target_protein_id"),

            orthology_type=r.get("orthology_type"),
            raw_score=r.get("raw_score"),
        )

        normalized.append(hit)

    return normalized


# ------------------------------------------------------------
# 4. TSV EXPORT
# ------------------------------------------------------------

def save_tsv(path: str, hits: List[OrthologHit]) -> None:
    """
    Save list of OrthologHit into TSV file.
    """

    os.makedirs(os.path.dirname(path), exist_ok=True)

    fieldnames = [
        "query_gene",
        "query_species",
        "query_protein_id",
        "target_gene",
        "target_species",
        "target_protein_id",
        "orthology_type",
        "raw_score",
        "source",
    ]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for h in hits:
            row = asdict(h)
            writer.writerow({k: row.get(k) for k in fieldnames})


# ------------------------------------------------------------
# 5. MAIN PIPELINE
# ------------------------------------------------------------

def run_orthology_fetch(
    input_uniprot: str,
    config: Dict[str, Any]
) -> None:
    """
    Full pipeline execution.
    """

    # 1. Map UniProt → Ensembl Gene ID
    ensembl_gene_id = map_uniprot_to_ensembl_gene(input_uniprot)

    query_species = config.get("query_species", "Homo sapiens")
    query_gene = config.get("query_gene", "UNKNOWN")
    query_protein_id = input_uniprot

    out_dir = config.get("output_dir", "data/raw/orthology")

    # --------------------------------------------------------
    # 2. Fetch raw data
    # --------------------------------------------------------

    ensembl_raw = fetch_ensembl_orthologs(ensembl_gene_id)
    #oma_raw = fetch_oma_orthologs(ensembl_gene_id)
    #orthodb_raw = fetch_orthodb_orthologs(ensembl_gene_id)
    oma_raw = []
    orthodb_raw = []
    # --------------------------------------------------------
    # 3. Normalize
    # --------------------------------------------------------

    ensembl_hits = normalize_to_standard_schema(
        ensembl_raw,
        source="Ensembl",
        query_gene=query_gene,
        query_species=query_species,
        query_protein_id=query_protein_id
    )

    oma_hits = normalize_to_standard_schema(
        oma_raw,
        source="OMA",
        query_gene=query_gene,
        query_species=query_species,
        query_protein_id=query_protein_id
    )

    orthodb_hits = normalize_to_standard_schema(
        orthodb_raw,
        source="OrthoDB",
        query_gene=query_gene,
        query_species=query_species,
        query_protein_id=query_protein_id
    )

    # --------------------------------------------------------
    # 4. Save TSVs
    # --------------------------------------------------------

    save_tsv(os.path.join(out_dir, "ensembl.tsv"), ensembl_hits)
    save_tsv(os.path.join(out_dir, "oma.tsv"), oma_hits)
    save_tsv(os.path.join(out_dir, "orthodb.tsv"), orthodb_hits)

if __name__ == "__main__":

    config = {
        "query_gene": "PHF10",
        "query_species": "Homo sapiens",
        "output_dir": "data/raw/orthology"
    }

    run_orthology_fetch("Q8WUB8", config)
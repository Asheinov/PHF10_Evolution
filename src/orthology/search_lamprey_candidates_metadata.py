"""
Search Ensembl metadata for lamprey PHF10 candidate ortholog records.

Purpose:
    Query the Ensembl REST homology endpoint with the human PHF10 gene
    symbol and target lamprey species from the missing-lamprey fetch plan.
    This is a metadata-only candidate discovery step for later manual
    orthology, domain, and readiness QC.

Input files:
    data/interim/orthology/taxonomy/missing_lamprey_fetch_plan.tsv

Output files:
    data/interim/orthology/taxonomy/lamprey_candidate_metadata.tsv
    data/raw/orthology/lamprey_candidate_metadata_raw.json

Pipeline position:
    missing_lamprey_fetch_plan
        ->
    search_lamprey_candidates_metadata

Biological limitation:
    This module does not fetch protein FASTA sequences into the PHF10
    working set, classify final orthology, or use phosphomodule presence,
    Ser-rich sequence, N-terminal extension, domains, IDR, or trait-level
    information as candidate criteria. It only discovers candidate records
    using source metadata and homology annotations.
"""

from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any


FETCH_PLAN_TSV = Path(
    "data/interim/orthology/taxonomy/missing_lamprey_fetch_plan.tsv"
)
OUTPUT_METADATA_TSV = Path(
    "data/interim/orthology/taxonomy/lamprey_candidate_metadata.tsv"
)
OUTPUT_RAW_JSON = Path(
    "data/raw/orthology/lamprey_candidate_metadata_raw.json"
)

ENSEMBL_HOMOLOGY_BASE_URL = "https://rest.ensembl.org/homology/symbol"
SOURCE = "Ensembl"
QUERY_SPECIES = "homo_sapiens"
QUERY_SYMBOL = "PHF10"
HOMOLOGY_TYPE = "orthologues"
USER_AGENT = "PHF10_Evolution/1.0 metadata-only lamprey candidate search"
REQUEST_TIMEOUT_SECONDS = 30
POLITE_SLEEP_SECONDS = 1.0
SEQUENCE_FETCH_STATUS = "NOT_FETCHED"
METADATA_ONLY_STATUS = "CANDIDATE_METADATA_ONLY"
NO_CANDIDATE_STATUS = "NO_METADATA_CANDIDATE"

FETCH_PLAN_REQUIRED_FIELDS = {
    "species",
    "taxon_group",
    "priority",
    "preferred_sources",
    "search_stage",
    "fetch_status",
}
OUTPUT_FIELDNAMES = [
    "source",
    "query_species",
    "query_symbol",
    "target_species",
    "candidate_gene_id",
    "candidate_protein_id",
    "candidate_transcript_id",
    "candidate_symbol",
    "homology_type",
    "orthology_confidence",
    "perc_id",
    "perc_pos",
    "protein_sequence_available",
    "sequence_fetch_status",
    "raw_record_index",
    "manual_review_status",
    "manual_note",
]


class EnsemblRequestError(RuntimeError):
    """
    Raised when Ensembl cannot be queried or returns a non-success status.
    """


class MalformedJsonError(RuntimeError):
    """
    Raised when Ensembl returns a response that is not parseable JSON.
    """


@dataclass(frozen=True)
class FetchPlanRow:
    """
    One lamprey target row from the metadata fetch plan.
    """

    species: str
    taxon_group: str
    priority: str
    preferred_sources: str
    search_stage: str
    fetch_status: str


@dataclass(frozen=True)
class CandidateMetadataRow:
    """
    One metadata-only candidate row for manual review.
    """

    source: str
    query_species: str
    query_symbol: str
    target_species: str
    candidate_gene_id: str
    candidate_protein_id: str
    candidate_transcript_id: str
    candidate_symbol: str
    homology_type: str
    orthology_confidence: str
    perc_id: str
    perc_pos: str
    protein_sequence_available: str
    sequence_fetch_status: str
    raw_record_index: str
    manual_review_status: str
    manual_note: str

    def to_dict(self) -> dict[str, str]:
        """
        Convert the candidate row to a DictWriter-compatible dictionary.
        """

        return {
            "source": self.source,
            "query_species": self.query_species,
            "query_symbol": self.query_symbol,
            "target_species": self.target_species,
            "candidate_gene_id": self.candidate_gene_id,
            "candidate_protein_id": self.candidate_protein_id,
            "candidate_transcript_id": self.candidate_transcript_id,
            "candidate_symbol": self.candidate_symbol,
            "homology_type": self.homology_type,
            "orthology_confidence": self.orthology_confidence,
            "perc_id": self.perc_id,
            "perc_pos": self.perc_pos,
            "protein_sequence_available": self.protein_sequence_available,
            "sequence_fetch_status": self.sequence_fetch_status,
            "raw_record_index": self.raw_record_index,
            "manual_review_status": self.manual_review_status,
            "manual_note": self.manual_note,
        }


def require_input_file(path: Path) -> None:
    """
    Fail clearly when an expected input file is missing.
    """

    if not path.exists():
        raise FileNotFoundError(f"Required input file is missing: {path}")


def validate_fields(
    path: Path,
    fieldnames: list[str] | None,
    required_fields: set[str],
) -> None:
    """
    Check that a TSV header contains all fields required by this module.
    """

    if fieldnames is None:
        raise ValueError(f"Input TSV has no header: {path}")

    missing_fields = sorted(required_fields.difference(fieldnames))
    if missing_fields:
        fields = ", ".join(missing_fields)
        raise ValueError(f"Input TSV is missing required fields: {path}: {fields}")


def load_fetch_plan(path: Path) -> list[FetchPlanRow]:
    """
    Read target lamprey species from the missing-lamprey fetch plan.
    """

    require_input_file(path)
    rows = []

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_fields(path, reader.fieldnames, FETCH_PLAN_REQUIRED_FIELDS)

        for row in reader:
            species = row["species"].strip()
            if not species:
                raise ValueError(f"Fetch plan row has empty species: {path}")

            rows.append(
                FetchPlanRow(
                    species=species,
                    taxon_group=row["taxon_group"].strip(),
                    priority=row["priority"].strip(),
                    preferred_sources=row["preferred_sources"].strip(),
                    search_stage=row["search_stage"].strip(),
                    fetch_status=row["fetch_status"].strip(),
                )
            )

    if not rows:
        raise ValueError(f"Input TSV has no lamprey target rows: {path}")

    return rows


def build_ensembl_homology_url(target_species: str) -> str:
    """
    Build an Ensembl REST homology URL for one target species.
    """

    species = urllib.parse.quote(QUERY_SPECIES, safe="")
    symbol = urllib.parse.quote(QUERY_SYMBOL, safe="")
    params = urllib.parse.urlencode(
        {
            "type": HOMOLOGY_TYPE,
            "target_species": target_species,
            "content-type": "application/json",
        }
    )
    return f"{ENSEMBL_HOMOLOGY_BASE_URL}/{species}/{symbol}?{params}"


def get_json(url: str, timeout: int = REQUEST_TIMEOUT_SECONDS) -> dict[str, Any]:
    """
    Fetch a JSON document with timeout, user-agent, and readable errors.
    """

    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        },
    )

    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            status = getattr(response, "status", "")
            body = response.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        context = exc.read().decode("utf-8", errors="replace")[:500]
        raise EnsemblRequestError(
            f"Ensembl request failed: status={exc.code} url={url} "
            f"context={context}"
        ) from exc
    except urllib.error.URLError as exc:
        raise EnsemblRequestError(
            f"Ensembl request failed: url={url} context={exc.reason}"
        ) from exc
    except TimeoutError as exc:
        raise EnsemblRequestError(
            f"Ensembl request timed out: url={url} timeout={timeout}"
        ) from exc

    try:
        parsed = json.loads(body)
    except json.JSONDecodeError as exc:
        raise MalformedJsonError(
            f"Ensembl response was not valid JSON: status={status} url={url}"
        ) from exc

    if not isinstance(parsed, dict):
        raise MalformedJsonError(
            f"Ensembl response JSON was not an object: status={status} url={url}"
        )

    return parsed


def first_string_value(record: dict[str, Any], keys: list[str]) -> str:
    """
    Return the first non-empty scalar value for the given keys.
    """

    for key in keys:
        value = record.get(key)
        if value is not None and not isinstance(value, (dict, list)):
            text = str(value).strip()
            if text:
                return text
    return ""


def has_protein_sequence_metadata(target: dict[str, Any]) -> bool:
    """
    Check whether the homology target explicitly includes protein metadata.
    """

    protein_keys = [
        "protein_id",
        "protein_stable_id",
        "translation_id",
        "translation_stable_id",
        "peptide_id",
        "protein_sequence",
        "sequence",
        "seq",
        "align_seq",
    ]
    return any(first_string_value(target, [key]) for key in protein_keys)


def iter_homology_targets(
    payload: dict[str, Any],
    target_species: str,
) -> list[tuple[int, dict[str, Any], dict[str, Any]]]:
    """
    Extract target homology records from an Ensembl homology response.
    """

    data = payload.get("data", [])
    if not isinstance(data, list):
        return []

    matches = []
    raw_record_index = 0

    for data_item in data:
        if not isinstance(data_item, dict):
            continue
        homologies = data_item.get("homologies", [])
        if not isinstance(homologies, list):
            continue

        for homology in homologies:
            if not isinstance(homology, dict):
                continue
            target = homology.get("target", {})
            if not isinstance(target, dict):
                continue
            if target.get("species", "") != target_species:
                continue

            raw_record_index += 1
            matches.append((raw_record_index, homology, target))

    return matches


def candidate_from_homology(
    target_species: str,
    raw_record_index: int,
    homology: dict[str, Any],
    target: dict[str, Any],
) -> CandidateMetadataRow:
    """
    Convert one Ensembl homology target into a metadata-only candidate row.
    """

    return CandidateMetadataRow(
        source=SOURCE,
        query_species=QUERY_SPECIES,
        query_symbol=QUERY_SYMBOL,
        target_species=target_species,
        candidate_gene_id=first_string_value(
            target,
            ["id", "gene_id", "gene_stable_id", "stable_id"],
        ),
        candidate_protein_id=first_string_value(
            target,
            [
                "protein_id",
                "protein_stable_id",
                "translation_id",
                "translation_stable_id",
                "peptide_id",
            ],
        ),
        candidate_transcript_id=first_string_value(
            target,
            ["transcript_id", "transcript_stable_id", "canonical_transcript"],
        ),
        candidate_symbol=first_string_value(
            target,
            ["display_id", "external_name", "gene_symbol", "symbol"],
        ),
        homology_type=first_string_value(homology, ["type"]),
        orthology_confidence=first_string_value(
            homology,
            ["confidence", "is_high_confidence", "orthology_confidence"],
        ),
        perc_id=first_string_value(target, ["perc_id", "percentage_identity"]),
        perc_pos=first_string_value(
            target,
            ["perc_pos", "percentage_positive"],
        ),
        protein_sequence_available=str(has_protein_sequence_metadata(target)),
        sequence_fetch_status=SEQUENCE_FETCH_STATUS,
        raw_record_index=str(raw_record_index),
        manual_review_status=METADATA_ONLY_STATUS,
        manual_note="",
    )


def no_candidate_row(
    target_species: str,
    manual_note: str = "",
) -> CandidateMetadataRow:
    """
    Build a placeholder row when no Ensembl candidate metadata is returned.
    """

    return CandidateMetadataRow(
        source=SOURCE,
        query_species=QUERY_SPECIES,
        query_symbol=QUERY_SYMBOL,
        target_species=target_species,
        candidate_gene_id="",
        candidate_protein_id="",
        candidate_transcript_id="",
        candidate_symbol="",
        homology_type="",
        orthology_confidence="",
        perc_id="",
        perc_pos="",
        protein_sequence_available="False",
        sequence_fetch_status=SEQUENCE_FETCH_STATUS,
        raw_record_index="",
        manual_review_status=NO_CANDIDATE_STATUS,
        manual_note=manual_note,
    )


def query_target_species(
    target_species: str,
) -> tuple[list[CandidateMetadataRow], dict[str, Any]]:
    """
    Query Ensembl for one target species and return candidate rows plus raw data.
    """

    url = build_ensembl_homology_url(target_species)

    try:
        payload = get_json(url)
    except EnsemblRequestError as exc:
        raw_entry = {
            "source": SOURCE,
            "query_species": QUERY_SPECIES,
            "query_symbol": QUERY_SYMBOL,
            "target_species": target_species,
            "url": url,
            "ok": False,
            "error": str(exc),
        }
        return (
            [
                no_candidate_row(
                    target_species,
                    "No Ensembl metadata candidate returned; inspect raw JSON.",
                )
            ],
            raw_entry,
        )

    matches = iter_homology_targets(payload, target_species)
    raw_entry = {
        "source": SOURCE,
        "query_species": QUERY_SPECIES,
        "query_symbol": QUERY_SYMBOL,
        "target_species": target_species,
        "url": url,
        "ok": True,
        "response": payload,
    }

    if not matches:
        return ([no_candidate_row(target_species)], raw_entry)

    rows = [
        candidate_from_homology(
            target_species=target_species,
            raw_record_index=raw_record_index,
            homology=homology,
            target=target,
        )
        for raw_record_index, homology, target in matches
    ]
    return rows, raw_entry


def search_lamprey_candidates(
    plan_rows: list[FetchPlanRow],
) -> tuple[list[CandidateMetadataRow], dict[str, Any]]:
    """
    Search Ensembl metadata for all target species in the fetch plan.
    """

    candidate_rows = []
    raw_entries = []

    for index, plan_row in enumerate(plan_rows):
        if index > 0:
            time.sleep(POLITE_SLEEP_SECONDS)

        rows, raw_entry = query_target_species(plan_row.species)
        candidate_rows.extend(rows)
        raw_entries.append(raw_entry)

    raw_bundle = {
        "query": {
            "source": SOURCE,
            "query_species": QUERY_SPECIES,
            "query_symbol": QUERY_SYMBOL,
            "homology_type": HOMOLOGY_TYPE,
        },
        "queries": raw_entries,
    }
    return candidate_rows, raw_bundle


def write_candidate_metadata(
    records: list[CandidateMetadataRow],
    output_path: Path,
) -> None:
    """
    Write candidate metadata rows as a tab-separated table.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_FIELDNAMES,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for record in records:
            writer.writerow(record.to_dict())


def write_raw_json(raw_data: dict[str, Any], output_path: Path) -> None:
    """
    Write complete raw Ensembl JSON responses and query metadata.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(raw_data, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> None:
    """
    Run the metadata-only lamprey candidate search.
    """

    plan_rows = load_fetch_plan(FETCH_PLAN_TSV)
    candidate_rows, raw_data = search_lamprey_candidates(plan_rows)
    write_candidate_metadata(candidate_rows, OUTPUT_METADATA_TSV)
    write_raw_json(raw_data, OUTPUT_RAW_JSON)

    print(f"lamprey target species read: {len(plan_rows)}")
    print(f"candidate metadata rows written: {len(candidate_rows)}")
    print(f"saved: {OUTPUT_METADATA_TSV}")
    print(f"saved: {OUTPUT_RAW_JSON}")


if __name__ == "__main__":
    main()

# PHF10 Evolution Project — instructions for Codex

You are working on a reproducible bioinformatics project for PHF10 evolution.

## Core rules

1. Do not invent existing code.

Before changing anything, inspect the repository files. If a file, function, class, or field is not present in the repository, treat it as absent.

2. Always distinguish:

* DONE — already implemented in code
* TODO — agreed but not implemented
* IDEA — conceptual discussion only

Do not describe already implemented code in conditional language.

3. Use the existing project structure.

Current preferred layout:

src/

* fetch/
* orthology/
* idr/
* msa_tree/
* motifs/
* phosphosites/
* null_model/
* utils/

Do not create folders or Python packages starting with digits, such as `00_fetch` or `01_orthology`.

4. One file = one responsibility.

Do not create modules that simultaneously fetch data, analyze data, filter data, and export final results.

5. First inspect, then patch.

Before making changes, report:

* files inspected
* current implementation found
* exact files to modify
* risks
* expected outputs

6. Keep changes minimal.

Do not refactor unrelated modules. Do not rename files unless the task explicitly requires it.

7. Do not delete working code unless explicitly instructed.

If code is problematic, prefer:

* isolate it
* mark it as TODO
* replace only the unsafe part

8. Do not make Git commits unless explicitly asked.

You may suggest commit commands, but do not commit automatically.

## Scientific rules

1. Orthology, phylogeny, traits, and evolution are separate layers.

Do not use traits such as N-terminal extension, Ser-rich regions, IDR presence, or domain length to define orthology.

2. Treat every signal as a possible artifact.

Especially:

* gain/loss events
* missing N-terminal extensions
* missing IDRs
* long proteins
* truncated proteins
* extreme sequence lengths
* singleton clusters

3. Do not use Ser-rich regions to identify orthologs.

This would create circular reasoning.

4. Do not use N-terminal extension to identify orthologs.

It is a downstream evolutionary trait, not an orthology criterion.

5. Domain architecture is QC, not orthology definition.

SAY/PHD domain detection may be used to flag unusual architecture, but not to decide orthology by itself.

## Coding rules

1. Prefer explicit inputs and outputs.

Avoid hidden state and hardcoded biological assumptions.

2. Avoid hardcoded paths where possible.

If a path is used in a script entry point, keep it visible and easy to override later.

3. Use type hints and docstrings.

Every module should explain:

* purpose
* input files
* output files
* pipeline position

4. Network code must be safe.

For REST API calls:

* use timeout
* handle requests exceptions
* report status code and URL on failure
* avoid infinite blocking
* do not assume the server is always available

5. Fetch modules should save raw or standardized local files.

Downstream analysis should use local files, not repeated API calls.

6. Scripts should be runnable from project root.

Preferred command style:

python -m src.package.module

## Current project state

DONE:

* `src/orthology/models.py`
  Contains orthology dataclasses.

* `src/fetch/orthology.py`
  Initial Ensembl orthology fetch pipeline.
  Fetches PHF10 ortholog records from Ensembl Compara and exports TSV.

* `src/orthology/export_fasta.py`
  Exports protein sequences for Ensembl orthologs to FASTA.

* `src/orthology/qc_fasta.py`
  Computes basic sequence length QC.

* `src/orthology/length_distribution.py`
  Clusters protein lengths with GMM.

* `src/orthology/length_cluster_inspection.py`
  Summarizes length clusters.

* `src/orthology/sort_by_length.py`
  Sorts ortholog proteins by length.

* `src/orthology/length_cluster_taxonomy_join.py`
  Joins length clusters with species information.

* `src/orthology/architecture_split.py`
  Splits sequences into core and extended groups based on length clusters.

TODO:

* Stabilize `src/fetch/orthology.py`.
* Add safe HTTP wrapper.
* Remove hardcoded `homo_sapiens`.
* Add source flags so OMA and OrthoDB stubs are not executed unless requested.
* Preserve `orthology_type` but do not filter by it during fetch.
* Add TSV schema validation layer.
* Add domain resource setup for HMMER/Pfam.
* Add domain scan module only after external dependencies are prepared.

## Output format for Codex responses

For every task, respond with:

Status:
DONE / TODO / PARTIAL

Files changed:

* path/to/file.py

What changed:

* concise list

Validation:

* commands run
* results
* failures, if any

Not changed:

* explicitly mention important things left untouched

Next recommended step:

* one concrete next task

# Macroalgae genome–environment analysis code

This repository contains the analysis and figure-support code for the manuscript *Earth-Observation and Hybrid Spatiotemporal Embeddings Reveal Genome–Environment Associations in Macroalgae*.

All genome–environment analyses reported in the manuscript use the authenticated 126-genome cohort: 70 Rhodophyta, 43 Ochrophyta, and 13 Chlorophyta. The code does not substitute simulated, synthetic, or reconstructed summary values for experimental or observational data.

## Repository scope

The current revision workflows are organized as follows:

- `ISCIENCE_REVISION_20260711/integrity/`: cohort reconciliation, raw HMMER/Pfam count reconstruction, HMM header audit, and UAE metadata audit.
- `ISCIENCE_REVISION_20260711/aef/`: year-explicit exact-ID AEF extraction plus the corrected 126-genome, 10,707-Pfam by 64-axis screen, global multiplicity correction, within-phylum specifications, phylum-centered analysis, and Figure 5 generator.
- `ISCIENCE_REVISION_20260711/gee_validation/`: exact-genome-ID Earth Engine extraction and the 87,123-test raw-count correlation workflow.
- `ISCIENCE_REVISION_20260711/analysis_stats/`: the seven post hoc robustness checks, including site aggregation, alternate denominators, assembly-quality adjustment, completeness thresholds, topology-aware adjustment, structured permutations, and bootstrap intervals.
- `ISCIENCE_REVISION_20260711/pf00092_midas/`: PF00092.35 reannotation and canonical MIDAS-position analysis using real peptide sequences.
- `ISCIENCE_REVISION_20260711/figure3_126/`: the corrected 42,828-test recorded-metadata screen and Figure 3 generator for the exact 126-genome cohort.
- `ISCIENCE_REVISION_20260711/figure4_126/`: the within-phylum Figure 4 analysis on the 126-genome cohort.
- `ISCIENCE_REVISION_20260711/spatial/` and `ISCIENCE_REVISION_20260711/figures/`: spatial sensitivity and Figure 1 layout support.
- `ISCIENCE_REVISION_20260711/supplement/`: builders for corrected Tables S1–S3 and Figure S3 support.
- `triangulation/scripts/`: retained rbcL/FastTree workflow.

## Reproducibility boundaries

Large sequence files, raw HMMER output, derived matrices, and submission workbooks are not duplicated in this code repository. Their expected locations, public accessions, and service requirements are listed in [DATA_INPUTS.md](DATA_INPUTS.md).

Several final figures include documented iTOL or Adobe Illustrator assembly. The exact computational and manual boundaries are listed in [MANUAL_STEPS.md](MANUAL_STEPS.md). Figures 3–5 have current generators tied to the authenticated 126-genome inputs.

The repository includes no manuscript correspondence, author contact data, credentials, private Earth Engine assets, legacy salinity-proxy scripts, synthetic-data demonstrations, invalid signed-matrix spectral biclustering, or superseded machine-learning analyses.

## Environment

The principal revision analyses were run in more than one recorded software environment. See [SOFTWARE_ENVIRONMENTS.md](SOFTWARE_ENVIRONMENTS.md) before reproducing a workflow. `requirements.txt` is a convenience installation list, not a claim that every historical workflow used one identical environment.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

External commands used by selected workflows include FastTree 2.1.11, Poppler utilities, `qpdf`, and `tectonic`. Google Earth Engine scripts require an authorized Earth Engine account and project.

## Running analyses

Place the authenticated inputs at the paths shown in [DATA_INPUTS.md](DATA_INPUTS.md), then run the relevant script from the repository root unless its header states otherwise. Current revision scripts accept explicit input or output arguments where implemented and stop when required data are absent. Output-producing revision scripts use timestamped or versioned names.

Examples:

```bash
python ISCIENCE_REVISION_20260711/gee_validation/run_exact_id_gee_correlation_validation_20260712_072151.py
python ISCIENCE_REVISION_20260711/aef/extract_exact_id_aef_embeddings_20260718_224936.py --manifest ISCIENCE_REVISION_20260711/integrity/reconciled_analysis_manifest_20260711_110650.csv --year 2024 --validate-only
python ISCIENCE_REVISION_20260711/aef/recompute_full_aef_pfam_analysis_20260718_222823.py --help
python ISCIENCE_REVISION_20260711/figure3_126/rebuild_figure3_recorded_metadata_20260718_223224.py --help
python ISCIENCE_REVISION_20260711/analysis_stats/run_robustness_20260711_085930.py --help
python ISCIENCE_REVISION_20260711/figure4_126/rebuild_figure4_126_20260715_232158.py --help
python ISCIENCE_REVISION_20260711/pf00092_midas/test_pf00092_midas_20260715_232919.py --help
python ISCIENCE_REVISION_20260711/supplement/build_corrected_table_s3_20260718_224354.py --help
```

## Provenance and citation

The timestamped inventory under `provenance/` records the SHA-256 digest of every released script and identifies whether it is a current analysis, an integrity utility, or figure/table support. Software and database references are listed in [SOFTWARE_CITATIONS.md](SOFTWARE_CITATIONS.md).

Please cite this repository using `CITATION.cff` and cite the manuscript when its bibliographic record is available.

## License

Author-owned code is released under the MIT License. External data, databases, imagery, software, and manuscript assets retain their original licenses and terms.

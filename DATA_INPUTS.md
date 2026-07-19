# Data inputs and external services

All analyses use observed or experimentally derived inputs. Scripts stop when required files are absent or fail authentication; they do not generate substitute scientific data. Large inputs and generated scientific outputs are not committed to this code-only repository.

## Public sequence resources

- UAE sequencing reads: NCBI BioProject `PRJNA929663`.
- UAE genome assemblies: Zenodo DOI [10.5281/zenodo.7758508](https://doi.org/10.5281/zenodo.7758508).
- Pfam models and accession metadata: InterPro/Pfam. The PF00092 workflow pins model `PF00092.35` and records the retrieved model hash.
- Structural controls for the PF00092 workflow: RCSB PDB entries `1SHU` and `1AO3`.

## Core authenticated inputs

Paths are relative to the repository root and mirror the analysis workspace.

| Input | Expected path or role |
|---|---|
| Raw HMMER table output | `AlphaEarth/TAGGED_HMMsearch-raw-out/` and `AF3/transfer_hmmsearch_tblout/` |
| Protein FASTA files | `AF3/base_seqs/` |
| Reconstructed raw Pfam matrix | `ISCIENCE_REVISION_20260711/analysis_stats/reconstructed_raw_pfam_counts_20260711_131706.csv.gz` |
| Reconciled cohort manifest | `ISCIENCE_REVISION_20260711/integrity/reconciled_analysis_manifest_20260711_110650.csv` |
| Peptide denominator manifest | `ISCIENCE_REVISION_20260711/analysis_stats/final_peptide_denominator_manifest_20260711_131706.csv` |
| Topology/phylogenetic eigenvectors | `ISCIENCE_REVISION_20260711/analysis_stats/topology_phylogenetic_eigenvectors_20260711_131706.csv` |
| Coordinate-confidence table | `ISCIENCE_REVISION_20260711/spatial/coordinate_confidence_audit_20260711_105047.csv` |
| Authoritative metadata | `AlphaEarth/CSV/Metadata_Table_macroalgae-published.csv` |
| Archived reported AEF embeddings | `AlphaEarth/CSV/alphaearth_embeddings_20251019_122918.csv` |
| Historical AEF extractor input | `pfam_counts_with_metadata_20251019.csv`, in the historical extractor's working directory |
| Retained rbcL alignment/tree | `DATA/rbcL_alignment.fa` and `DATA/rbcL_phylogenetic_tree_20251118_173538.nwk` |
| Genome/phylogeny crosswalk | `MACROALGAE_PHYLOGENIES.csv` |

## Primary GEE workflow products

The following are generated in order and then consumed by the sensitivity, annotation, Figure 5, and table builders:

| Product | Expected path |
|---|---|
| Exact-ID environmental table | `ISCIENCE_REVISION_20260711/gee_validation/exact_id_gee_environmental_extraction_20260712_071838.csv` |
| Environmental extraction manifest | `ISCIENCE_REVISION_20260711/gee_validation/exact_id_gee_environmental_extraction_manifest_20260712_071838.json` |
| Full 87,123-pair discovery table | `ISCIENCE_REVISION_20260711/gee_validation/exact_id_gee_raw_pfam_correlations_20260712_072151.csv.gz` |
| Discovery-supported 84-pair table | `ISCIENCE_REVISION_20260711/gee_validation/exact_id_gee_raw_pfam_fdr05_20260712_072151.csv` |
| Discovery manifest | `ISCIENCE_REVISION_20260711/gee_validation/exact_id_gee_correlation_validation_manifest_20260712_072151.json` |
| Initial selected-family sensitivity table | `ISCIENCE_REVISION_20260711/gee_validation/GEE_primary_selected84_sensitivity_20260719_205317.csv` |
| Initial selected-family null/candidate tables | `GEE_primary_selected84_structured_null_20260719_205317.csv` and `GEE_primary_selected84_candidate_summary_20260719_205317.csv` in the same directory |
| Refined 99,999-permutation null | `ISCIENCE_REVISION_20260711/gee_validation/GEE_primary_selected84_structured_null_refined99999_20260719_210049.csv` |
| Refined 84-pair candidate table | `ISCIENCE_REVISION_20260711/gee_validation/GEE_primary_selected84_candidate_summary_refined99999_20260719_210049.csv` |
| Refined run manifest | `ISCIENCE_REVISION_20260711/gee_validation/GEE_primary_selected84_manifest_refined99999_20260719_210049.json` |

Pinned filenames identify the authenticated manuscript run. To reproduce with a new run ID, either retain the corresponding timestamped outputs or update downstream input arguments/constants and their expected hashes transparently.

## Secondary AEF products

The historical extractor and explicit-year extractor are intentionally distinct:

- The reported archived embeddings came from `aef/archived_reported_extractor/extract_aef_embeddings_20251019.py`, which mosaics the unfiltered annual collection and does not select a year.
- `aef/extract_exact_id_aef_embeddings_20260718_224936.py` requires an explicit `--year`; any output from it is a new year-specific dataset and does not retroactively assign a year to the archived table.

Secondary downstream products include the full corrected AEF run directory `ISCIENCE_REVISION_20260711/aef/full_aef_corrected_run_20260718_154213/`, the 832-pair AEF–GEE alignment tables/manifests under `gee_validation/`, and corrected selected-AEF null/candidate tables under `analysis_stats/`. These are generated scientific results and remain outside the public code repository.

## Supplemental-workbook inputs

The Table S2/S3 and Supplemental Information builders consume prior versioned workbooks, completed analysis tables, integrity JSON files, and the final source DOCX. Those submission artifacts are not code and are not redistributed here. The retained builder chain and pinned SHA-256 checks document the exact dependency sequence. Figure S4 is generated by the full AEF workflow; the final supplemental builder verifies and copies its PDF/SVG byte-for-byte.

## Google Earth Engine

The Earth Engine workflows use:

- `NASA/OCEANDATA/MODIS-Aqua/L3SMI`
- `NOAA/NGDC/ETOPO1`
- `GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL`

Authenticate with the Earth Engine CLI and supply a project to which the user has access. Credentials and private asset identifiers must remain outside the repository.

```bash
earthengine authenticate
```

## Network-dependent annotations

`fetch_gee_robust_interpro_20260719_212508.py` reads the authenticated refined candidate table, selects the 30 Pfams represented by the 49 all-seven pairs, and records every InterPro response and response hash. Missing annotations remain `annotation not available`. The PF00092 workflow also accesses InterPro/Pfam and RCSB PDB endpoints.

## Data rights

Code licensing does not grant redistribution rights for third-party sequence data, satellite products, databases, fonts, or submission assets. Obtain each input under the source's current terms.

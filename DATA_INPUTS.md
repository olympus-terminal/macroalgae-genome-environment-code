# Data inputs and external services

All analyses use observed or experimentally derived inputs. The scripts stop when required files are absent; they do not generate substitute scientific data.

## Public sequence resources

- UAE sequencing reads: NCBI BioProject `PRJNA929663`.
- UAE genome assemblies: Zenodo DOI [10.5281/zenodo.7758508](https://doi.org/10.5281/zenodo.7758508).
- Pfam models and accession metadata: InterPro/Pfam. The PF00092 workflow pins model `PF00092.35` and records the retrieved model hash.
- Structural controls for the PF00092 workflow: RCSB PDB entries `1SHU` and `1AO3`.

## Expected local inputs

The following large or submission-associated inputs are not committed here. Paths are relative to the repository root and mirror the analysis layout.

| Input | Expected path or role |
|---|---|
| Raw HMMER table output | `AlphaEarth/TAGGED_HMMsearch-raw-out/` and `AF3/transfer_hmmsearch_tblout/` |
| Protein FASTA files | `AF3/base_seqs/` |
| Reconstructed raw Pfam matrix | `ISCIENCE_REVISION_20260711/analysis_stats/reconstructed_raw_pfam_counts_20260711_131706.csv.gz` |
| Reconciled cohort manifest | `ISCIENCE_REVISION_20260711/integrity/reconciled_analysis_manifest_20260711_110650.csv` |
| Coordinate-confidence table | `ISCIENCE_REVISION_20260711/spatial/coordinate_confidence_audit_20260711_105047.csv` |
| Archived AEF embeddings used in the reported analysis | `AlphaEarth/CSV/alphaearth_embeddings_20251019_122918.csv` |
| Retained rbcL alignment | `DATA/rbcL_alignment.fa` |
| Retained rbcL tree | `DATA/rbcL_phylogenetic_tree_20251118_173538.nwk` |
| Genome/phylogeny crosswalk | `MACROALGAE_PHYLOGENIES.csv` |
| Revised source workbooks | submission Tables S1–S3 at the paths expected by the builders |

Input and output hashes from the manuscript revision runs are recorded in the corresponding run manifests supplied with the revision package. This repository’s `provenance/` directory hashes the released code itself.

## Google Earth Engine

The exact-ID GEE workflow uses:

- `NASA/OCEANDATA/MODIS-Aqua/L3SMI`
- `NOAA/NGDC/ETOPO1`
- `GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL`

Authenticate with the Earth Engine CLI and supply a project for which the user has access. Credentials and private asset identifiers must remain outside the repository.

```bash
earthengine authenticate
```

The exact-ID environmental extraction retains source-product missingness at the recorded coordinates. The corrected AEF extractor selects the 126 safe rows from the SHA-256-pinned reconciled manifest, requires one explicit calendar year, preserves canonical genome IDs, and samples all 64 axes at a 10-m request scale. Its output is a new year-specific dataset; it does not retroactively assign a year to the archived embedding table used in the reported analysis.

## Network-dependent annotation workflows

`fetch_priority_interpro_20260711_114039.py` records the InterPro API response for each requested accession. The PF00092 workflow also accesses InterPro/Pfam and RCSB PDB endpoints and records unavailable annotations explicitly.

## Data rights

Code licensing does not grant redistribution rights for third-party sequence data, satellite products, databases, fonts, or figure assets. Obtain each input from its cited source under the source’s current terms.

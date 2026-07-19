# Computational and manual figure boundaries

This file records where final manuscript artwork was assembled outside a single analysis script.

| Item | Computational source | Final assembly |
|---|---|---|
| Figure 1 | Sampling-map and layout scripts under `ISCIENCE_REVISION_20260711/spatial/` and `ISCIENCE_REVISION_20260711/figures/` | Specimen panels and final layout were assembled manually. |
| Figure 2 | Retained rbcL alignment/tree workflow | Tree annotation and the standalone phylogeny display were completed in iTOL. |
| Figure 3 | `ISCIENCE_REVISION_20260711/figure3_126/rebuild_figure3_recorded_metadata_20260718_223224.py` | Generated from exact-ID reconstructed Pfam counts and recorded metadata for 126 genomes. |
| Figure 4 | `ISCIENCE_REVISION_20260711/figure4_126/rebuild_figure4_126_20260715_232158.py` | Generated computationally from the authenticated 126-genome inputs. |
| Figure 5 | `ISCIENCE_REVISION_20260711/figure5_gee_primary/build_figure5_primary_gee_sensitivity_20260719_210207.py` | Generated from the authenticated 84-pair primary GEE sensitivity and refined 99,999-permutation null tables; all 84 pairs are shown. |
| Figure S1 | Retained source artwork | The exact final-panel generator was not retained. |
| Figure S3 | `ISCIENCE_REVISION_20260711/supplement/build_supplemental_legends_v8_20260719_214546.py` | Generated from the retained AEF non-null results and corrected 99,999-permutation null/candidate tables; the builder validates displayed values against those inputs. |
| Figure S4 | `ISCIENCE_REVISION_20260711/aef/recompute_full_aef_pfam_analysis_20260718_222823.py` | The complete 10,707-Pfam × 64-axis AEF correlation-profile landscape is generated computationally and copied byte-for-byte by the final supplemental builder. |
| Tables S2/S3 | Versioned builders under `ISCIENCE_REVISION_20260711/supplement/` | Completed scientific tables and prior audited workbook versions are consumed in the documented chain; each builder writes an integrity record and preserves prior versions. |

The starting rbcL retrieval and alignment commands were not retained. The repository therefore supplies the tree-building script but does not present an unrecoverable retrieval/alignment history as a fully automated workflow.

Editorial manuscript and response-letter assembly utilities are outside this scientific-code release. They alter submission documents but do not generate or validate scientific results.

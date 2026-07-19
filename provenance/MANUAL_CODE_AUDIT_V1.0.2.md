# v1.0.2 manual scientific-code audit

The final timestamped automated release audit in this directory passed syntax, required-file, private-path, credential, file-type/size, and prohibited synthetic-random-API checks. The remaining flagged calls were inspected manually:

| File | Flag | Disposition |
|---|---|---|
| `analysis_stats/run_robustness_20260711_085930.py` | `np.random.Generator` / `default_rng` | Used for bootstrap resampling and permutation of loaded observed units; the seed is recorded. |
| `analysis_stats/correct_aef_structured_null_conditional_tail_20260719_205621.py` | `SeedSequence` / `default_rng` | Used only to permute intact observed coordinate sites within observed phylum-composition strata. |
| `gee_validation/run_gee_primary_sensitivity_20260719_203058.py` | `default_rng` | Passed to observed-site bootstrap and structured-permutation functions; no values are generated as observations. |
| `gee_validation/run_gee_structured_null_refinement_20260719_205905.py` | `default_rng` | Permutes observed environmental ranks among observed sites within strata. |
| `figure4_126/rebuild_figure4_126_20260715_232158.py` | `np.linspace` | Creates histogram bins from loaded-data minima/maxima and symmetric plot ticks from a computed axis limit; it does not create scientific values. |

No `torch.rand*`, unscoped `np.random.rand*`, dummy-value generator, placeholder-result return, private dataset, credential, or generated scientific result file is present. Figure/table builders read completed authenticated tables; fixed counts and caption values act as checked expectations and presentation text, not substitutes for analysis.

Read-only trace checks against the retained workspace inputs confirmed that the released workflows point to the materialized 84-pair selected GEE family, 56 structured-null-supported pairs, 49 all-seven pairs spanning 30 Pfams, all 13 Bonferroni pairs retained, the 832-pair AEF–GEE alignment (224 BH and 44 Bonferroni associations), the two corrected all-seven AEF pairs, and 30 successful InterPro/Pfam records. These scientific result files were used for verification only and are not included in this repository.

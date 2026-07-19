# Scientific integrity rules for this release

- Scientific outputs must be computed from loaded observational or experimental data.
- Synthetic, simulated, randomly generated, or reconstructed summary values must not substitute for missing publication data.
- Random-number generation is limited to documented resampling, permutation, or coordinate-jitter sensitivity analyses using observed inputs and recorded seeds.
- Annotation names and functions must come from cited source records. Unavailable annotations remain unavailable.
- Scripts must fail on missing or invalid required inputs instead of returning plausible placeholder results.
- Numerical outputs must retain input, parameter, software, and code provenance where implemented.
- The authenticated analysis cohort contains 126 genomes; cohort joins are enforced by genome identifier.
- Expected result counts embedded in figure/table builders are assertions checked against loaded authenticated result tables; they are not substitutes for analysis.
- The historical unfiltered AEF extractor is retained verbatim and explicitly labeled archival. The later explicit-year extractor is a separate prospective workflow.

The timestamped release audit under `provenance/` records automated checks for absolute private paths, credential-like strings, syntax errors, oversized/private-input file types, and prohibited synthetic-data patterns. Random APIs are inventoried for manual confirmation that they are limited to observed-data resampling or permutation.

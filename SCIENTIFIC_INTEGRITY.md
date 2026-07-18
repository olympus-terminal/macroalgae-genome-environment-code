# Scientific integrity rules for this release

- Scientific outputs must be computed from loaded observational or experimental data.
- Synthetic, simulated, randomly generated, or reconstructed summary values must not substitute for missing publication data.
- Random-number generation is limited to documented resampling, permutation, or coordinate-jitter sensitivity analyses using observed inputs and recorded seeds.
- Annotation names and functions must come from cited source records. Unavailable annotations remain unavailable.
- Scripts must fail on missing or invalid required inputs instead of returning plausible placeholder results.
- Numerical outputs must retain input, parameter, software, and code provenance where implemented.
- The authenticated analysis cohort contains 126 genomes; cohort joins are enforced by genome identifier.

The timestamped release audit under `provenance/` records automated checks for absolute private paths, credential-like strings, syntax errors, and prohibited synthetic-data patterns.

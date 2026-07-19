# Recorded software environments

Exact versions are reported only where they were retained in local run provenance.

## Revision robustness workflow

Canonical run `20260711_131706` recorded:

- Python 3.9.23
- NumPy 2.0.2
- pandas 2.3.3
- SciPy 1.13.1
- statsmodels 0.14.5
- Biopython 1.85

## Exact-ID GEE extraction and correlations

The extraction run recorded Python 3.10.9 and Earth Engine Python API 1.6.12. The correlation run recorded:

- Python 3.10.9
- NumPy 1.22.0
- pandas 1.5.3
- SciPy 1.11.4
- statsmodels 0.14.5

## Primary GEE sensitivity and refined structured null

The completed selected-family runs recorded Python, NumPy, pandas, SciPy, and statsmodels versions in their JSON manifests. The final workflow uses the same version families as the exact-ID correlation and revision-robustness environments above. Random operations are restricted to observed-site bootstrap resampling and intact-site structured permutations, with seed `20260719` recorded by the scripts.

The final Figure 5 builder additionally uses Matplotlib and exports vector PDF/SVG with 6-point Arial text, 0.25-point strokes, and a transparent background.

## PF00092 MIDAS workflow

The retained PF00092 run recorded:

- Python 3.10.9
- pandas 1.5.3
- PyHMMER 0.11.0

The workflow pins Pfam model `PF00092.35` and verifies retrieved resources by SHA-256.

## Figure 4 workflow

The retained Figure 4 run recorded:

- Python 3.10.9
- NumPy 1.22.0
- pandas 1.5.3
- SciPy 1.11.4
- statsmodels 0.14.5
- Matplotlib 3.8.4

## Corrected AEF screen and Figure S4

Canonical corrected run `20260718_154213` recorded:

- Python 3.10.9
- NumPy 1.22.0
- pandas 1.5.3
- SciPy 1.11.4
- statsmodels 0.14.5
- Matplotlib 3.8.4

In the final manuscript hierarchy, this AEF landscape is Figure S4; the computation and canonical run are unchanged.

## AEF–GEE alignment and corrected selected-AEF null

The same-site AEF–GEE alignment uses Python, NumPy, pandas, SciPy, and statsmodels. The corrected selected-AEF null additionally uses NumPy's random generator only for permutation of observed coordinate-site units within observed phylum-composition strata and records the seed and software versions in its manifest.

## Corrected Figure 3

Canonical corrected run `20260718_223934` recorded the same Python, NumPy, pandas, SciPy, statsmodels, and Matplotlib versions as the corrected AEF screen.

## Spatial figure workflow

The retained spatial manifests recorded:

- Cartopy 0.23.0
- Matplotlib 3.9.4
- NumPy 2.0.2
- pandas 2.3.3
- SciPy 1.13.1
- Earth Engine Python API 1.6.12 for the network-dependent re-extraction script

## Phylogenetic tree

The retained tree-construction log records FastTree 2.1.11, double precision, with JTT+CAT20 settings and SH-like local support.

## Corrected Table S3 builder

The retained builder uses Python 3.10.9, NumPy 1.22.0, pandas 1.5.3, SciPy 1.11.4, statsmodels 0.14.5, and openpyxl 3.1.5.

## Final supplemental builders

The final Table S2/S3 and Supplemental Information builders use Python 3, openpyxl, pandas, NumPy, SciPy, statsmodels, Matplotlib, python-docx, and lxml as applicable. Rendering and validation additionally call LibreOffice, Poppler (`pdfinfo`/`pdffonts`), `qpdf`, and Tectonic. These builders consume completed authenticated result tables and do not recompute inferential analyses.

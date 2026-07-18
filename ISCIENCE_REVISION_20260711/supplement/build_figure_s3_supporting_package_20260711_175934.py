#!/usr/bin/env python3
"""Build a traceable Figure S3 package from the completed robustness run.

This program does not recompute or alter scientific results. It reads the
canonical 20260711_131706 outputs, validates the selected-set scope, copies the
existing vector figure without modification, and creates caption/provenance
materials plus a one-page vector PDF. No synthetic data are used.

Builder version: 2026-07-11.1
Required software for the one-page PDF: Tectonic 0.15 or compatible.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


PACKAGE_ID = "20260711_175934"
SOURCE_RUN_ID = "20260711_131706"
BUILDER_VERSION = "2026-07-11.1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def one(rows: list[dict[str, str]], **criteria: str) -> dict[str, str]:
    hits = [row for row in rows if all(row.get(key) == value for key, value in criteria.items())]
    if len(hits) != 1:
        raise ValueError(f"Expected one row for {criteria}, found {len(hits)}")
    return hits[0]


def tex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
        "≥": ">=",
        "−": "--",
        "–": "--",
        "×": "x",
        "ρ": "rho",
    }
    return "".join(replacements.get(char, char) for char in text)


def fmt(value: str, digits: int = 3) -> str:
    return f"{float(value):.{digits}f}"


def validate_inputs(root: Path) -> dict[str, object]:
    stats = root / "ISCIENCE_REVISION_20260711" / "analysis_stats"
    paths = {
        "source_figure_pdf": stats / f"Figure_Robustness_Sensitivity_{SOURCE_RUN_ID}.pdf",
        "source_figure_svg": stats / f"Figure_Robustness_Sensitivity_{SOURCE_RUN_ID}.svg",
        "run_manifest": stats / f"run_manifest_{SOURCE_RUN_ID}.json",
        "candidate_table": stats / f"robust_candidate_table_{SOURCE_RUN_ID}.csv",
        "aef_results": stats / f"archived_AEF_priority_robustness_{SOURCE_RUN_ID}.csv",
        "aef_null": stats / f"archived_AEF_within_phylum_null_{SOURCE_RUN_ID}.csv",
        "quantitative_summary": stats / f"manuscript_ready_quantitative_summary_{SOURCE_RUN_ID}.md",
        "analysis_script": stats / "run_robustness_20260711_085930.py",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required inputs: " + ", ".join(missing))

    manifest = json.loads(paths["run_manifest"].read_text(encoding="utf-8"))
    if manifest.get("status") != "complete":
        raise ValueError("Source robustness run is not marked complete")
    parameters = manifest.get("parameters", {})
    expected_parameters = {"permutations": 9999, "bootstraps": 999, "seed": 20260711}
    if any(parameters.get(key) != value for key, value in expected_parameters.items()):
        raise ValueError(f"Unexpected source parameters: {parameters}")
    counts = manifest.get("data_counts", {})
    expected_counts = {
        "canonical_genomes": 126,
        "unique_coordinate_sites": 90,
        "repeated_coordinate_sites": 19,
        "max_genomes_per_coordinate_site": 7,
        "traceable_peptide_denominators": 121,
        "tree_mapped_AEF_genomes": 112,
    }
    if any(counts.get(key) != value for key, value in expected_counts.items()):
        raise ValueError(f"Unexpected source data counts: {counts}")

    candidates = read_csv(paths["candidate_table"])
    if len(candidates) != 68:
        raise ValueError(f"Expected 68 selected pairs, found {len(candidates)}")
    pf00092_axes = {row["latent_axis"] for row in candidates if row["pfam"] == "PF00092"}
    if pf00092_axes != {f"A{i:02d}" for i in range(64)}:
        raise ValueError("PF00092 selected set does not contain all 64 AEF axes")
    non_pf00092 = {(row["pfam"], row["latent_axis"]) for row in candidates if row["pfam"] != "PF00092"}
    expected_priority = {
        ("PF13411", "A18"),
        ("PF01638", "A52"),
        ("PF01638", "A53"),
        ("PF10988", "A36"),
    }
    if non_pf00092 != expected_priority:
        raise ValueError(f"Unexpected highlighted-pair set: {sorted(non_pf00092)}")
    robust = {
        (row["pfam"], row["latent_axis"])
        for row in candidates
        if row["robust_candidate_all_required_checks"] == "True"
    }
    if robust != {("PF01638", "A52"), ("PF01638", "A53")}:
        raise ValueError(f"Unexpected all-check robust set: {sorted(robust)}")

    results = read_csv(paths["aef_results"])
    nulls = read_csv(paths["aef_null"])
    raw_site = sorted(
        (
            row for row in results
            if row["representation"] == "raw_count" and row["method"] == "site_mean_spearman"
        ),
        key=lambda row: float(row["p_value"]),
    )
    displayed = [row for row in raw_site if row["pfam"] != "PF00092"]
    displayed += [row for row in raw_site if row["pfam"] == "PF00092"][:8]
    displayed = displayed[:12]
    if len(displayed) != 12 or len({(row["pfam"], row["latent_axis"]) for row in displayed}) != 12:
        raise ValueError("Could not reconstruct the 12-pair display subset")

    key_rows: dict[str, dict[str, str]] = {}
    key_nulls: dict[str, dict[str, str]] = {}
    for axis in ("A52", "A53"):
        key_rows[axis] = one(
            results,
            pfam="PF01638",
            latent_axis=axis,
            representation="per_total_pfam_hit",
            method="site_mean_spearman",
        )
        key_nulls[axis] = one(
            nulls,
            pfam="PF01638",
            latent_axis=axis,
            representation="per_total_pfam_hit",
        )
    return {
        "root": root,
        "paths": paths,
        "manifest": manifest,
        "counts": counts,
        "candidates": candidates,
        "displayed": displayed,
        "key_rows": key_rows,
        "key_nulls": key_nulls,
    }


def build_caption(validated: dict[str, object]) -> str:
    key_rows = validated["key_rows"]
    key_nulls = validated["key_nulls"]
    a52 = key_rows["A52"]
    a53 = key_rows["A53"]
    n52 = key_nulls["A52"]
    n53 = key_nulls["A53"]
    return (
        "Figure S3. Post hoc selected-set sensitivity analysis of archived Pfam–AEF "
        "latent-feature associations. (A) Rank-based effects for 12 displayed priority pairs: "
        "the four previously highlighted discovery pairs and eight PF00092 axes selected for "
        "display by the smallest raw-count site-mean P values. These 12 are a display subset of "
        "68 pairs tested (four highlighted pairs plus PF00092 against all 64 AEF axes). R, raw "
        "count; T, count per total Pfam hit; P, count per final peptide record. Qual, "
        "coordinate-clustered rank model with continuous quality/search-space covariates and "
        "phylum; B50, BUSCO ≥50% site-mean sensitivity; Phy, broken-stick topological "
        "eigenvectors plus coordinate-clustered covariance. Stars denote selected-set "
        "Benjamini–Hochberg q<0.05 within each representation–method family across all 68 pairs. "
        "(B) Observed total-hit site-mean ρ and the 95% interval from 9,999 permutations of intact "
        "sites within site phylum-composition strata; red points have nominal empirical P<0.05. "
        "(C) Number of the 12 displayed pairs with nominal P<0.05 per check. Only PF01638/A52 "
        "and PF01638/A53 met direction consistency and all seven required selected-set gates "
        "(7/7). Their total-hit site-mean effects were ρ=" + fmt(a52["effect"]) +
        " (95% bootstrap CI " + fmt(a52["ci_low"]) + "–" + fmt(a52["ci_high"]) +
        "; q=" + f"{float(a52['selected_AEF_set_bh_q']):.3g}" + ") and ρ=" +
        fmt(a53["effect"]) + " (95% CI " + fmt(a53["ci_low"]) + " to " +
        fmt(a53["ci_high"]) + "; q=" + f"{float(a53['selected_AEF_set_bh_q']):.3g}" +
        "). Both structured-null tests had empirical P=" +
        f"{float(n52['empirical_two_sided_p']):.4f}" + " and selected-set q=" +
        f"{float(n52['selected_AEF_empirical_bh_q']):.4f}" +
        ". PF10988/A36 (6/7) and PF13411/A18 (5/7) did not pass all gates. AEF axes "
        "are unitless latent descriptors; these post hoc selected-set results do not decode "
        "physical variables or establish adaptation or mechanism."
    )


def build_method_note(validated: dict[str, object], caption: str) -> str:
    root = validated["root"]
    paths = validated["paths"]
    manifest = validated["manifest"]
    counts = validated["counts"]
    displayed = validated["displayed"]
    displayed_text = ", ".join(f"{row['pfam']}/{row['latent_axis']}" for row in displayed)
    command = " ".join(str(value) for value in manifest["command"])
    input_lines = "\n".join(
        f"- {name}: {path.resolve().relative_to(root)}\n  SHA-256: {sha256(path)}"
        for name, path in paths.items()
    )
    return f"""FIGURE S3 CAPTION AND METHOD NOTE

ANALYSIS STATUS
POST HOC; SELECTED SET. This is a sensitivity analysis of archived candidate associations, not an independent confirmation, a discovery-wide rerun, or evidence of causality.

CAPTION
{caption}

METHOD/INTERPRETATION NOTE
- Cohort: {counts['canonical_genomes']} genomes at {counts['unique_coordinate_sites']} exact coordinate sites; {counts['repeated_coordinate_sites']} sites were repeated (maximum {counts['max_genomes_per_coordinate_site']} genomes per site).
- Tested set: 68 pairs, comprising four previously highlighted Pfam–AEF pairs and PF00092 against all 64 AEF axes.
- Display subset (12): {displayed_text}.
- Representations: raw HMM-search counts; counts divided by total Pfam hits; and counts divided by final peptide records. The peptide-record denominator was directly traceable for {counts['traceable_peptide_denominators']}/{counts['canonical_genomes']} genomes and is a search-space proxy, not a documented predicted-gene count.
- Site analysis: Spearman rank correlation of Pfam site means with the site AEF vector. Formula: rho = cor(rank(mean_site(Pfam representation)), rank(AEF axis)).
- Quality analysis: rank-based regression with continuous assembly size, BUSCO completeness, total Pfam hits, peptide-search-space size, and phylum, using covariance clustered by exact coordinate site.
- BUSCO sensitivity: site-mean analysis after BUSCO completeness ≥50%; a ≥70% threshold was also computed but was not one of the seven required gates.
- Topology sensitivity: rank model with outcome-blind broken-stick topological eigenvectors and coordinate-clustered covariance. The tree-mapped AEF set contained {counts['tree_mapped_AEF_genomes']} genomes. This is not Gaussian PGLS.
- Structured empirical null: 9,999 permutations of intact site AEF vectors within site phylum-composition strata. Formula: P_empirical = (1 + number(|rho_perm| ≥ |rho_observed|)) / (9,999 + 1).
- Intervals: 999 bootstrap resamples of unique site means; random seed 20260711. Randomness was used only for resampling real observations.
- Multiple testing: Benjamini–Hochberg adjustment within each representation–method family across the 68 post hoc selected pairs. This selected-set q value is not discovery-wide error control.
- Seven required gates: selected-set q<0.05 for raw, total-hit, and peptide-record site means; total-hit quality model; total-hit BUSCO≥50% site mean; total-hit broken-stick topology model; and total-hit structured empirical null. The effect direction also had to agree across all estimable checks.
- GEE variables were not included because the retained raw export did not overlap the 126-row derived table; publication inference remains gated pending reproducible re-extraction.
- Interpretation boundary: AEF axes are unitless latent geospatial descriptors. Passing the selected-set checks supports statistical stability within this post hoc set only; it does not identify a physical environmental driver or establish selection, adaptation, or molecular mechanism.

SOURCE ANALYSIS
Command: {command}
Source run status: {manifest['status']}
Source run ID: {manifest['run_id']}
Source script version: {manifest['script_version']}

INPUTS AND SHA-256
{input_lines}
"""


def build_tex(caption: str, figure_filename: str) -> str:
    return rf"""\documentclass[10pt]{{article}}
\usepackage[T1]{{fontenc}}
\usepackage{{helvet}}
\renewcommand{{\familydefault}}{{\sfdefault}}
\usepackage[a4paper,margin=12mm]{{geometry}}
\usepackage{{graphicx}}
\usepackage{{microtype}}
\usepackage{{ragged2e}}
\pagestyle{{empty}}
\setlength{{\parindent}}{{0pt}}
\setlength{{\parskip}}{{3pt}}
\setlength{{\emergencystretch}}{{1em}}
\begin{{document}}
\textbf{{Figure S3}}\\[-1pt]
\begin{{center}}
\includegraphics[width=\textwidth,height=0.58\textheight,keepaspectratio]{{{tex_escape(figure_filename)}}}
\end{{center}}
\vspace{{-5pt}}
\small\RaggedRight {tex_escape(caption)}
\end{{document}}
"""


def pdf_page_count(path: Path) -> int:
    result = subprocess.run(
        ["pdfinfo", str(path)], check=True, capture_output=True, text=True
    )
    for line in result.stdout.splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    raise ValueError(f"Could not determine page count for {path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    validated = validate_inputs(root)
    caption = build_caption(validated)
    method_note = build_method_note(validated, caption)
    if args.validate_only:
        print(json.dumps({
            "status": "validation_complete_outputs_not_written",
            "source_run_id": SOURCE_RUN_ID,
            "selected_pairs": len(validated["candidates"]),
            "displayed_pairs": len(validated["displayed"]),
            "all_check_robust_pairs": ["PF01638/A52", "PF01638/A53"],
        }, indent=2))
        return 0

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_paths = {
        "submission_figure_pdf": output_dir / f"Figure_S3_posthoc_selected_set_{PACKAGE_ID}.pdf",
        "submission_figure_svg": output_dir / f"Figure_S3_posthoc_selected_set_{PACKAGE_ID}.svg",
        "caption_method_note": output_dir / f"Figure_S3_caption_method_note_{PACKAGE_ID}.txt",
        "supporting_page_tex": output_dir / f"Figure_S3_supporting_page_{PACKAGE_ID}.tex",
        "supporting_page_pdf": output_dir / f"Figure_S3_supporting_page_{PACKAGE_ID}.pdf",
        "integrity_manifest": output_dir / f"Figure_S3_supporting_integrity_{PACKAGE_ID}.json",
    }
    existing = [str(path) for path in output_paths.values() if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing package files: " + ", ".join(existing))

    source_paths = validated["paths"]
    shutil.copyfile(source_paths["source_figure_pdf"], output_paths["submission_figure_pdf"])
    shutil.copyfile(source_paths["source_figure_svg"], output_paths["submission_figure_svg"])
    output_paths["caption_method_note"].write_text(method_note, encoding="utf-8")
    output_paths["supporting_page_tex"].write_text(
        build_tex(caption, output_paths["submission_figure_pdf"].name), encoding="utf-8"
    )

    subprocess.run(
        [
            "tectonic",
            "--only-cached",
            "--chatter",
            "minimal",
            "--outdir",
            str(output_dir),
            str(output_paths["supporting_page_tex"]),
        ],
        check=True,
    )
    subprocess.run(["qpdf", "--check", str(output_paths["submission_figure_pdf"])], check=True)
    subprocess.run(["qpdf", "--check", str(output_paths["supporting_page_pdf"])], check=True)
    if pdf_page_count(output_paths["submission_figure_pdf"]) != 1:
        raise ValueError("Submission Figure S3 is not one page")
    if pdf_page_count(output_paths["supporting_page_pdf"]) != 1:
        raise ValueError("Figure S3 supporting page is not one page")
    if sha256(output_paths["submission_figure_pdf"]) != sha256(source_paths["source_figure_pdf"]):
        raise ValueError("Copied PDF differs from canonical source")
    if sha256(output_paths["submission_figure_svg"]) != sha256(source_paths["source_figure_svg"]):
        raise ValueError("Copied SVG differs from canonical source")

    manifest_path = output_paths.pop("integrity_manifest")
    integrity = {
        "package_id": PACKAGE_ID,
        "builder_version": BUILDER_VERSION,
        "builder": str(Path(__file__).resolve().relative_to(root)),
        "builder_sha256": sha256(Path(__file__).resolve()),
        "generated_at": datetime.now().astimezone().isoformat(),
        "source_run_id": SOURCE_RUN_ID,
        "analysis_status": "post hoc; selected set; not discovery-wide",
        "scientific_data_recomputed": False,
        "synthetic_or_placeholder_data_used": False,
        "source_parameters": validated["manifest"]["parameters"],
        "source_counts": validated["counts"],
        "selected_pairs_tested": len(validated["candidates"]),
        "pairs_displayed": len(validated["displayed"]),
        "displayed_pairs": [
            f"{row['pfam']}/{row['latent_axis']}" for row in validated["displayed"]
        ],
        "all_check_robust_pairs": ["PF01638/A52", "PF01638/A53"],
        "inputs": {
            name: {
                "path": str(path.resolve().relative_to(root)),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
            for name, path in source_paths.items()
        },
        "outputs": {
            name: {
                "path": str(path.resolve().relative_to(root)),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
            for name, path in output_paths.items()
        },
        "validations": {
            "candidate_rows_equal_68": True,
            "pf00092_axes_equal_all_64_AEF_axes": True,
            "highlighted_non_PF00092_pairs_equal_expected_four": True,
            "all_check_robust_set_equal_PF01638_A52_A53": True,
            "display_subset_reconstructed_from_source_code_rule": True,
            "copied_pdf_sha256_equals_canonical": True,
            "copied_svg_sha256_equals_canonical": True,
            "submission_figure_pdf_pages": 1,
            "supporting_page_pdf_pages": 1,
            "qpdf_checks_passed": True,
            "vector_source_preserved": True,
        },
        "interpretation_boundary": (
            "AEF axes are unitless latent descriptors; results are post hoc selected-set "
            "stability checks and do not establish a physical driver, adaptation, or mechanism."
        ),
    }
    manifest_path.write_text(json.dumps(integrity, indent=2), encoding="utf-8")
    print(json.dumps({
        "status": "complete",
        "outputs": {name: str(path.resolve()) for name, path in output_paths.items()},
        "integrity_manifest": str(manifest_path.resolve()),
    }, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise

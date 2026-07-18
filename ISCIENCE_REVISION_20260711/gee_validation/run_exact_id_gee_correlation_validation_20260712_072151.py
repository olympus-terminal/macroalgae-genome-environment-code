#!/usr/bin/env python3
"""Correct the original GEE–Pfam correlation workflow using exact genome IDs.

This bounded validation repeats the submitted raw-count discovery specification
after replacing the non-unique species-name merge with an exact 126-ID join.
It uses the newly re-extracted GEE values and retained raw HMM counts. It does
not add new variables, tune a threshold to significance, impute missing values,
or generate synthetic data.

Generated: 2026-07-12
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata, spearmanr, t as student_t
import statsmodels
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "ISCIENCE_REVISION_20260711" / "gee_validation"
ENV = OUTDIR / "exact_id_gee_environmental_extraction_20260712_071838.csv"
ENV_MANIFEST = OUTDIR / "exact_id_gee_environmental_extraction_manifest_20260712_071838.json"
RAW = (
    ROOT
    / "ISCIENCE_REVISION_20260711"
    / "analysis_stats"
    / "reconstructed_raw_pfam_counts_20260711_131706.csv.gz"
)
OUTPUT_ALL = OUTDIR / "exact_id_gee_raw_pfam_correlations_20260712_072151.csv.gz"
OUTPUT_SIGNIFICANT = OUTDIR / "exact_id_gee_raw_pfam_fdr05_20260712_072151.csv"
OUTPUT_SUMMARY = OUTDIR / "exact_id_gee_raw_pfam_summary_20260712_072151.csv"
MANIFEST = OUTDIR / "exact_id_gee_correlation_validation_manifest_20260712_072151.json"

EXPECTED_RAW_SHA256 = "683228342a90d2ecf2897930cd7c147f23973de0444155bfc162d50b09dd22bb"
PREVALENCE_FRACTION = 0.05
MIN_VALID_PAIRS = 10

ENVIRONMENTAL_COLUMNS = [
    "sst_mean_c",
    "sst_max_c",
    "sst_min_c",
    "sst_annual_range_c",
    "sst_summer_c",
    "sst_winter_c",
    "chlorophyll_mean_mg_m3",
    "chlorophyll_max_mg_m3",
    "chlorophyll_std_mg_m3",
    "poc_mean_mg_m3",
    "depth_meters",
    "distance_coast_km",
    "water_clarity_ratio",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    for path in [ENV, ENV_MANIFEST, RAW]:
        if not path.is_file():
            raise FileNotFoundError(path)
    env_manifest = json.loads(ENV_MANIFEST.read_text(encoding="utf-8"))
    if sha256(ENV) != env_manifest["output"]["sha256"]:
        raise RuntimeError("GEE extraction hash differs from its provenance manifest")
    if sha256(RAW) != EXPECTED_RAW_SHA256:
        raise RuntimeError("Raw Pfam matrix hash differs from the validated run")

    env = pd.read_csv(ENV, low_memory=False)
    raw = pd.read_csv(RAW, compression="gzip", low_memory=False)
    if len(env) != 126 or env["genome_id"].nunique() != 126:
        raise RuntimeError("Environmental input is not 126 unique exact genome IDs")
    if len(raw) != 126 or raw["Genome"].nunique() != 126:
        raise RuntimeError("Raw Pfam input is not 126 unique exact genome IDs")
    env_ids = set(env["genome_id"])
    raw_ids = set(raw["Genome"])
    if env_ids != raw_ids:
        raise RuntimeError(
            f"Exact-ID sets differ: missing raw={sorted(env_ids-raw_ids)}, "
            f"missing environment={sorted(raw_ids-env_ids)}"
        )
    raw = raw.set_index("Genome").loc[env["genome_id"]].reset_index()
    if raw["Genome"].tolist() != env["genome_id"].tolist():
        raise RuntimeError("Exact-ID reorder failed")
    return env, raw, env_manifest


def vector_spearman(matrix: np.ndarray, outcome: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return rho, asymptotic two-sided p, and nonconstant mask by column."""
    n = len(outcome)
    if n < MIN_VALID_PAIRS:
        raise ValueError("Insufficient observations")
    ranked_x = rankdata(matrix, axis=0, method="average")
    ranked_y = rankdata(outcome, method="average")
    centered_x = ranked_x - ranked_x.mean(axis=0)
    centered_y = ranked_y - ranked_y.mean()
    sumsq_x = np.einsum("ij,ij->j", centered_x, centered_x)
    sumsq_y = float(np.dot(centered_y, centered_y))
    estimable = (sumsq_x > 0) & (sumsq_y > 0)
    rho = np.full(matrix.shape[1], np.nan, dtype=float)
    rho[estimable] = (centered_x[:, estimable].T @ centered_y) / np.sqrt(
        sumsq_x[estimable] * sumsq_y
    )
    rho[estimable] = np.clip(rho[estimable], -1.0, 1.0)
    p_value = np.full(matrix.shape[1], np.nan, dtype=float)
    interior = estimable & (np.abs(rho) < 1.0)
    statistic = rho[interior] * np.sqrt((n - 2) / ((1.0 + rho[interior]) * (1.0 - rho[interior])))
    p_value[interior] = 2.0 * student_t.sf(np.abs(statistic), df=n - 2)
    p_value[estimable & (np.abs(rho) == 1.0)] = 0.0
    return rho, p_value, estimable


def compute(env: pd.DataFrame, raw: pd.DataFrame) -> tuple[pd.DataFrame, int, list[str]]:
    pfam_columns = [column for column in raw.columns if re.fullmatch(r"PF\d{5}", str(column))]
    if not pfam_columns:
        raise RuntimeError("No Pfam columns in raw matrix")
    raw_values = raw[pfam_columns].apply(pd.to_numeric, errors="raise")
    if raw_values.isna().any().any():
        raise RuntimeError("Raw count matrix contains missing values")
    if (raw_values < 0).any().any():
        raise RuntimeError("Raw count matrix contains negative values")

    prevalence_minimum = math.ceil(PREVALENCE_FRACTION * len(raw_values))
    prevalence_counts = (raw_values > 0).sum(axis=0)
    eligible = prevalence_counts[prevalence_counts >= prevalence_minimum].index.tolist()
    matrix_all = raw_values[eligible].to_numpy(dtype=float, copy=True)
    results: list[pd.DataFrame] = []

    for variable in ENVIRONMENTAL_COLUMNS:
        outcome = pd.to_numeric(env[variable], errors="coerce")
        valid = outcome.notna().to_numpy()
        n = int(valid.sum())
        if n < MIN_VALID_PAIRS:
            continue
        matrix = matrix_all[valid, :]
        y = outcome.to_numpy(dtype=float)[valid]
        rho, p_value, estimable = vector_spearman(matrix, y)
        indices = np.flatnonzero(estimable)
        block = pd.DataFrame(
            {
                "pfam": [eligible[index] for index in indices],
                "env_var": variable,
                "spearman_r": rho[indices],
                "p_value": p_value[indices],
                "n_samples": n,
                "nonzero_samples_full_cohort": [
                    int(prevalence_counts[eligible[index]]) for index in indices
                ],
                "prevalence_full_cohort": [
                    float(prevalence_counts[eligible[index]] / len(raw_values))
                    for index in indices
                ],
            }
        )
        results.append(block)

    if not results:
        raise RuntimeError("No estimable correlations")
    correlations = pd.concat(results, ignore_index=True)
    correlations["direction"] = np.where(correlations["spearman_r"] > 0, "positive", "negative")
    n_tests = len(correlations)
    correlations["p_bonferroni"] = np.minimum(correlations["p_value"] * n_tests, 1.0)
    correlations["p_fdr"] = multipletests(correlations["p_value"].to_numpy(), method="fdr_bh")[1]
    correlations["sig_p001"] = correlations["p_value"] < 0.001
    correlations["sig_fdr05"] = correlations["p_fdr"] < 0.05
    correlations["sig_bonferroni05"] = correlations["p_bonferroni"] < 0.05
    correlations["filter_criterion"] = (
        f"raw Pfam count >0 in >= {prevalence_minimum}/126 genomes "
        f"({PREVALENCE_FRACTION:.0%}) before environmental testing"
    )
    correlations["n_total_tests"] = n_tests
    correlations = correlations.sort_values(
        ["p_value", "env_var", "pfam"], kind="mergesort"
    ).reset_index(drop=True)

    # Directly verify vectorized values against scipy.stats.spearmanr for traced pairs.
    verification_indices = sorted({0, len(correlations) // 2, len(correlations) - 1})
    for index in verification_indices:
        row = correlations.iloc[index]
        y_series = pd.to_numeric(env[row["env_var"]], errors="coerce")
        x_series = pd.to_numeric(raw[row["pfam"]], errors="raise")
        valid = y_series.notna() & x_series.notna()
        expected_r, expected_p = spearmanr(x_series[valid], y_series[valid])
        if abs(expected_r - row["spearman_r"]) > 1e-12:
            raise RuntimeError(f"Vectorized rho mismatch at row {index}")
        if abs(expected_p - row["p_value"]) > 1e-12:
            raise RuntimeError(f"Vectorized p mismatch at row {index}")
    return correlations, prevalence_minimum, eligible


def write_deterministic_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    with path.open("wb") as binary:
        with gzip.GzipFile(filename="", mode="wb", fileobj=binary, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text_handle:
                frame.to_csv(text_handle, index=False, lineterminator="\n")


def summarize(correlations: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for variable in [*ENVIRONMENTAL_COLUMNS, "TOTAL"]:
        subset = correlations if variable == "TOTAL" else correlations[correlations["env_var"] == variable]
        if subset.empty:
            continue
        top = subset.sort_values("p_value", kind="mergesort").iloc[0]
        rows.append(
            {
                "variable": variable,
                "n_samples_min": int(subset["n_samples"].min()),
                "n_samples_max": int(subset["n_samples"].max()),
                "total_tests": len(subset),
                "p_lt_0_001": int(subset["sig_p001"].sum()),
                "fdr_lt_0_05": int(subset["sig_fdr05"].sum()),
                "bonferroni_lt_0_05": int(subset["sig_bonferroni05"].sum()),
                "unique_fdr_pfams": int(subset.loc[subset["sig_fdr05"], "pfam"].nunique()),
                "top_pfam": top["pfam"],
                "top_spearman_r": float(top["spearman_r"]),
                "top_p_value": float(top["p_value"]),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc)
    env, raw, env_manifest = load_inputs()
    correlations, prevalence_minimum, eligible = compute(env, raw)
    summary = summarize(correlations)
    significant = correlations[correlations["sig_fdr05"]].copy()

    write_deterministic_gzip_csv(correlations, OUTPUT_ALL)
    significant.to_csv(OUTPUT_SIGNIFICANT, index=False, lineterminator="\n")
    summary.to_csv(OUTPUT_SUMMARY, index=False, lineterminator="\n")

    # Re-read every materialized output.
    reread_all = pd.read_csv(OUTPUT_ALL, compression="gzip", low_memory=False)
    reread_sig = pd.read_csv(OUTPUT_SIGNIFICANT, low_memory=False)
    reread_summary = pd.read_csv(OUTPUT_SUMMARY, low_memory=False)
    if len(reread_all) != len(correlations):
        raise RuntimeError("All-results row count changed on materialization")
    if len(reread_sig) != int(correlations["sig_fdr05"].sum()):
        raise RuntimeError("FDR-results row count changed on materialization")
    if reread_summary.iloc[-1]["variable"] != "TOTAL":
        raise RuntimeError("Summary lacks TOTAL row")

    completed = datetime.now(timezone.utc)
    total = summary[summary["variable"] == "TOTAL"].iloc[0].to_dict()
    manifest = {
        "purpose": "bounded exact-genome-ID correction/validation of original raw-count GEE correlations",
        "not_a_new_discovery_analysis": True,
        "generated_by": str(Path(__file__).relative_to(ROOT)),
        "generator_sha256": sha256(Path(__file__)),
        "started_utc": started.isoformat(),
        "completed_utc": completed.isoformat(),
        "runtime_seconds": (completed - started).total_seconds(),
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "statsmodels": statsmodels.__version__,
        },
        "inputs": {
            str(ENV.relative_to(ROOT)): sha256(ENV),
            str(ENV_MANIFEST.relative_to(ROOT)): sha256(ENV_MANIFEST),
            str(RAW.relative_to(ROOT)): sha256(RAW),
        },
        "specification": {
            "join_key": "exact canonical genome_id",
            "cohort_rows": len(env),
            "prevalence_fraction": PREVALENCE_FRACTION,
            "prevalence_minimum_nonzero_genomes": prevalence_minimum,
            "eligible_pfams_before_variable_missingness": len(eligible),
            "minimum_valid_pairs": MIN_VALID_PAIRS,
            "effect": "Spearman rank correlation of raw reconstructed HMM count and GEE variable",
            "multiplicity": "Benjamini-Hochberg and Bonferroni across the full materialized raw-count GEE test family",
            "environmental_extraction_manifest_output_sha256": env_manifest["output"]["sha256"],
        },
        "results": {
            "total_tests": int(total["total_tests"]),
            "p_lt_0_001": int(total["p_lt_0_001"]),
            "fdr_lt_0_05": int(total["fdr_lt_0_05"]),
            "bonferroni_lt_0_05": int(total["bonferroni_lt_0_05"]),
            "unique_fdr_pfams": int(total["unique_fdr_pfams"]),
        },
        "outputs": {
            str(OUTPUT_ALL.relative_to(ROOT)): sha256(OUTPUT_ALL),
            str(OUTPUT_SIGNIFICANT.relative_to(ROOT)): sha256(OUTPUT_SIGNIFICANT),
            str(OUTPUT_SUMMARY.relative_to(ROOT)): sha256(OUTPUT_SUMMARY),
        },
        "integrity_checks": {
            "exact_126_id_sets_equal": True,
            "no_species_name_join": True,
            "no_imputation": True,
            "no_synthetic_values": True,
            "prevalence_filter_applied_before_environmental_testing": True,
            "three_vectorized_results_verified_against_scipy_spearmanr": True,
            "outputs_reopened_and_row_counts_verified": True,
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"summary": summary.to_dict(orient="records"), "manifest": str(MANIFEST)}, indent=2))


if __name__ == "__main__":
    main()

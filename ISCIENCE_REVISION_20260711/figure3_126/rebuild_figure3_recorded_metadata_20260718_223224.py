#!/usr/bin/env python3
"""Recompute Figure 3 from the authenticated AEF-126 cohort and raw Pfam counts.

The analysis uses only exact canonical genome identifiers, the reconciled
metadata manifest, and Pfam counts reconstructed from per-genome HMM-search
outputs.  The recorded environment field is an ordered category.  Its
prespecified ordinal coding

    Freshwater = 0, Brackish = 1, Marine = 2

is used solely to represent the recorded order in a Spearman rank test.  No
unrecorded numeric environment values are assigned or inferred.

All four pooled genome-level discovery variables form one multiplicity family:
temperature, latitude, longitude, and recorded environment category
(Freshwater < Brackish < Marine).  Benjamini--Hochberg and Bonferroni
corrections are applied globally across every estimable strict-Pfam test.

No synthetic, simulated, imputed, subsampled, or hardcoded result values are
used.  Functional descriptions are deliberately omitted: plotted rows are
identified by exact Pfam accession only.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import math
import platform
import re
import sys
from collections import OrderedDict
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import statsmodels
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Patch
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from statsmodels.stats.multitest import multipletests


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]

DEFAULT_MATRIX = (
    PROJECT_ROOT
    / "ISCIENCE_REVISION_20260711"
    / "analysis_stats"
    / "reconstructed_raw_pfam_counts_20260711_131706.csv.gz"
)
DEFAULT_MANIFEST = (
    PROJECT_ROOT
    / "ISCIENCE_REVISION_20260711"
    / "integrity"
    / "reconciled_analysis_manifest_20260711_110650.csv"
)
DEFAULT_AEF = (
    PROJECT_ROOT
    / "AlphaEarth"
    / "CSV"
    / "alphaearth_embeddings_20251019_122918.csv"
)
DEFAULT_OUTPUT_DIR = SCRIPT_PATH.parent

EXPECTED_INPUT_SHA256 = {
    "matrix": "683228342a90d2ecf2897930cd7c147f23973de0444155bfc162d50b09dd22bb",
    "manifest": "5880f930192d4cbb11a7563825c853e458fb1690b01c9cf7cc83323e7541bd67",
    "aef": "e0b05e727aec4a5b45565c9026de21c46561a18c06786f235eae371771a4cf87",
}
EXPECTED_COHORT_SIZE = 126
EXPECTED_PHYLA = OrderedDict(
    [("Rhodophyta", 70), ("Ochrophyta", 43), ("Chlorophyta", 13)]
)
PFAM_PATTERN = re.compile(r"^PF\d{5}$")

RECORDED_ENVIRONMENT_VARIABLE = (
    "recorded environment category (Freshwater < Brackish < Marine)"
)
RECORDED_ENVIRONMENT_ORDER = OrderedDict(
    [("Freshwater", 0), ("Brackish", 1), ("Marine", 2)]
)
VARIABLES = OrderedDict(
    [
        ("temperature", "Temperature (C)"),
        ("latitude", "Latitude"),
        ("longitude", "Longitude"),
        (RECORDED_ENVIRONMENT_VARIABLE, RECORDED_ENVIRONMENT_VARIABLE),
    ]
)
P001_THRESHOLD = 0.001
FDR_THRESHOLD = 0.05
BONFERRONI_ALPHA = 0.05
TOP_CATEGORY_PFAMS = 10

ENVIRONMENT_COLORS = {
    "Freshwater": "#0072B2",
    "Brackish": "#E69F00",
    "Marine": "#009E73",
}
PHYLUM_COLORS = {
    "Rhodophyta": "#CC79A7",
    "Ochrophyta": "#E69F00",
    "Chlorophyta": "#009E73",
}
NEGATIVE_COLOR = "#0072B2"
POSITIVE_COLOR = "#D55E00"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--aef", type=Path, default=DEFAULT_AEF)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--run-id",
        default=None,
        help="Optional YYYYMMDD_HHMMSS output tag; defaults to local current time.",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Run every full-data computation and render in memory without writing outputs.",
    )
    return parser.parse_args()


def run_stamp() -> str:
    return datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def project_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(PROJECT_ROOT.resolve()))
    except ValueError:
        return str(path)


def file_record(path: Path) -> dict[str, Any]:
    return {
        "path": project_path(path),
        "bytes": int(path.stat().st_size),
        "sha256": sha256_file(path),
    }


def parse_boolean(series: pd.Series, name: str) -> pd.Series:
    true_values = {"true", "1", "yes"}
    false_values = {"false", "0", "no"}
    normalized = series.astype(str).str.strip().str.lower()
    unknown = ~normalized.isin(true_values | false_values)
    if unknown.any():
        raise ValueError(f"{name} contains invalid booleans: {series[unknown].head().tolist()}")
    return normalized.isin(true_values)


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica"],
            "font.size": 6,
            "axes.labelsize": 6,
            "axes.titlesize": 6,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "axes.linewidth": 0.25,
            "xtick.major.width": 0.25,
            "ytick.major.width": 0.25,
            "xtick.major.size": 2,
            "ytick.major.size": 2,
            "axes.labelpad": 1,
            "xtick.major.pad": 1,
            "ytick.major.pad": 1,
            "axes.unicode_minus": False,
        }
    )


def style_axis(ax: mpl.axes.Axes, grid: bool = False) -> None:
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    ax.tick_params(which="both", width=0.25, labelsize=6)
    if grid:
        ax.grid(axis="y", linewidth=0.25, alpha=0.20, color="#777777")
        ax.set_axisbelow(True)


def load_and_validate(
    matrix_path: Path, manifest_path: Path, aef_path: Path
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict[str, Any]]:
    paths = {"matrix": matrix_path, "manifest": manifest_path, "aef": aef_path}
    for key, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        observed = sha256_file(path)
        expected = EXPECTED_INPUT_SHA256[key]
        if observed != expected:
            raise RuntimeError(
                f"Pinned {key} input hash mismatch: observed {observed}, expected {expected}"
            )

    manifest = pd.read_csv(manifest_path, low_memory=False)
    required_manifest = [
        "Genome",
        "Species",
        "Phylum",
        "Temperature (°C)",
        "Environment",
        "DD latitude",
        "DD longitude",
        "raw_pfam_hit_total",
        "safe_for_aef_pfam_analysis",
    ]
    missing_manifest = [column for column in required_manifest if column not in manifest]
    if missing_manifest:
        raise ValueError(f"Manifest is missing required fields: {missing_manifest}")

    safe = parse_boolean(
        manifest["safe_for_aef_pfam_analysis"], "safe_for_aef_pfam_analysis"
    )
    cohort = manifest.loc[safe, required_manifest].copy()
    if len(cohort) != EXPECTED_COHORT_SIZE or cohort["Genome"].nunique() != EXPECTED_COHORT_SIZE:
        raise ValueError(
            f"Expected {EXPECTED_COHORT_SIZE} unique safe AEF rows; observed {len(cohort)} rows "
            f"and {cohort['Genome'].nunique()} IDs"
        )
    if cohort[required_manifest].isna().any().any():
        bad = cohort.columns[cohort.isna().any()].tolist()
        raise ValueError(f"Selected cohort contains missing required values: {bad}")

    aef_ids = pd.read_csv(aef_path, usecols=["Genome"])["Genome"].astype(str)
    if len(aef_ids) != EXPECTED_COHORT_SIZE or aef_ids.nunique() != EXPECTED_COHORT_SIZE:
        raise ValueError("Pinned AEF table is not 126 unique genome IDs")
    cohort_ids = set(cohort["Genome"].astype(str))
    if cohort_ids != set(aef_ids):
        raise ValueError(
            f"Manifest/AEF ID mismatch: manifest-only={sorted(cohort_ids-set(aef_ids))[:10]}, "
            f"AEF-only={sorted(set(aef_ids)-cohort_ids)[:10]}"
        )

    cohort["Genome"] = cohort["Genome"].astype(str)
    cohort = cohort.set_index("Genome").loc[aef_ids].reset_index()
    phylum_counts = cohort["Phylum"].value_counts().to_dict()
    if phylum_counts != dict(EXPECTED_PHYLA):
        raise ValueError(
            f"Unexpected phylum composition: {phylum_counts}; expected {dict(EXPECTED_PHYLA)}"
        )
    observed_environments = set(cohort["Environment"])
    if observed_environments != set(RECORDED_ENVIRONMENT_ORDER):
        raise ValueError(
            f"Unexpected recorded environment categories: {observed_environments}"
        )

    for column in ["Temperature (°C)", "DD latitude", "DD longitude", "raw_pfam_hit_total"]:
        numeric = pd.to_numeric(cohort[column], errors="coerce")
        if numeric.isna().any() or not np.isfinite(numeric.to_numpy()).all():
            raise ValueError(f"Cohort field is not complete finite numeric data: {column}")
        cohort[column] = numeric
    cohort[RECORDED_ENVIRONMENT_VARIABLE] = cohort["Environment"].map(
        RECORDED_ENVIRONMENT_ORDER
    )

    raw = pd.read_csv(matrix_path, low_memory=False)
    if "Genome" not in raw:
        raise ValueError("Raw matrix lacks Genome")
    if raw["Genome"].isna().any() or raw["Genome"].duplicated().any():
        raise ValueError("Raw matrix Genome identifiers must be complete and unique")
    raw["Genome"] = raw["Genome"].astype(str)
    if set(raw["Genome"]) != cohort_ids:
        raise ValueError(
            f"Raw matrix ID mismatch: raw-only={sorted(set(raw.Genome)-cohort_ids)[:10]}, "
            f"manifest-only={sorted(cohort_ids-set(raw.Genome))[:10]}"
        )
    raw = raw.set_index("Genome").loc[aef_ids]
    pfams = sorted(column for column in raw.columns if PFAM_PATTERN.fullmatch(str(column)))
    if not pfams:
        raise ValueError("Raw matrix contains no strict ^PF\\d{5}$ columns")
    prefixed_non_accessions = [
        column for column in raw.columns if str(column).startswith("PF") and column not in pfams
    ]
    if prefixed_non_accessions:
        raise ValueError(f"Non-accession PF-prefixed columns are prohibited: {prefixed_non_accessions}")

    counts = raw[pfams]
    if any(not pd.api.types.is_numeric_dtype(counts[column]) for column in pfams):
        raise ValueError("At least one strict Pfam count column is not numeric")
    values = counts.to_numpy(dtype=np.float64, copy=False)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("Raw Pfam counts must be complete, finite, and nonnegative")
    if not np.equal(values, np.floor(values)).all():
        raise ValueError("Raw Pfam counts must be integers")
    row_sums = counts.sum(axis=1).to_numpy(dtype=np.int64)
    expected_sums = cohort.set_index("Genome").loc[aef_ids, "raw_pfam_hit_total"].to_numpy(
        dtype=np.int64
    )
    if not np.array_equal(row_sums, expected_sums):
        raise ValueError("Strict-Pfam row sums differ from manifest raw_pfam_hit_total")
    if np.any(row_sums <= 0):
        raise ValueError("At least one authenticated genome has an all-zero Pfam vector")

    validation = {
        "cohort_rows": int(len(cohort)),
        "unique_genome_ids": int(cohort["Genome"].nunique()),
        "strict_pfam_columns": int(len(pfams)),
        "phylum_counts": {key: int(phylum_counts[key]) for key in EXPECTED_PHYLA},
        "recorded_environment_counts": {
            key: int((cohort["Environment"] == key).sum())
            for key in RECORDED_ENVIRONMENT_ORDER
        },
        "aef_id_set_exact": True,
        "raw_matrix_id_set_exact": True,
        "strict_pfam_row_sums_match_manifest": True,
        "all_raw_pfam_rows_nonzero": True,
    }
    return cohort, counts, pfams, validation


def vector_spearman(
    ranked_counts_centered: np.ndarray,
    count_sum_squares: np.ndarray,
    outcome: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ranked_outcome = stats.rankdata(outcome, method="average")
    ranked_outcome = ranked_outcome - ranked_outcome.mean()
    outcome_ss = float(np.dot(ranked_outcome, ranked_outcome))
    estimable = (count_sum_squares > 0) & (outcome_ss > 0)
    r = np.full(len(count_sum_squares), np.nan, dtype=float)
    r[estimable] = (
        ranked_counts_centered[:, estimable].T @ ranked_outcome
    ) / np.sqrt(count_sum_squares[estimable] * outcome_ss)
    r[estimable] = np.clip(r[estimable], -1.0, 1.0)
    p = np.full(len(count_sum_squares), np.nan, dtype=float)
    n = len(outcome)
    interior = estimable & (np.abs(r) < 1.0)
    statistic = r[interior] * np.sqrt(
        (n - 2) / ((1.0 + r[interior]) * (1.0 - r[interior]))
    )
    p[interior] = 2.0 * stats.t.sf(np.abs(statistic), df=n - 2)
    p[estimable & (np.abs(r) == 1.0)] = 0.0
    return r, p, estimable


def compute_correlations(
    cohort: pd.DataFrame, counts: pd.DataFrame, pfams: list[str]
) -> tuple[pd.DataFrame, dict[str, Any]]:
    count_values = counts.to_numpy(dtype=np.float64, copy=True)
    ranked_counts = stats.rankdata(count_values, axis=0, method="average")
    ranked_counts -= ranked_counts.mean(axis=0, keepdims=True)
    count_ss = np.einsum("ij,ij->j", ranked_counts, ranked_counts)

    outcomes = OrderedDict(
        [
            ("temperature", cohort["Temperature (°C)"].to_numpy(dtype=float)),
            ("latitude", cohort["DD latitude"].to_numpy(dtype=float)),
            ("longitude", cohort["DD longitude"].to_numpy(dtype=float)),
            (
                RECORDED_ENVIRONMENT_VARIABLE,
                cohort[RECORDED_ENVIRONMENT_VARIABLE].to_numpy(dtype=float),
            ),
        ]
    )
    frames: list[pd.DataFrame] = []
    max_r_error = 0.0
    max_p_error = 0.0
    checks = 0
    for variable, outcome in outcomes.items():
        r, p, estimable = vector_spearman(ranked_counts, count_ss, outcome)
        indices = np.flatnonzero(estimable)
        coding = (
            "prespecified ordinal 0/1/2 solely representing Freshwater < Brackish < Marine"
            if variable == RECORDED_ENVIRONMENT_VARIABLE
            else "continuous recorded numeric value"
        )
        frames.append(
            pd.DataFrame(
                {
                    "pfam": [pfams[index] for index in indices],
                    "variable": variable,
                    "spearman_r": r[indices],
                    "p_value": p[indices],
                    "n_samples": len(outcome),
                    "variable_coding": coding,
                }
            )
        )

        check_indices = sorted({0, len(indices) // 2, len(indices) - 1})
        for relative_index in check_indices:
            column_index = indices[relative_index]
            reference_r, reference_p = stats.spearmanr(
                count_values[:, column_index], outcome
            )
            max_r_error = max(max_r_error, abs(reference_r - r[column_index]))
            max_p_error = max(max_p_error, abs(reference_p - p[column_index]))
            checks += 1

    result = pd.concat(frames, ignore_index=True)
    if len(result) != len(pfams) * len(outcomes):
        raise ValueError(
            f"Expected every strict Pfam to be estimable for four variables; observed {len(result)} "
            f"rather than {len(pfams) * len(outcomes)}"
        )
    if max_r_error > 1e-12 or max_p_error > 1e-12:
        raise RuntimeError(
            "Vectorized Spearman results differ from scipy.stats.spearmanr: "
            f"max r error={max_r_error}, max p error={max_p_error}"
        )

    n_tests = len(result)
    result["p_bonferroni"] = np.minimum(result["p_value"] * n_tests, 1.0)
    result["p_fdr_bh"] = multipletests(
        result["p_value"].to_numpy(), method="fdr_bh"
    )[1]
    result["direction"] = np.select(
        [result["spearman_r"] > 0, result["spearman_r"] < 0],
        ["positive", "negative"],
        default="zero",
    )
    result["sig_p001"] = result["p_value"] < P001_THRESHOLD
    result["sig_fdr05"] = result["p_fdr_bh"] < FDR_THRESHOLD
    result["sig_bonferroni05"] = result["p_bonferroni"] < BONFERRONI_ALPHA
    result["global_test_family_size"] = n_tests
    order = {variable: index for index, variable in enumerate(outcomes)}
    result["_variable_order"] = result["variable"].map(order)
    result = result.sort_values(
        ["_variable_order", "p_value", "pfam"], kind="mergesort"
    ).drop(columns="_variable_order")
    result = result.reset_index(drop=True)

    crosscheck = {
        "comparisons": checks,
        "maximum_absolute_r_error": float(max_r_error),
        "maximum_absolute_p_error": float(max_p_error),
    }
    return result, crosscheck


def summarize(correlations: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for variable in VARIABLES:
        subset = correlations.loc[correlations["variable"].eq(variable)]
        top = subset.sort_values(["p_value", "pfam"], kind="mergesort").iloc[0]
        p001 = subset["sig_p001"]
        rows.append(
            {
                "variable": variable,
                "n_samples": int(subset["n_samples"].min()),
                "total_tests": int(len(subset)),
                "p_lt_0_001": int(p001.sum()),
                "positive_p_lt_0_001": int((p001 & subset["spearman_r"].gt(0)).sum()),
                "negative_p_lt_0_001": int((p001 & subset["spearman_r"].lt(0)).sum()),
                "bh_q_lt_0_05": int(subset["sig_fdr05"].sum()),
                "bonferroni_p_lt_0_05": int(subset["sig_bonferroni05"].sum()),
                "top_pfam": top["pfam"],
                "top_spearman_r": float(top["spearman_r"]),
                "top_p_value": float(top["p_value"]),
                "top_bh_q": float(top["p_fdr_bh"]),
            }
        )
    top = correlations.sort_values(["p_value", "variable", "pfam"], kind="mergesort").iloc[0]
    p001 = correlations["sig_p001"]
    rows.append(
        {
            "variable": "TOTAL_GLOBAL_FAMILY",
            "n_samples": int(correlations["n_samples"].min()),
            "total_tests": int(len(correlations)),
            "p_lt_0_001": int(p001.sum()),
            "positive_p_lt_0_001": int(
                (p001 & correlations["spearman_r"].gt(0)).sum()
            ),
            "negative_p_lt_0_001": int(
                (p001 & correlations["spearman_r"].lt(0)).sum()
            ),
            "bh_q_lt_0_05": int(correlations["sig_fdr05"].sum()),
            "bonferroni_p_lt_0_05": int(correlations["sig_bonferroni05"].sum()),
            "top_pfam": top["pfam"],
            "top_spearman_r": float(top["spearman_r"]),
            "top_p_value": float(top["p_value"]),
            "top_bh_q": float(top["p_fdr_bh"]),
        }
    )
    return pd.DataFrame(rows)


def prepare_panel_data(
    cohort: pd.DataFrame, counts: pd.DataFrame, correlations: pd.DataFrame
) -> dict[str, Any]:
    category_results = correlations.loc[
        correlations["variable"].eq(RECORDED_ENVIRONMENT_VARIABLE)
    ].sort_values(["p_value", "pfam"], kind="mergesort")
    top = category_results.head(TOP_CATEGORY_PFAMS).copy()
    top_pfams = top["pfam"].tolist()

    group_rows: list[dict[str, Any]] = []
    for pfam in top_pfams:
        stats_row = top.loc[top["pfam"].eq(pfam)].iloc[0]
        for environment in RECORDED_ENVIRONMENT_ORDER:
            mask = cohort["Environment"].eq(environment).to_numpy()
            raw_values = counts.loc[mask, pfam].to_numpy(dtype=float)
            mean_raw = float(raw_values.mean())
            group_rows.append(
                {
                    "pfam": pfam,
                    "recorded_environment_category": environment,
                    "ordinal_code": RECORDED_ENVIRONMENT_ORDER[environment],
                    "n_genomes": int(mask.sum()),
                    "mean_raw_count": mean_raw,
                    "log2_mean_raw_count_plus_1": float(np.log2(mean_raw + 1.0)),
                    "spearman_r": float(stats_row["spearman_r"]),
                    "p_value": float(stats_row["p_value"]),
                    "p_fdr_bh_global": float(stats_row["p_fdr_bh"]),
                    "p_bonferroni_global": float(stats_row["p_bonferroni"]),
                }
            )
    top_groups = pd.DataFrame(group_rows)

    supported = category_results.loc[category_results["sig_fdr05"]].copy()
    if supported.empty:
        raise RuntimeError(
            f"No BH-supported results for {RECORDED_ENVIRONMENT_VARIABLE}; panel C is undefined"
        )
    panel_pfams = supported["pfam"].tolist()
    raw_matrix = counts.loc[:, panel_pfams].T
    log_matrix = np.log2(raw_matrix.astype(float) + 1.0)
    row_std = log_matrix.std(axis=1, ddof=0)
    if row_std.le(0).any():
        raise RuntimeError("A panel-C Pfam profile is constant")
    clustering_matrix = log_matrix.sub(log_matrix.mean(axis=1), axis=0).div(
        row_std, axis=0
    )
    row_linkage = linkage(
        clustering_matrix.to_numpy(),
        method="average",
        metric="euclidean",
        optimal_ordering=True,
    )
    column_linkage = linkage(
        clustering_matrix.to_numpy().T,
        method="average",
        metric="euclidean",
        optimal_ordering=True,
    )
    row_order = leaves_list(row_linkage)
    column_order = leaves_list(column_linkage)
    ordered_pfams = [panel_pfams[index] for index in row_order]
    ordered_genomes = [counts.index[index] for index in column_order]
    ordered_raw = raw_matrix.loc[ordered_pfams, ordered_genomes]
    ordered_log = log_matrix.loc[ordered_pfams, ordered_genomes]
    ordered_stats = supported.set_index("pfam").loc[ordered_pfams].reset_index()
    ordered_metadata = cohort.set_index("Genome").loc[ordered_genomes].reset_index()
    ordered_metadata["clustered_sample_position"] = np.arange(len(ordered_metadata))
    ordered_stats["clustered_pfam_position"] = np.arange(len(ordered_stats))

    return {
        "top_results": top,
        "top_groups": top_groups,
        "supported": supported,
        "ordered_raw": ordered_raw,
        "ordered_log": ordered_log,
        "ordered_stats": ordered_stats,
        "ordered_metadata": ordered_metadata,
        "row_linkage": row_linkage,
        "column_linkage": column_linkage,
    }


def add_discrete_track(
    ax: mpl.axes.Axes,
    labels: list[str],
    categories: list[str],
    colors: dict[str, str],
    ylabel: str,
) -> None:
    codes = np.asarray([categories.index(label) for label in labels], dtype=float)[None, :]
    cmap = ListedColormap([colors[category] for category in categories])
    norm = BoundaryNorm(np.arange(-0.5, len(categories) + 0.5), cmap.N)
    ax.pcolormesh(
        np.arange(len(labels) + 1),
        np.arange(2),
        codes,
        cmap=cmap,
        norm=norm,
        shading="flat",
        edgecolors="none",
        linewidth=0,
    )
    ax.set_xlim(0, len(labels))
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel(ylabel, rotation=0, ha="right", va="center", labelpad=2)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)


def add_continuous_track(
    ax: mpl.axes.Axes, values: np.ndarray, ylabel: str, cmap: mpl.colors.Colormap
) -> None:
    minimum = float(np.min(values))
    maximum = float(np.max(values))
    if not maximum > minimum:
        raise RuntimeError(f"Continuous track is constant: {ylabel}")
    normalized = ((values - minimum) / (maximum - minimum))[None, :]
    ax.pcolormesh(
        np.arange(len(values) + 1),
        np.arange(2),
        normalized,
        cmap=cmap,
        vmin=0,
        vmax=1,
        shading="flat",
        edgecolors="none",
        linewidth=0,
    )
    ax.set_xlim(0, len(values))
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel(ylabel, rotation=0, ha="right", va="center", labelpad=2)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)


def build_figure(
    summary: pd.DataFrame,
    panel: dict[str, Any],
    family_size: int,
) -> mpl.figure.Figure:
    configure_matplotlib()
    abundance_cmap = LinearSegmentedColormap.from_list(
        "pfam_abundance",
        [
            (1.0, 1.0, 1.0),
            (0.95, 0.95, 0.95),
            (0.10, 0.10, 0.10),
            (0.40, 0.00, 0.00),
            (0.70, 0.00, 0.00),
            (0.90, 0.20, 0.00),
            (1.00, 0.70, 0.00),
        ],
    )
    continuous_cmap = LinearSegmentedColormap.from_list(
        "metadata_continuous",
        ["#132B43", "#1D6996", "#73AF48", "#FDE725"],
    )

    fig = plt.figure(figsize=(8.27, 10.8))
    fig.patch.set_alpha(0.0)
    outer = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.18, 4.25],
        width_ratios=[1.0, 1.34],
        left=0.075,
        right=0.975,
        bottom=0.055,
        top=0.975,
        hspace=0.24,
        wspace=0.23,
    )

    # Panel A: exact counts from the global family.
    ax_a = fig.add_subplot(outer[0, 0])
    rows = summary.loc[summary["variable"].isin(VARIABLES)].set_index("variable").loc[
        list(VARIABLES)
    ]
    x = np.arange(len(rows), dtype=float)
    negative = rows["negative_p_lt_0_001"].to_numpy(dtype=float)
    positive = rows["positive_p_lt_0_001"].to_numpy(dtype=float)
    bh = rows["bh_q_lt_0_05"].to_numpy(dtype=float)
    bonf = rows["bonferroni_p_lt_0_05"].to_numpy(dtype=float)
    width = 0.25
    ax_a.bar(
        x - width,
        negative,
        width=width,
        color=NEGATIVE_COLOR,
        edgecolor="none",
        label="Negative, p<0.001",
    )
    ax_a.bar(
        x - width,
        positive,
        bottom=negative,
        width=width,
        color=POSITIVE_COLOR,
        edgecolor="none",
        label="Positive, p<0.001",
    )
    ax_a.bar(
        x + 0.03,
        bh,
        width=width,
        color="#777777",
        edgecolor="none",
        label="Global BH q<0.05",
    )
    ax_a.bar(
        x + width + 0.06,
        bonf,
        width=width,
        color="#111111",
        edgecolor="none",
        label="Global Bonferroni p<0.05",
    )
    for index, total in enumerate(negative + positive):
        ax_a.text(index - width, total + 3, f"{int(total)}", ha="center", va="bottom")
    for index, value in enumerate(bh):
        ax_a.text(index + 0.03, value + 3, f"{int(value)}", ha="center", va="bottom")
    for index, value in enumerate(bonf):
        ax_a.text(index + width + 0.06, value + 3, f"{int(value)}", ha="center", va="bottom")
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(
        [
            "temperature",
            "latitude",
            "longitude",
            "recorded environment category\n(Freshwater < Brackish < Marine)",
        ],
        rotation=28,
        ha="right",
    )
    ax_a.set_ylabel("Number of Pfam-variable associations")
    ax_a.set_title(
        f"A  Pooled genome-level discovery family ({family_size:,} tests)",
        loc="left",
        fontweight="bold",
        pad=2,
    )
    ax_a.legend(loc="upper left", frameon=False, ncol=1)
    style_axis(ax_a, grid=True)

    # Panel B: top ten recorded-environment-category results.
    grid_b = outer[0, 1].subgridspec(
        1, 4, width_ratios=[3.2, 1.0, 1.0, 0.16], wspace=0.12
    )
    ax_b = fig.add_subplot(grid_b[0, 0])
    top_results = panel["top_results"].set_index("pfam")
    top_pfams = panel["top_results"]["pfam"].tolist()
    top_groups = panel["top_groups"]
    group_matrix = (
        top_groups.pivot(
            index="pfam",
            columns="recorded_environment_category",
            values="log2_mean_raw_count_plus_1",
        )
        .loc[top_pfams, list(RECORDED_ENVIRONMENT_ORDER)]
        .to_numpy(dtype=float)
    )
    mesh_b = ax_b.pcolormesh(
        np.arange(4),
        np.arange(len(top_pfams) + 1),
        group_matrix,
        cmap=abundance_cmap,
        shading="flat",
        edgecolors="#DDDDDD",
        linewidth=0.15,
    )
    ax_b.set_xlim(0, 3)
    ax_b.set_ylim(len(top_pfams), 0)
    ax_b.set_xticks(np.arange(3) + 0.5)
    category_counts = panel["ordered_metadata"]["Environment"].value_counts()
    ax_b.set_xticklabels(
        [f"{category}\n(n={int(category_counts[category])})" for category in RECORDED_ENVIRONMENT_ORDER],
        rotation=30,
        ha="right",
    )
    ax_b.set_yticks(np.arange(len(top_pfams)) + 0.5)
    ax_b.set_yticklabels(top_pfams)
    ax_b.set_title(
        "B  Ten strongest recorded environment category\n(Freshwater < Brackish < Marine) associations",
        loc="left",
        fontweight="bold",
        pad=2,
    )
    style_axis(ax_b)

    ax_b_r = fig.add_subplot(grid_b[0, 1])
    r_values = top_results.loc[top_pfams, "spearman_r"].to_numpy(dtype=float)
    y = np.arange(len(top_pfams)) + 0.5
    ax_b_r.barh(
        y,
        r_values,
        height=0.66,
        color=[POSITIVE_COLOR if value > 0 else NEGATIVE_COLOR for value in r_values],
        edgecolor="none",
    )
    ax_b_r.axvline(0, color="black", linewidth=0.25)
    ax_b_r.set_ylim(len(top_pfams), 0)
    limit = max(0.1, float(np.ceil(np.max(np.abs(r_values)) * 10) / 10))
    ax_b_r.set_xlim(-limit, limit)
    ax_b_r.set_yticks([])
    ax_b_r.set_xlabel("Spearman r")
    ax_b_r.set_title("Effect", fontweight="bold", pad=2)
    style_axis(ax_b_r)

    ax_b_q = fig.add_subplot(grid_b[0, 2])
    q_values = top_results.loc[top_pfams, "p_fdr_bh"].to_numpy(dtype=float)
    neg_log_q = -np.log10(q_values)
    ax_b_q.barh(y, neg_log_q, height=0.66, color="#777777", edgecolor="none")
    ax_b_q.axvline(-math.log10(FDR_THRESHOLD), color="black", linewidth=0.25, linestyle="--")
    ax_b_q.set_ylim(len(top_pfams), 0)
    ax_b_q.set_yticks([])
    ax_b_q.set_xlabel("-log10(global BH q)")
    ax_b_q.set_title("Multiplicity", fontweight="bold", pad=2)
    style_axis(ax_b_q)

    cax_b = fig.add_subplot(grid_b[0, 3])
    colorbar_b = fig.colorbar(mesh_b, cax=cax_b)
    colorbar_b.set_label("Log2(mean raw count + 1)")
    colorbar_b.outline.set_linewidth(0.25)
    colorbar_b.ax.tick_params(width=0.25, labelsize=6)
    colorbar_b.solids.set_rasterized(False)

    # Panel C: all BH-supported recorded-environment-category results.
    grid_c = outer[1, :].subgridspec(
        8,
        6,
        height_ratios=[0.74, 0.14, 0.14, 0.14, 0.14, 0.14, 5.05, 0.45],
        width_ratios=[0.92, 5.02, 0.54, 0.76, 0.88, 0.18],
        hspace=0.025,
        wspace=0.045,
    )
    ax_col_dendro = fig.add_subplot(grid_c[0, 1])
    dendrogram(
        panel["column_linkage"],
        ax=ax_col_dendro,
        no_labels=True,
        color_threshold=0,
        above_threshold_color="black",
        link_color_func=lambda _: "black",
    )
    ax_col_dendro.set_xticks([])
    ax_col_dendro.set_yticks([])
    for spine in ax_col_dendro.spines.values():
        spine.set_visible(False)
    ax_col_dendro.set_title(
        f"C  All {len(panel['ordered_stats'])} global-BH-supported {RECORDED_ENVIRONMENT_VARIABLE} Pfams",
        loc="left",
        fontweight="bold",
        pad=2,
    )

    metadata = panel["ordered_metadata"]
    environment_ax = fig.add_subplot(grid_c[1, 1])
    add_discrete_track(
        environment_ax,
        metadata["Environment"].tolist(),
        list(RECORDED_ENVIRONMENT_ORDER),
        ENVIRONMENT_COLORS,
        "recorded environment category\n(Freshwater < Brackish < Marine)",
    )
    phylum_ax = fig.add_subplot(grid_c[2, 1])
    add_discrete_track(
        phylum_ax,
        metadata["Phylum"].tolist(),
        list(EXPECTED_PHYLA),
        PHYLUM_COLORS,
        "phylum",
    )
    temperature_ax = fig.add_subplot(grid_c[3, 1])
    add_continuous_track(
        temperature_ax,
        metadata["Temperature (°C)"].to_numpy(dtype=float),
        f"temperature ({metadata['Temperature (°C)'].min():.1f}-{metadata['Temperature (°C)'].max():.1f} C)",
        continuous_cmap,
    )
    latitude_ax = fig.add_subplot(grid_c[4, 1])
    add_continuous_track(
        latitude_ax,
        metadata["DD latitude"].to_numpy(dtype=float),
        f"latitude ({metadata['DD latitude'].min():.1f}-{metadata['DD latitude'].max():.1f})",
        continuous_cmap,
    )
    longitude_ax = fig.add_subplot(grid_c[5, 1])
    add_continuous_track(
        longitude_ax,
        metadata["DD longitude"].to_numpy(dtype=float),
        f"longitude ({metadata['DD longitude'].min():.1f}-{metadata['DD longitude'].max():.1f})",
        continuous_cmap,
    )

    ax_row_dendro = fig.add_subplot(grid_c[6, 0])
    dendrogram(
        panel["row_linkage"],
        ax=ax_row_dendro,
        orientation="left",
        no_labels=True,
        color_threshold=0,
        above_threshold_color="black",
        link_color_func=lambda _: "black",
    )
    ax_row_dendro.invert_yaxis()
    ax_row_dendro.set_xticks([])
    ax_row_dendro.set_yticks([])
    for spine in ax_row_dendro.spines.values():
        spine.set_visible(False)

    ax_heatmap = fig.add_subplot(grid_c[6, 1])
    heatmap_values = panel["ordered_log"].to_numpy(dtype=float)
    mesh_c = ax_heatmap.pcolormesh(
        np.arange(heatmap_values.shape[1] + 1),
        np.arange(heatmap_values.shape[0] + 1),
        heatmap_values,
        cmap=abundance_cmap,
        shading="flat",
        edgecolors="none",
        linewidth=0,
    )
    ax_heatmap.set_xlim(0, heatmap_values.shape[1])
    ax_heatmap.set_ylim(heatmap_values.shape[0], 0)
    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks([])
    ax_heatmap.set_xlabel("126 authenticated genomes (hierarchically ordered)")
    style_axis(ax_heatmap)

    ax_c_labels = fig.add_subplot(grid_c[6, 2])
    ax_c_labels.set_xlim(0, 1)
    ax_c_labels.set_ylim(heatmap_values.shape[0], 0)
    ax_c_labels.axis("off")
    for index, pfam in enumerate(panel["ordered_stats"]["pfam"]):
        ax_c_labels.text(0.02, index + 0.5, pfam, ha="left", va="center")

    ax_c_r = fig.add_subplot(grid_c[6, 3])
    c_r = panel["ordered_stats"]["spearman_r"].to_numpy(dtype=float)
    c_y = np.arange(len(c_r)) + 0.5
    ax_c_r.barh(
        c_y,
        c_r,
        height=0.75,
        color=[POSITIVE_COLOR if value > 0 else NEGATIVE_COLOR for value in c_r],
        edgecolor="none",
    )
    ax_c_r.axvline(0, color="black", linewidth=0.25)
    ax_c_r.set_ylim(len(c_r), 0)
    c_limit = max(0.1, float(np.ceil(np.max(np.abs(c_r)) * 10) / 10))
    ax_c_r.set_xlim(-c_limit, c_limit)
    ax_c_r.set_yticks([])
    ax_c_r.set_xlabel("Spearman r")
    style_axis(ax_c_r)

    ax_c_q = fig.add_subplot(grid_c[6, 4])
    c_q = panel["ordered_stats"]["p_fdr_bh"].to_numpy(dtype=float)
    c_neg_log_q = -np.log10(c_q)
    ax_c_q.barh(c_y, c_neg_log_q, height=0.75, color="#777777", edgecolor="none")
    ax_c_q.axvline(-math.log10(FDR_THRESHOLD), color="black", linewidth=0.25, linestyle="--")
    ax_c_q.set_ylim(len(c_q), 0)
    ax_c_q.set_yticks([])
    ax_c_q.set_xlabel("-log10(global BH q)")
    style_axis(ax_c_q)

    cax_c = fig.add_subplot(grid_c[6, 5])
    colorbar_c = fig.colorbar(mesh_c, cax=cax_c)
    colorbar_c.set_label("Log2(raw Pfam count + 1)")
    colorbar_c.outline.set_linewidth(0.25)
    colorbar_c.ax.tick_params(width=0.25, labelsize=6)
    colorbar_c.solids.set_rasterized(False)

    legend_ax = fig.add_subplot(grid_c[7, :])
    legend_ax.axis("off")
    handles = [
        Patch(facecolor=ENVIRONMENT_COLORS[key], edgecolor="none", label=key)
        for key in RECORDED_ENVIRONMENT_ORDER
    ] + [
        Patch(facecolor=PHYLUM_COLORS[key], edgecolor="none", label=key)
        for key in EXPECTED_PHYLA
    ]
    legend_ax.legend(
        handles=handles,
        loc="center",
        ncol=len(handles),
        frameon=False,
        columnspacing=1.2,
        handlelength=1.4,
        handletextpad=0.35,
    )
    legend_ax.text(
        0.0,
        0.02,
        "Ordinal coding 0/1/2 solely represents the prespecified recorded order; all plotted row labels are exact Pfam accessions.",
        transform=legend_ax.transAxes,
        ha="left",
        va="bottom",
    )

    fig.canvas.draw()
    for text in fig.findobj(match=mpl.text.Text):
        if text.get_text() and abs(float(text.get_fontsize()) - 6.0) > 1e-9:
            raise RuntimeError(
                f"Figure protocol violation: text {text.get_text()!r} is {text.get_fontsize()} pt"
            )
    return fig


def write_deterministic_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    with path.open("wb") as binary:
        with gzip.GzipFile(filename="", mode="wb", fileobj=binary, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text_handle:
                frame.to_csv(text_handle, index=False, lineterminator="\n")


def create_output_paths(output_dir: Path, run_id: str) -> dict[str, Path]:
    return {
        "figure_pdf": output_dir / f"Figure3_recorded_metadata_126_{run_id}.pdf",
        "figure_svg": output_dir / f"Figure3_recorded_metadata_126_{run_id}.svg",
        "all_correlations": output_dir
        / f"Figure3_recorded_metadata_all_correlations_{run_id}.csv.gz",
        "significant_correlations": output_dir
        / f"Figure3_recorded_metadata_BH_q05_{run_id}.csv",
        "summary": output_dir / f"Figure3_recorded_metadata_summary_{run_id}.csv",
        "panel_b": output_dir
        / f"Figure3_recorded_environment_category_top10_group_means_{run_id}.csv",
        "panel_c_counts": output_dir
        / f"Figure3_recorded_environment_category_BH_heatmap_raw_counts_{run_id}.csv.gz",
        "panel_c_pfams": output_dir
        / f"Figure3_recorded_environment_category_BH_pfam_order_{run_id}.csv",
        "panel_c_samples": output_dir
        / f"Figure3_recorded_environment_category_sample_order_{run_id}.csv",
        "provenance": output_dir / f"Figure3_recorded_metadata_provenance_{run_id}.json",
    }


def write_outputs(
    paths: dict[str, Path],
    figure: mpl.figure.Figure,
    correlations: pd.DataFrame,
    summary: pd.DataFrame,
    panel: dict[str, Any],
    validation: dict[str, Any],
    crosscheck: dict[str, Any],
    args: argparse.Namespace,
    run_id: str,
) -> None:
    for path in paths.values():
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    paths["figure_pdf"].parent.mkdir(parents=True, exist_ok=True)

    figure.savefig(
        paths["figure_pdf"], format="pdf", bbox_inches="tight", transparent=True, edgecolor="none"
    )
    figure.savefig(
        paths["figure_svg"], format="svg", bbox_inches="tight", transparent=True, edgecolor="none"
    )
    write_deterministic_gzip_csv(correlations, paths["all_correlations"])
    correlations.loc[correlations["sig_fdr05"]].to_csv(
        paths["significant_correlations"], index=False, lineterminator="\n"
    )
    summary.to_csv(paths["summary"], index=False, lineterminator="\n")
    panel["top_groups"].to_csv(paths["panel_b"], index=False, lineterminator="\n")

    ordered_counts = panel["ordered_raw"].copy()
    ordered_counts.insert(0, "pfam", ordered_counts.index)
    write_deterministic_gzip_csv(ordered_counts.reset_index(drop=True), paths["panel_c_counts"])
    panel["ordered_stats"].to_csv(paths["panel_c_pfams"], index=False, lineterminator="\n")
    sample_columns = [
        "clustered_sample_position",
        "Genome",
        "Species",
        "Phylum",
        "Temperature (°C)",
        "Environment",
        RECORDED_ENVIRONMENT_VARIABLE,
        "DD latitude",
        "DD longitude",
        "raw_pfam_hit_total",
    ]
    panel["ordered_metadata"][sample_columns].to_csv(
        paths["panel_c_samples"], index=False, lineterminator="\n"
    )

    # Reopen every tabular output and validate computed row counts.
    reread_all = pd.read_csv(paths["all_correlations"], compression="gzip", low_memory=False)
    reread_sig = pd.read_csv(paths["significant_correlations"], low_memory=False)
    reread_summary = pd.read_csv(paths["summary"], low_memory=False)
    reread_panel_b = pd.read_csv(paths["panel_b"], low_memory=False)
    reread_panel_c = pd.read_csv(paths["panel_c_counts"], compression="gzip", low_memory=False)
    reread_pfams = pd.read_csv(paths["panel_c_pfams"], low_memory=False)
    reread_samples = pd.read_csv(paths["panel_c_samples"], low_memory=False)
    if len(reread_all) != len(correlations):
        raise RuntimeError("All-correlation output row count changed on materialization")
    if len(reread_sig) != int(correlations["sig_fdr05"].sum()):
        raise RuntimeError("BH output row count changed on materialization")
    if reread_summary.iloc[-1]["variable"] != "TOTAL_GLOBAL_FAMILY":
        raise RuntimeError("Summary output lacks TOTAL_GLOBAL_FAMILY")
    if len(reread_panel_b) != TOP_CATEGORY_PFAMS * len(RECORDED_ENVIRONMENT_ORDER):
        raise RuntimeError("Panel-B group output has an unexpected row count")
    if len(reread_panel_c) != len(panel["ordered_stats"]):
        raise RuntimeError("Panel-C matrix output has an unexpected row count")
    if len(reread_pfams) != len(panel["ordered_stats"]):
        raise RuntimeError("Panel-C Pfam-order output has an unexpected row count")
    if len(reread_samples) != EXPECTED_COHORT_SIZE:
        raise RuntimeError("Panel-C sample-order output has an unexpected row count")

    total = summary.loc[summary["variable"].eq("TOTAL_GLOBAL_FAMILY")].iloc[0]
    output_records = {
        key: file_record(path) for key, path in paths.items() if key != "provenance"
    }
    provenance = {
        "run_id": run_id,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "purpose": "corrected Figure 3 analysis on authenticated AEF-126 IDs and reconstructed raw Pfam counts",
        "generator": file_record(SCRIPT_PATH),
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "statsmodels": statsmodels.__version__,
            "matplotlib": mpl.__version__,
        },
        "inputs": {
            "matrix": file_record(args.matrix),
            "manifest": file_record(args.manifest),
            "aef_id_table": file_record(args.aef),
        },
        "validation": validation,
        "statistical_specification": {
            "cohort": "safe_for_aef_pfam_analysis rows; exact ID equality with the pinned AEF table",
            "pfam_columns": "strict regular expression ^PF\\d{5}$",
            "effect": "two-sided pooled genome-level Spearman rank correlation of raw Pfam count and recorded variable",
            "variables": list(VARIABLES),
            "recorded_environment_category_coding": {
                "mapping": dict(RECORDED_ENVIRONMENT_ORDER),
                "purpose": "prespecified ordinal 0/1/2 solely representing Freshwater < Brackish < Marine",
            },
            "multiplicity": "Benjamini-Hochberg and Bonferroni globally across every estimable Pfam-variable pair",
            "p001_threshold": P001_THRESHOLD,
            "fdr_threshold": FDR_THRESHOLD,
            "bonferroni_alpha": BONFERRONI_ALPHA,
            "panel_b_selection": "ten smallest raw p values for the recorded environment category; deterministic accession tie-break",
            "panel_c_selection": "all recorded environment category pairs with global BH q < 0.05",
            "clustering": "average-linkage Euclidean clustering of row-z-scored log2(raw count + 1) profiles; optimal ordering; no randomness",
        },
        "spearman_crosscheck": crosscheck,
        "results": {
            "global_test_family_size": int(total["total_tests"]),
            "p_lt_0_001": int(total["p_lt_0_001"]),
            "bh_q_lt_0_05": int(total["bh_q_lt_0_05"]),
            "bonferroni_p_lt_0_05": int(total["bonferroni_p_lt_0_05"]),
            "recorded_environment_category_bh_supported_pfams": int(
                len(panel["ordered_stats"])
            ),
        },
        "integrity": {
            "synthetic_or_simulated_data_used": False,
            "imputation_used": False,
            "subsampling_used": False,
            "functional_annotations_used": False,
            "unrecorded_numeric_environment_values_used": False,
            "outputs_reopened_and_row_counts_verified": True,
        },
        "outputs": output_records,
    }
    paths["provenance"].write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main() -> int:
    args = parse_args()
    run_id = args.run_id or run_stamp()
    if not re.fullmatch(r"\d{8}_\d{6}", run_id):
        raise ValueError("--run-id must use YYYYMMDD_HHMMSS")
    args.matrix = args.matrix.resolve()
    args.manifest = args.manifest.resolve()
    args.aef = args.aef.resolve()
    args.output_dir = args.output_dir.resolve()

    cohort, counts, pfams, validation = load_and_validate(
        args.matrix, args.manifest, args.aef
    )
    correlations, crosscheck = compute_correlations(cohort, counts, pfams)
    summary = summarize(correlations)
    panel = prepare_panel_data(cohort, counts, correlations)
    figure = build_figure(summary, panel, len(correlations))

    if args.validate_only:
        figure.canvas.draw()
        plt.close(figure)
        print(
            json.dumps(
                {
                    "status": "VALIDATION_PASSED_FULL_REAL_DATA_NO_OUTPUTS_WRITTEN",
                    "validation": validation,
                    "crosscheck": crosscheck,
                    "summary": summary.to_dict(orient="records"),
                    "panel_c_pfams": int(len(panel["ordered_stats"])),
                },
                indent=2,
            )
        )
        return 0

    paths = create_output_paths(args.output_dir, run_id)
    write_outputs(
        paths,
        figure,
        correlations,
        summary,
        panel,
        validation,
        crosscheck,
        args,
        run_id,
    )
    plt.close(figure)
    print(
        json.dumps(
            {
                "status": "COMPLETE",
                "run_id": run_id,
                "summary": summary.to_dict(orient="records"),
                "outputs": {key: project_path(path) for key, path in paths.items()},
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

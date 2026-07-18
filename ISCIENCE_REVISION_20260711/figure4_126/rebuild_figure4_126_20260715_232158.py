#!/usr/bin/env python3
"""Recompute the legacy Figure 4 A--P analysis on the authenticated 126 genomes.

Data provenance
---------------
The analysis uses the authenticated 126-genome Pfam count matrix (70
Rhodophyta, 43 Ochrophyta, and 13 Chlorophyta). The metadata input is the
reconciled source-inventory manifest; a validated one-to-one join on ``Genome``
selects those same 126 genomes. Only columns matching ``^PF\d{5}$`` are treated
as Pfam counts.

The statistical workflow intentionally reproduces the legacy analysis:

* global Pfam filter: sample variance >= 0.1 and prevalence >= 0.1;
* within-phylum Spearman rank correlations with recorded temperature, latitude,
  and longitude;
* signed two-sided Z scores combined by Stouffer's method with sqrt(n) weights;
* inverse-variance pooling of Fisher-transformed correlations, Cochran's Q, and
  I2 heterogeneity;
* nominal shared-direction class = Stouffer p < 0.05, same sign, and I2 < 50%;
* ERS = mean(r**2) over nominally significant environmental correlations.

No synthetic, simulated, randomly generated, subsampled, or reconstructed
summary data are used.  ``--validate-only`` executes the full-data in-memory
analysis and figure draw without writing result files, and reports observed
runtime and memory estimates for the complete run.

Created: 2026-07-15 23:21:58 +07:00
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import resource
import subprocess
import sys
import time
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
from matplotlib.collections import LineCollection, PathCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.text import Text
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import statsmodels


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
DEFAULT_MATRIX = (
    PROJECT_ROOT
    / "ISCIENCE_REVISION_20260711"
    / "analysis_stats"
    / "reconstructed_raw_pfam_counts_20260711_131706.csv.gz"
)
DEFAULT_METADATA = (
    PROJECT_ROOT
    / "ISCIENCE_REVISION_20260711"
    / "integrity"
    / "reconciled_analysis_manifest_20260711_110650.csv"
)
DEFAULT_OUTPUT_DIR = SCRIPT_PATH.parent

PFAM_PATTERN = re.compile(r"^PF\d{5}$")
EXPECTED_JOINED_N = 126
EXPECTED_PHYLA = OrderedDict(
    [("Rhodophyta", 70), ("Ochrophyta", 43), ("Chlorophyta", 13)]
)
PHYLA = list(EXPECTED_PHYLA)
ENVIRONMENTS = OrderedDict(
    [
        ("latitude", "DD latitude"),
        ("longitude", "DD longitude"),
        ("temperature", "Temperature (°C)"),
    ]
)
MIN_VARIANCE = 0.1
MIN_PREVALENCE = 0.1
MIN_SAMPLES = 8
MIN_PHYLA = 2
N_TOP = 12
N_HEATMAP = 15

PHYLUM_COLORS = {
    # Okabe--Ito colors; text labels provide redundant phylum encoding.
    "Rhodophyta": "#CC79A7",  # reddish purple
    "Ochrophyta": "#E69F00",  # orange
    "Chlorophyta": "#009E73",  # bluish green
}
ENV_COLORS = {
    # Okabe--Ito colors; lat/lon/temp are also named in text keys.
    "latitude": "#0072B2",  # blue
    "longitude": "#E69F00",  # orange
    "temperature": "#D55E00",  # vermillion
}
INTERPRETATION_COLORS = {
    "NOMINAL_SHARED_DIRECTION_LOW_HETEROGENEITY": "#009E73",
    "NOMINAL_SHARED_DIRECTION_HIGH_HETEROGENEITY": "#E69F00",
    "NOMINAL_MIXED_DIRECTION": "#CC79A7",
    "NOT_NOMINAL": "#999999",
}
INTERPRETATION_ORDER = [
    "NOMINAL_SHARED_DIRECTION_LOW_HETEROGENEITY",
    "NOMINAL_SHARED_DIRECTION_HIGH_HETEROGENEITY",
    "NOMINAL_MIXED_DIRECTION",
    "NOT_NOMINAL",
]

DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "correlation_okabe_ito_blue_white_vermillion",
    [
        "#0072B2",
        "#79B8D8",
        "#F2F2F2",
        "#ECA882",
        "#D55E00",
    ],
)
DIVERGING_CMAP.set_bad("white", alpha=1.0)


def configure_matplotlib() -> None:
    """Apply FIGURE_PROTOCOL.md typography and line settings."""
    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.family": "Arial",
            "font.sans-serif": ["Arial"],
            "font.size": 6,
            "axes.labelsize": 6,
            "axes.titlesize": 6,
            "axes.linewidth": 0.25,
            "axes.labelpad": 1,
            "axes.titlepad": 2,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "xtick.major.width": 0.25,
            "ytick.major.width": 0.25,
            "xtick.minor.width": 0.25,
            "ytick.minor.width": 0.25,
            "xtick.major.size": 2,
            "ytick.major.size": 2,
            "xtick.major.pad": 1,
            "ytick.major.pad": 1,
            "legend.fontsize": 6,
            "legend.frameon": False,
            "legend.handlelength": 1.0,
            "legend.handletextpad": 0.3,
            "legend.borderaxespad": 0.2,
            "lines.linewidth": 0.25,
            "lines.markeredgewidth": 0.25,
            "patch.linewidth": 0.25,
            "grid.linewidth": 0.25,
            "hatch.linewidth": 0.25,
            "savefig.transparent": True,
            "savefig.edgecolor": "none",
            "savefig.facecolor": "none",
        }
    )


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def realpath(path: Path) -> str:
    return str(path.expanduser().resolve())


def fail(message: str) -> None:
    raise RuntimeError(message)


def dataframe_bytes(frame: pd.DataFrame) -> int:
    return int(frame.memory_usage(index=True, deep=True).sum())


def peak_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes; Linux reports KiB.
    if sys.platform.startswith("darwin"):
        return int(value)
    return int(value * 1024)


def iso_now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def load_and_validate(
    matrix_path: Path, metadata_path: Path
) -> tuple[pd.DataFrame, list[str], pd.DataFrame, dict[str, Any]]:
    """Load complete inputs and enforce the authenticated 126-genome join."""
    if not matrix_path.is_file():
        fail(f"Count matrix not found: {matrix_path}")
    if not metadata_path.is_file():
        fail(f"Metadata manifest not found: {metadata_path}")

    matrix = pd.read_csv(matrix_path, low_memory=False)
    metadata = pd.read_csv(metadata_path, low_memory=False)

    if "Genome" not in matrix.columns or "Genome" not in metadata.columns:
        fail("Both inputs must contain a Genome column")
    if matrix["Genome"].isna().any() or metadata["Genome"].isna().any():
        fail("Genome identifiers must not be missing")
    if matrix["Genome"].duplicated().any():
        duplicates = matrix.loc[matrix["Genome"].duplicated(), "Genome"].tolist()
        fail(f"Count matrix Genome identifiers are not unique: {duplicates[:5]}")
    if metadata["Genome"].duplicated().any():
        duplicates = metadata.loc[metadata["Genome"].duplicated(), "Genome"].tolist()
        fail(f"Metadata Genome identifiers are not unique: {duplicates[:5]}")

    pfam_columns = [column for column in matrix.columns if PFAM_PATTERN.fullmatch(column)]
    excluded_matrix_columns = [
        column
        for column in matrix.columns
        if column != "Genome" and not PFAM_PATTERN.fullmatch(column)
    ]
    if not pfam_columns:
        fail("No count columns match the required ^PF\\d{5}$ pattern")
    if len(set(pfam_columns)) != len(pfam_columns):
        fail("Strict Pfam columns are not unique")

    counts = matrix.loc[:, pfam_columns]
    non_numeric = [column for column in pfam_columns if not pd.api.types.is_numeric_dtype(counts[column])]
    if non_numeric:
        fail(f"Strict Pfam columns must be numeric: {non_numeric[:5]}")
    count_values = counts.to_numpy(dtype=np.float64, copy=False)
    if not np.isfinite(count_values).all():
        fail("Raw Pfam count matrix contains missing or non-finite values")
    if np.any(count_values < 0):
        fail("Raw Pfam count matrix contains negative values")
    if not np.equal(count_values, np.floor(count_values)).all():
        fail("Raw Pfam count matrix contains non-integer values")

    required_metadata = [
        "Genome",
        "Phylum",
        *ENVIRONMENTS.values(),
        "raw_pfam_hit_total",
    ]
    missing_metadata = [column for column in required_metadata if column not in metadata.columns]
    if missing_metadata:
        fail(f"Metadata is missing required columns: {missing_metadata}")

    matrix_ids = set(matrix["Genome"])
    metadata_ids = set(metadata["Genome"])
    missing_ids = sorted(matrix_ids - metadata_ids)
    if missing_ids:
        fail(f"Matrix genomes absent from metadata: {missing_ids[:5]}")

    joined = matrix.loc[:, ["Genome", *pfam_columns]].merge(
        metadata.loc[:, required_metadata],
        on="Genome",
        how="inner",
        validate="one_to_one",
        sort=False,
    )
    if len(joined) != EXPECTED_JOINED_N:
        fail(f"Expected {EXPECTED_JOINED_N} joined genomes; observed {len(joined)}")

    phylum_counts = joined["Phylum"].value_counts().to_dict()
    if phylum_counts != dict(EXPECTED_PHYLA):
        fail(f"Expected phylum counts {dict(EXPECTED_PHYLA)}; observed {phylum_counts}")

    for env_column in ENVIRONMENTS.values():
        numeric_env = pd.to_numeric(joined[env_column], errors="coerce")
        if numeric_env.isna().any() or not np.isfinite(numeric_env.to_numpy()).all():
            fail(f"Joined environmental column is incomplete or non-numeric: {env_column}")
        if numeric_env.nunique() < 2:
            fail(f"Joined environmental column is constant: {env_column}")
        joined[env_column] = numeric_env.astype(float)

    observed_totals = joined.loc[:, pfam_columns].sum(axis=1).to_numpy(dtype=np.int64)
    expected_totals = pd.to_numeric(joined["raw_pfam_hit_total"], errors="raise").to_numpy(
        dtype=np.int64
    )
    if not np.array_equal(observed_totals, expected_totals):
        fail("Strict-Pfam row sums do not equal metadata raw_pfam_hit_total")

    extra_metadata_ids = sorted(metadata_ids - matrix_ids)
    summary = {
        "matrix_rows": int(len(matrix)),
        "metadata_rows": int(len(metadata)),
        "joined_rows": int(len(joined)),
        "metadata_rows_not_in_count_matrix": int(len(extra_metadata_ids)),
        "metadata_genomes_not_in_count_matrix": extra_metadata_ids,
        "phylum_counts": {key: int(phylum_counts[key]) for key in PHYLA},
        "strict_pfam_columns": int(len(pfam_columns)),
        "matrix_non_pfam_columns_excluded": excluded_matrix_columns,
        "all_environment_values_complete": True,
        "strict_pfam_row_sums_match_raw_hit_totals": True,
        "matrix_memory_bytes": dataframe_bytes(matrix),
        "metadata_memory_bytes": dataframe_bytes(metadata),
        "joined_memory_bytes": dataframe_bytes(joined),
    }
    return joined, pfam_columns, metadata, summary


def filter_pfams(
    joined: pd.DataFrame, pfam_columns: list[str]
) -> tuple[list[str], pd.DataFrame]:
    counts = joined.loc[:, pfam_columns]
    variance = counts.var(axis=0, ddof=1)
    prevalence = counts.gt(0).mean(axis=0)
    means = counts.mean(axis=0)
    keep = variance.ge(MIN_VARIANCE) & prevalence.ge(MIN_PREVALENCE)
    stats_frame = pd.DataFrame(
        {
            "pfam": pfam_columns,
            "sample_variance": variance.reindex(pfam_columns).to_numpy(),
            "nonzero_prevalence": prevalence.reindex(pfam_columns).to_numpy(),
            "mean_raw_count": means.reindex(pfam_columns).to_numpy(),
            "passes_legacy_filter": keep.reindex(pfam_columns).to_numpy(),
        }
    )
    filtered = stats_frame.loc[stats_frame["passes_legacy_filter"], "pfam"].tolist()
    if not filtered:
        fail("No Pfams pass the legacy global variance/prevalence filter")
    return filtered, stats_frame


def vectorized_spearman(
    x: np.ndarray, ranked_y_centered: np.ndarray, y_ss: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Spearman r and its legacy asymptotic two-sided p value for every column."""
    ranked_x = stats.rankdata(x, method="average")
    x_centered = ranked_x - ranked_x.mean()
    x_ss = float(np.dot(x_centered, x_centered))
    denom = np.sqrt(x_ss * y_ss)
    with np.errstate(divide="ignore", invalid="ignore"):
        r_values = np.dot(x_centered, ranked_y_centered) / denom
    r_values = np.clip(r_values, -1.0, 1.0)
    dof = len(x) - 2
    with np.errstate(divide="ignore", invalid="ignore"):
        t_values = r_values * np.sqrt(dof / ((1.0 + r_values) * (1.0 - r_values)))
    p_values = 2.0 * stats.t.sf(np.abs(t_values), df=dof)
    return r_values, p_values


def compute_correlations(
    joined: pd.DataFrame, filtered_pfams: list[str]
) -> tuple[pd.DataFrame, dict[str, float]]:
    frames: list[pd.DataFrame] = []
    crosscheck_r_error = 0.0
    crosscheck_p_error = 0.0

    for phylum in PHYLA:
        subset = joined.loc[joined["Phylum"].eq(phylum)].reset_index(drop=True)
        n_samples = len(subset)
        if n_samples < MIN_SAMPLES:
            fail(f"{phylum} has fewer than {MIN_SAMPLES} samples")

        y_values = subset.loc[:, filtered_pfams].to_numpy(dtype=np.float64, copy=True)
        ranked_y = stats.rankdata(y_values, axis=0, method="average")
        ranked_y_centered = ranked_y - ranked_y.mean(axis=0, keepdims=True)
        y_ss = np.einsum("ij,ij->j", ranked_y_centered, ranked_y_centered)
        variable = y_ss > 0
        variable_pfams = np.asarray(filtered_pfams, dtype=object)[variable]
        ranked_variable = ranked_y_centered[:, variable]
        variable_ss = y_ss[variable]

        for env_name, env_column in ENVIRONMENTS.items():
            x_values = subset[env_column].to_numpy(dtype=np.float64)
            r_values, p_values = vectorized_spearman(
                x_values, ranked_variable, variable_ss
            )
            frame = pd.DataFrame(
                {
                    "group": phylum,
                    "pfam": variable_pfams,
                    "env_var": env_name,
                    "r": r_values,
                    "p": p_values,
                    "r_squared": np.square(r_values),
                    "n": n_samples,
                }
            )
            if frame[["r", "p"]].isna().any().any():
                fail(f"Non-finite correlation result for {phylum}/{env_name}")
            frames.append(frame)

            # Deterministic validation against scipy.stats.spearmanr on real data.
            check_indices = sorted(set([0, len(variable_pfams) // 2, len(variable_pfams) - 1]))
            for index in check_indices:
                ref_r, ref_p = stats.spearmanr(x_values, y_values[:, np.flatnonzero(variable)[index]])
                crosscheck_r_error = max(crosscheck_r_error, abs(ref_r - r_values[index]))
                crosscheck_p_error = max(crosscheck_p_error, abs(ref_p - p_values[index]))

    correlations = pd.concat(frames, ignore_index=True)
    correlations["p_fdr"] = np.nan
    for phylum in PHYLA:
        mask = correlations["group"].eq(phylum)
        _, corrected = fdrcorrection(correlations.loc[mask, "p"].to_numpy())
        correlations.loc[mask, "p_fdr"] = corrected
    correlations["significant_nominal"] = correlations["p"].lt(0.05)
    correlations["significant_fdr"] = correlations["p_fdr"].lt(0.05)
    correlations = correlations.sort_values(
        ["group", "env_var", "p", "pfam"], kind="stable"
    ).reset_index(drop=True)

    if crosscheck_r_error > 1e-12 or crosscheck_p_error > 1e-12:
        fail(
            "Vectorized Spearman implementation did not match scipy.stats.spearmanr: "
            f"max r error={crosscheck_r_error}, max p error={crosscheck_p_error}"
        )
    crosscheck = {
        "max_absolute_r_error": float(crosscheck_r_error),
        "max_absolute_p_error": float(crosscheck_p_error),
        "comparisons": int(len(PHYLA) * len(ENVIRONMENTS) * 3),
    }
    return correlations, crosscheck


def meta_analyze(correlations: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    phylum_rank = {phylum: index for index, phylum in enumerate(PHYLA)}

    for (pfam, env_var), group in correlations.groupby(["pfam", "env_var"], sort=True):
        group = group.assign(_order=group["group"].map(phylum_rank)).sort_values("_order")
        if len(group) < MIN_PHYLA:
            continue

        r_values = group["r"].to_numpy(dtype=float)
        p_values = group["p"].to_numpy(dtype=float)
        n_values = group["n"].to_numpy(dtype=float)

        # Numerically stable equivalent of sign(r) * Phi^-1(1 - p/2).
        clipped_p = np.clip(p_values, np.finfo(float).tiny, 1.0)
        signed_z = np.sign(r_values) * stats.norm.isf(clipped_p / 2.0)
        stouffer_weights = np.sqrt(n_values)
        stouffer_z = float(
            np.dot(signed_z, stouffer_weights)
            / np.sqrt(np.dot(stouffer_weights, stouffer_weights))
        )
        stouffer_p = float(2.0 * stats.norm.sf(abs(stouffer_z)))

        fisher_z = np.arctanh(np.clip(r_values, -0.999, 0.999))
        fisher_weights = n_values - 3.0
        pooled_z = float(np.dot(fisher_z, fisher_weights) / fisher_weights.sum())
        pooled_r = float(np.tanh(pooled_z))
        q_stat = float(np.dot(fisher_weights, np.square(fisher_z - pooled_z)))
        q_df = len(group) - 1
        q_p = float(stats.chi2.sf(q_stat, q_df)) if q_df > 0 else math.nan
        i_squared = float(max(0.0, (q_stat - q_df) / q_stat)) if q_stat > 0 else 0.0

        n_positive = int(np.count_nonzero(r_values > 0))
        n_negative = int(np.count_nonzero(r_values < 0))
        same_direction = n_positive == 0 or n_negative == 0
        if stouffer_p < 0.05 and same_direction and i_squared < 0.5:
            interpretation = "NOMINAL_SHARED_DIRECTION_LOW_HETEROGENEITY"
        elif stouffer_p < 0.05 and same_direction:
            interpretation = "NOMINAL_SHARED_DIRECTION_HIGH_HETEROGENEITY"
        elif stouffer_p < 0.05:
            interpretation = "NOMINAL_MIXED_DIRECTION"
        else:
            interpretation = "NOT_NOMINAL"

        details = {row.group: row for row in group.itertuples(index=False)}
        output: dict[str, Any] = {
            "pfam": pfam,
            "env_var": env_var,
            "n_phyla": int(len(group)),
            "phyla_tested": ",".join(group["group"]),
            "pooled_r": pooled_r,
            "pooled_r_squared": pooled_r * pooled_r,
            "stouffer_z": stouffer_z,
            "stouffer_p": stouffer_p,
            "heterogeneity_Q": q_stat,
            "heterogeneity_df": int(q_df),
            "heterogeneity_p": q_p,
            "I_squared": i_squared,
            "same_direction": bool(same_direction),
            "n_positive": n_positive,
            "n_negative": n_negative,
            "interpretation": interpretation,
            "mean_within_phylum_r": float(np.mean(r_values)),
            "mean_within_phylum_r_squared": float(np.mean(np.square(r_values))),
        }
        for phylum in PHYLA:
            item = details.get(phylum)
            output[f"{phylum}_r"] = float(item.r) if item is not None else math.nan
            output[f"{phylum}_p"] = float(item.p) if item is not None else math.nan
            output[f"{phylum}_n"] = int(item.n) if item is not None else math.nan
        rows.append(output)

    meta = pd.DataFrame(rows)
    if meta.empty:
        fail("No Pfam-environment pairs were eligible for meta-analysis")
    _, meta["stouffer_p_fdr"] = fdrcorrection(meta["stouffer_p"].to_numpy())
    meta["nominal_stouffer_p_lt_0_05"] = meta["stouffer_p"].lt(0.05)
    meta["fdr_supported"] = meta["stouffer_p_fdr"].lt(0.05)
    meta["fdr_supported_shared_direction"] = meta["fdr_supported"] & meta[
        "interpretation"
    ].eq("NOMINAL_SHARED_DIRECTION_LOW_HETEROGENEITY")
    return meta.sort_values(["stouffer_p", "pfam", "env_var"], kind="stable").reset_index(
        drop=True
    )


def identify_specific(correlations: pd.DataFrame) -> pd.DataFrame:
    significant = correlations.loc[correlations["p"].lt(0.05)].copy()
    n_significant = significant.groupby(["pfam", "env_var"])["group"].transform("size")
    specific = significant.loc[n_significant.eq(1), ["pfam", "env_var", "group", "r", "p", "n"]]
    specific = specific.rename(columns={"group": "specific_to"})
    return specific.sort_values(["specific_to", "env_var", "p", "pfam"], kind="stable").reset_index(
        drop=True
    )


def calculate_ers(correlations: pd.DataFrame) -> pd.DataFrame:
    work = correlations.copy()
    work["significant_r_squared"] = work["r_squared"].where(work["p"].lt(0.05))
    work["abs_r"] = work["r"].abs()
    grouped = work.groupby(["group", "pfam"], sort=True)
    ers = grouped.agg(
        ers=("significant_r_squared", "mean"),
        n_sig_correlations=("significant_r_squared", "count"),
        n_tested=("env_var", "size"),
        max_abs_r=("abs_r", "max"),
    ).reset_index()
    ers["ers"] = ers["ers"].fillna(0.0)
    best_indices = grouped["abs_r"].idxmax().to_numpy()
    best = work.loc[best_indices, ["group", "pfam", "env_var"]].rename(
        columns={"env_var": "best_env_var"}
    )
    ers = ers.merge(best, on=["group", "pfam"], how="left", validate="one_to_one")
    return ers.rename(columns={"group": "phylum"}).sort_values(
        ["phylum", "ers", "pfam"], ascending=[True, False, True], kind="stable"
    ).reset_index(drop=True)


def axis_title(ax: mpl.axes.Axes, letter: str, descriptor: str = "") -> None:
    title = letter if not descriptor else f"{letter}  {descriptor}"
    ax.set_title(title, loc="left", fontweight="bold", fontsize=6, pad=2)


def style_axis(ax: mpl.axes.Axes, grid: bool = False) -> None:
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    ax.tick_params(which="both", width=0.25, labelsize=6)
    if grid:
        ax.grid(True, linewidth=0.25, alpha=0.18, color="#666666")
        ax.set_axisbelow(True)


def add_environment_scatter(
    ax: mpl.axes.Axes,
    frame: pd.DataFrame,
    x_column: str,
    y_column: str,
    size: float,
    alpha: float,
    labels: bool = False,
) -> None:
    for env_name in ENVIRONMENTS:
        subset = frame.loc[frame["env_var"].eq(env_name)]
        ax.scatter(
            subset[x_column],
            subset[y_column],
            s=size,
            alpha=alpha,
            color=ENV_COLORS[env_name],
            edgecolors="none",
            linewidths=0,
            marker="o",
            label=env_name[:3] if labels else None,
            rasterized=False,
        )


def build_figure(
    correlations: pd.DataFrame,
    meta: pd.DataFrame,
    ers: pd.DataFrame,
    specific: pd.DataFrame,
) -> mpl.figure.Figure:
    """Build the complete A--P figure entirely from computed real-data results."""
    configure_matplotlib()
    fig = plt.figure(figsize=(8.27, 10.15), layout="constrained")
    fig.patch.set_alpha(0.0)
    grid = fig.add_gridspec(
        5,
        4,
        height_ratios=[1.25, 1.05, 1.35, 1.35, 1.48],
        width_ratios=[1.05, 1.15, 1.15, 1.02],
        hspace=0.12,
        wspace=0.20,
    )

    # A: nonzero ERS distributions.
    ax = fig.add_subplot(grid[0, 0])
    nonzero_ers = ers.loc[ers["ers"].gt(0)]
    combined = nonzero_ers["ers"].to_numpy()
    bins = np.linspace(float(combined.min()), float(combined.max()), 26)
    for phylum in PHYLA:
        values = nonzero_ers.loc[nonzero_ers["phylum"].eq(phylum), "ers"]
        ax.hist(
            values,
            bins=bins,
            color=PHYLUM_COLORS[phylum],
            alpha=0.58,
            edgecolor="none",
            label=phylum[:5],
        )
    ax.set_xlabel("ERS (nominal p<0.05)")
    ax.set_ylabel("PFAMs")
    axis_title(ax, "A")
    ax.legend(loc="upper right")
    style_axis(ax)

    # B/C: strongest nominal shared-direction temperature and latitude results.
    shared_direction = meta.loc[
        meta["interpretation"].eq("NOMINAL_SHARED_DIRECTION_LOW_HETEROGENEITY")
    ]
    for column_index, (letter, env_name, descriptor) in enumerate(
        [
            ("B", "temperature", "Temp: nominal shared-dir."),
            ("C", "latitude", "Lat: nominal shared-dir."),
        ],
        start=1,
    ):
        ax = fig.add_subplot(grid[0, column_index])
        top = shared_direction.loc[shared_direction["env_var"].eq(env_name)].nsmallest(
            N_TOP, "stouffer_p"
        )
        y_positions = np.arange(len(top), dtype=float)
        for phylum_index, phylum in enumerate(PHYLA):
            values = top[f"{phylum}_r"].to_numpy(dtype=float)
            valid = np.isfinite(values)
            ax.barh(
                y_positions[valid] + (phylum_index - 1) * 0.23,
                values[valid],
                height=0.22,
                color=PHYLUM_COLORS[phylum],
                edgecolor="none",
                label=phylum[:5] if letter == "B" else None,
            )
        ax.set_yticks(y_positions)
        ax.set_yticklabels(top["pfam"].str.removeprefix("PF"))
        ax.invert_yaxis()
        ax.axvline(0, color="black", linewidth=0.25)
        ax.set_xlabel("r")
        axis_title(ax, letter, descriptor)
        style_axis(ax)

    # D: top nominal shared-direction cross-phylum heatmap.
    ax = fig.add_subplot(grid[0:2, 3])
    top_heatmap = shared_direction.nsmallest(N_HEATMAP, "stouffer_p").copy()
    heatmap = top_heatmap[[f"{phylum}_r" for phylum in PHYLA]].to_numpy(dtype=float)
    masked = np.ma.masked_invalid(heatmap)
    mesh = ax.pcolormesh(
        np.arange(4),
        np.arange(len(top_heatmap) + 1),
        masked,
        cmap=DIVERGING_CMAP,
        vmin=-0.6,
        vmax=0.6,
        shading="flat",
        edgecolors="none",
        linewidth=0,
        rasterized=False,
    )
    labels = [
        f"{row.pfam.removeprefix('PF')}_{row.env_var[:3]}"
        for row in top_heatmap.itertuples(index=False)
    ]
    ax.set_xticks(np.arange(3) + 0.5)
    ax.set_xticklabels(["Rhodo", "Ochro", "Chloro"], rotation=90)
    ax.set_yticks(np.arange(len(labels)) + 0.5)
    ax.set_yticklabels(labels)
    ax.set_ylim(len(labels), 0)
    axis_title(ax, "D", "Nominal shared-dir.")
    colorbar = fig.colorbar(mesh, ax=ax, fraction=0.08, pad=0.03)
    colorbar.set_label("r", rotation=0)
    colorbar.outline.set_linewidth(0.25)
    colorbar.dividers.set_linewidth(0.25)
    # Matplotlib otherwise auto-rasterizes colorbars with >50 gradient cells.
    colorbar.solids.set_rasterized(False)
    colorbar.ax.tick_params(width=0.25, labelsize=6)
    style_axis(ax)

    # E: fraction of tested Pfams with at least one nominal response.
    ax = fig.add_subplot(grid[1, 0])
    fractions = []
    for phylum in PHYLA:
        subset = ers.loc[ers["phylum"].eq(phylum)]
        fractions.append(float(100.0 * subset["ers"].gt(0).mean()))
    bars = ax.bar(
        np.arange(3),
        fractions,
        color=[PHYLUM_COLORS[phylum] for phylum in PHYLA],
        width=0.62,
        edgecolor="none",
    )
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels([phylum[:5] for phylum in PHYLA])
    ax.set_ylabel("% nominally responsive")
    ax.set_ylim(0, max(fractions) * 1.22)
    for bar, fraction in zip(bars, fractions):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(fractions) * 0.035,
            f"{fraction:.1f}%",
            ha="center",
            va="bottom",
        )
    axis_title(ax, "E")
    style_axis(ax)

    # F: nominal shared-direction effect magnitude versus heterogeneity.
    ax = fig.add_subplot(grid[1, 1])
    plot_frame = shared_direction.assign(
        abs_pooled_r=shared_direction["pooled_r"].abs(),
        I2_pct=100 * shared_direction["I_squared"],
    )
    add_environment_scatter(ax, plot_frame, "abs_pooled_r", "I2_pct", 3.0, 0.38, True)
    ax.axhline(50, color="black", linestyle="--", linewidth=0.25)
    ax.set_xlabel("|Pooled r|")
    ax.set_ylabel("I2 (%)")
    ax.set_ylim(-1, 52)
    axis_title(ax, "F", "Nominal shared-dir.")
    style_axis(ax)

    # G: signed Stouffer Z distribution for nominal shared-direction results.
    ax = fig.add_subplot(grid[1, 2])
    ax.hist(
        shared_direction["stouffer_z"],
        bins=40,
        color=INTERPRETATION_COLORS["NOMINAL_SHARED_DIRECTION_LOW_HETEROGENEITY"],
        alpha=0.78,
        edgecolor="none",
    )
    ax.axvline(-1.96, color="black", linestyle="--", linewidth=0.25)
    ax.axvline(1.96, color="black", linestyle="--", linewidth=0.25)
    ax.set_xlabel("Stouffer Z")
    ax.set_ylabel("Frequency")
    axis_title(ax, "G", "Nominal shared-dir.")
    style_axis(ax)

    # H--J: cross-phylum concordance of within-phylum r values.
    pairs = [
        ("Rhodophyta", "Ochrophyta"),
        ("Rhodophyta", "Chlorophyta"),
        ("Ochrophyta", "Chlorophyta"),
    ]
    pairwise_correlations: list[float] = []
    for index, (phylum_x, phylum_y) in enumerate(pairs):
        ax = fig.add_subplot(grid[2, index])
        col_x, col_y = f"{phylum_x}_r", f"{phylum_y}_r"
        subset = meta.loc[meta[col_x].notna() & meta[col_y].notna(), [col_x, col_y]]
        pearson_r = float(stats.pearsonr(subset[col_x], subset[col_y]).statistic)
        pairwise_correlations.append(pearson_r)
        limit = max(0.8, math.ceil(float(subset.abs().to_numpy().max()) * 10) / 10)
        limit = min(limit, 1.0)
        ax.scatter(
            subset[col_x],
            subset[col_y],
            s=1.1,
            alpha=0.16,
            color="#444444",
            edgecolors="none",
            linewidths=0,
            rasterized=False,
        )
        ax.plot([-limit, limit], [-limit, limit], color="#555555", linestyle="--", linewidth=0.25)
        ax.axhline(0, color="#777777", linewidth=0.25, alpha=0.5)
        ax.axvline(0, color="#777777", linewidth=0.25, alpha=0.5)
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        concordance_ticks = np.linspace(-limit, limit, 5)
        ax.set_xticks(concordance_ticks)
        # The lower-left x and y endpoint labels otherwise touch at 6 pt.
        ax.set_xticklabels(["", *[f"{value:.1f}" for value in concordance_ticks[1:]]])
        ax.set_yticks(concordance_ticks)
        ax.set_yticklabels(["", *[f"{value:.1f}" for value in concordance_ticks[1:]]])
        ax.set_xlabel(f"{phylum_x[:5]} r")
        ax.set_ylabel(f"{phylum_y[:5]} r")
        ax.text(0.04, 0.94, f"r={pearson_r:.2f}", transform=ax.transAxes, va="top")
        axis_title(ax, chr(ord("H") + index))
        style_axis(ax)

    # K--M: nominal within-phylum volcano panels.
    for index, phylum in enumerate(PHYLA):
        ax = fig.add_subplot(grid[3, index])
        subset = correlations.loc[correlations["group"].eq(phylum)].copy()
        subset["neg_log10_p"] = -np.log10(
            np.clip(subset["p"].to_numpy(dtype=float), np.finfo(float).tiny, 1.0)
        )
        add_environment_scatter(ax, subset, "r", "neg_log10_p", 1.2, 0.32, labels=index == 0)
        ax.axhline(-math.log10(0.05), color="black", linestyle="--", linewidth=0.25)
        ax.axvline(0, color="#777777", linewidth=0.25, alpha=0.5)
        ax.set_xlim(-1.0, 1.0)
        ax.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
        ax.set_xticklabels(["", "-0.5", "0.0", "0.5", "1.0"])
        ax.set_ylim(bottom=0.0)
        n_sig = int(subset["p"].lt(0.05).sum())
        ax.set_xlabel("r")
        ax.set_ylabel("-log10(p)" if index == 0 else "")
        axis_title(
            ax,
            chr(ord("K") + index),
            f"{phylum[:5]}: nominal n={n_sig:,}",
        )
        if index == 0:
            # The middle above the U-shaped point cloud is data-free.
            ax.legend(loc="upper center")
        style_axis(ax)

    # N: nominal associations significant in exactly one phylum.
    ax = fig.add_subplot(grid[3, 3])
    counts = (
        specific.groupby(["specific_to", "env_var"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=PHYLA, columns=list(ENVIRONMENTS), fill_value=0)
    )
    y = np.arange(len(PHYLA), dtype=float)
    for index, env_name in enumerate(ENVIRONMENTS):
        ax.barh(
            y + (index - 1) * 0.23,
            counts[env_name].to_numpy(),
            height=0.22,
            color=ENV_COLORS[env_name],
            edgecolor="none",
            label=env_name[:3],
        )
    ax.set_yticks(y)
    ax.set_yticklabels([phylum[:5] for phylum in PHYLA])
    ax.set_xlabel("Nominal detections\nin one tested phylum")
    axis_title(ax, "N")
    style_axis(ax)

    # O: full meta-analysis volcano and computed summary.
    ax = fig.add_subplot(grid[4, 0:3])
    volcano = meta.assign(
        neg_log10_p=-np.log10(
            np.clip(meta["stouffer_p"].to_numpy(dtype=float), np.finfo(float).tiny, 1.0)
        )
    )
    add_environment_scatter(ax, volcano, "pooled_r", "neg_log10_p", 1.6, 0.36, True)
    ax.axhline(-math.log10(0.05), color="black", linestyle="--", linewidth=0.25)
    ax.axvline(0, color="#777777", linewidth=0.25, alpha=0.5)
    significant_meta = meta.loc[meta["stouffer_p"].lt(0.05)]
    fdr_meta = meta.loc[meta["stouffer_p_fdr"].lt(0.05)]
    fdr_shared_direction = shared_direction.loc[
        shared_direction["stouffer_p_fdr"].lt(0.05)
    ]
    summary_text = (
        f"Eligible meta-tests (>=2 phyla): {len(meta):,}\n"
        f"Nominal p<0.05: {len(significant_meta):,} ({100 * len(significant_meta) / len(meta):.1f}%)\n"
        f"BH q<0.05: {len(fdr_meta):,} ({100 * len(fdr_meta) / len(meta):.1f}%)\n"
        f"Nominal shared-direction: {len(shared_direction):,} "
        f"({100 * len(shared_direction) / len(meta):.1f}%)\n"
        f"FDR-supported shared-direction: {len(fdr_shared_direction):,} "
        f"({100 * len(fdr_shared_direction) / len(meta):.1f}%)\n"
        f"Median |r| (nominal): {significant_meta['pooled_r'].abs().median():.3f}\n"
        f"Median I2 (nominal): {100 * significant_meta['I_squared'].median():.1f}%"
    )
    ax.text(0.27, 0.96, summary_text, transform=ax.transAxes, ha="left", va="top")
    ax.set_xlabel("Pooled r")
    ax.set_ylabel("-log10(p)")
    axis_title(ax, "O")
    style_axis(ax)

    # P: interpretation counts by recorded environmental variable.
    ax = fig.add_subplot(grid[4, 3])
    interpretation_counts = (
        meta.groupby(["env_var", "interpretation"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=list(ENVIRONMENTS), columns=INTERPRETATION_ORDER, fill_value=0)
    )
    x = np.arange(len(ENVIRONMENTS), dtype=float)
    bottom = np.zeros(len(ENVIRONMENTS), dtype=float)
    short_labels = {
        "NOMINAL_SHARED_DIRECTION_LOW_HETEROGENEITY": "Nom shared (I2<50%)",
        "NOMINAL_SHARED_DIRECTION_HIGH_HETEROGENEITY": "Nom shared (I2>=50%)",
        "NOMINAL_MIXED_DIRECTION": "Nom mixed-dir.",
        "NOT_NOMINAL": "Not nominal",
    }
    for category in INTERPRETATION_ORDER:
        values = interpretation_counts[category].to_numpy(dtype=float)
        ax.bar(
            x,
            values,
            bottom=bottom,
            width=0.64,
            color=INTERPRETATION_COLORS[category],
            edgecolor="none",
            label=short_labels[category],
        )
        bottom += values
    ax.set_xticks(x)
    ax.set_xticklabels(["Lat", "Lon", "Temp"])
    ax.set_ylabel("Eligible pairs (>=2 phyla)")
    axis_title(ax, "P")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, columnspacing=0.6)
    style_axis(ax)

    return fig


def bbox_intersection_area(a: mpl.transforms.Bbox, b: mpl.transforms.Bbox) -> float:
    width = max(0.0, min(a.x1, b.x1) - max(a.x0, b.x0))
    height = max(0.0, min(a.y1, b.y1) - max(a.y0, b.y0))
    return width * height


def audit_figure(fig: mpl.figure.Figure) -> dict[str, Any]:
    """Enforce font/line rules and detect text-text overlap after layout."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts: list[Text] = []
    font_files: set[str] = set()
    for item in fig.findobj(match=Text):
        if not item.get_visible() or not item.get_text().strip():
            continue
        if not math.isclose(float(item.get_fontsize()), 6.0, abs_tol=1e-9):
            fail(f"Figure text is not 6 pt: {item.get_text()!r} ({item.get_fontsize()})")
        font_file = mpl.font_manager.findfont(item.get_fontproperties(), fallback_to_default=False)
        font_files.add(realpath(Path(font_file)))
        if "arial" not in Path(font_file).name.lower():
            fail(f"Figure text did not resolve to Arial: {item.get_text()!r} -> {font_file}")
        texts.append(item)

    bad_lines: list[float] = []
    for item in fig.findobj(match=Line2D):
        if item.get_visible() and item.get_linewidth() not in (0, 0.25):
            bad_lines.append(float(item.get_linewidth()))
    for item in fig.findobj(match=LineCollection):
        if item.get_visible():
            for width in item.get_linewidths():
                if width not in (0, 0.25):
                    bad_lines.append(float(width))
    for item in fig.findobj(match=PathCollection):
        if item.get_visible():
            for width in item.get_linewidths():
                if width not in (0, 0.25):
                    bad_lines.append(float(width))
    if bad_lines:
        fail(f"Visible line widths other than 0 or 0.25 pt were found: {sorted(set(bad_lines))}")

    boxes: list[tuple[Text, mpl.transforms.Bbox]] = []
    for item in texts:
        box = item.get_window_extent(renderer=renderer)
        if box.width > 0 and box.height > 0:
            boxes.append((item, box))
    overlaps: list[dict[str, Any]] = []
    for first in range(len(boxes)):
        text_a, box_a = boxes[first]
        for second in range(first + 1, len(boxes)):
            text_b, box_b = boxes[second]
            area = bbox_intersection_area(box_a, box_b)
            if area > 0.5:
                overlaps.append(
                    {
                        "text_a": text_a.get_text(),
                        "text_b": text_b.get_text(),
                        "axes_a": text_a.axes.get_title() if text_a.axes is not None else None,
                        "axes_b": text_b.axes.get_title() if text_b.axes is not None else None,
                        "bbox_a": [float(value) for value in box_a.extents],
                        "bbox_b": [float(value) for value in box_b.extents],
                        "intersection_pixels_squared": float(area),
                    }
                )
    if overlaps:
        examples = overlaps[:8]
        fail(f"Detected {len(overlaps)} text-text overlaps; examples: {examples}")

    return {
        "text_elements_checked": int(len(texts)),
        "text_text_overlaps": 0,
        "font_size_points": 6.0,
        "resolved_font_files": sorted(font_files),
        "allowed_visible_line_widths_points": [0.0, 0.25],
    }


def output_paths(output_dir: Path, run_id: str) -> dict[str, Path]:
    return {
        "figure_pdf": output_dir / f"Figure4_126_A-P_{run_id}.pdf",
        "figure_svg": output_dir / f"Figure4_126_A-P_{run_id}.svg",
        "pfam_filter_stats": output_dir / f"Figure4_126_pfam_filter_stats_{run_id}.csv.gz",
        "within_phylum_correlations": output_dir
        / f"Figure4_126_within_phylum_correlations_{run_id}.csv.gz",
        "meta_analysis": output_dir / f"Figure4_126_meta_analysis_{run_id}.csv.gz",
        "ers": output_dir / f"Figure4_126_ers_{run_id}.csv.gz",
        "phylum_specific": output_dir / f"Figure4_126_phylum_specific_{run_id}.csv.gz",
        "provenance": output_dir / f"Figure4_126_provenance_{run_id}.json",
    }


def assert_non_overwriting(paths: dict[str, Path]) -> None:
    existing = [realpath(path) for path in paths.values() if path.exists()]
    if existing:
        fail(f"Refusing to overwrite existing outputs: {existing}")


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    if path.exists():
        fail(f"Refusing to overwrite: {path}")
    frame.to_csv(path, index=False, float_format="%.12g", compression="gzip")


def inspect_vector_outputs(pdf_path: Path, svg_path: Path) -> dict[str, Any]:
    svg_text = svg_path.read_text(encoding="utf-8")
    svg_image_elements = len(re.findall(r"<image\b", svg_text))
    if svg_image_elements:
        fail(f"SVG contains {svg_image_elements} raster image element(s)")

    pdfimages = subprocess.run(
        ["pdfimages", "-list", str(pdf_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    image_rows = [
        line
        for line in pdfimages.stdout.splitlines()
        if re.match(r"^\s*\d+\s+\d+\s+", line)
    ]
    if image_rows:
        fail(f"PDF contains raster images: {image_rows[:3]}")

    pdffonts = subprocess.run(
        ["pdffonts", str(pdf_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    font_rows = [
        line
        for line in pdffonts.stdout.splitlines()
        if line.strip() and not line.startswith("name") and not line.startswith("---")
    ]
    if not font_rows or any("Arial" not in line or " yes " not in f" {line} " for line in font_rows):
        fail(f"PDF font audit failed: {font_rows}")

    return {
        "pdf_raster_image_count": 0,
        "svg_image_element_count": 0,
        "pdf_font_rows": font_rows,
        "transparent_export_requested": True,
        "vector_formats": ["PDF", "SVG"],
    }


def safe_float(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not math.isfinite(float(value)) else float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, dict):
        return {str(key): safe_float(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe_float(item) for item in value]
    return value


def records(frame: pd.DataFrame, columns: list[str], n: int = N_TOP) -> list[dict[str, Any]]:
    return safe_float(frame.loc[:, columns].head(n).to_dict(orient="records"))


def summarize_results(
    correlations: pd.DataFrame,
    meta: pd.DataFrame,
    ers: pd.DataFrame,
    specific: pd.DataFrame,
) -> dict[str, Any]:
    within: dict[str, Any] = {}
    for phylum in PHYLA:
        subset = correlations.loc[correlations["group"].eq(phylum)]
        within[phylum] = {
            "tests": int(len(subset)),
            "distinct_pfams_tested": int(subset["pfam"].nunique()),
            "nominal_p_lt_0_05": int(subset["p"].lt(0.05).sum()),
            "fdr_q_lt_0_05": int(subset["p_fdr"].lt(0.05).sum()),
            "by_environment": {
                env_name: {
                    "tests": int(len(subset.loc[subset["env_var"].eq(env_name)])),
                    "nominal_p_lt_0_05": int(
                        subset.loc[subset["env_var"].eq(env_name), "p"].lt(0.05).sum()
                    ),
                    "fdr_q_lt_0_05": int(
                        subset.loc[subset["env_var"].eq(env_name), "p_fdr"].lt(0.05).sum()
                    ),
                }
                for env_name in ENVIRONMENTS
            },
        }

    nominal_meta = meta.loc[meta["stouffer_p"].lt(0.05)]
    shared_direction = meta.loc[
        meta["interpretation"].eq("NOMINAL_SHARED_DIRECTION_LOW_HETEROGENEITY")
    ]
    interpretation_counts = {
        category: int(meta["interpretation"].eq(category).sum())
        for category in INTERPRETATION_ORDER
    }
    fdr_cross_tab = {
        category: {
            "BH_q_lt_0_05": int(
                (meta["interpretation"].eq(category) & meta["stouffer_p_fdr"].lt(0.05)).sum()
            ),
            "BH_q_ge_0_05": int(
                (meta["interpretation"].eq(category) & meta["stouffer_p_fdr"].ge(0.05)).sum()
            ),
        }
        for category in INTERPRETATION_ORDER
    }
    by_environment: dict[str, Any] = {}
    for env_name in ENVIRONMENTS:
        subset = meta.loc[meta["env_var"].eq(env_name)]
        by_environment[env_name] = {
            "tests": int(len(subset)),
            "nominal_stouffer_p_lt_0_05": int(subset["stouffer_p"].lt(0.05).sum()),
            "BH_q_lt_0_05": int(subset["stouffer_p_fdr"].lt(0.05).sum()),
            "nominal_shared_direction_low_heterogeneity": int(
                subset["interpretation"]
                .eq("NOMINAL_SHARED_DIRECTION_LOW_HETEROGENEITY")
                .sum()
            ),
            "fdr_supported_shared_direction": int(
                (
                    subset["interpretation"].eq(
                        "NOMINAL_SHARED_DIRECTION_LOW_HETEROGENEITY"
                    )
                    & subset["stouffer_p_fdr"].lt(0.05)
                ).sum()
            ),
            "interpretation_counts": {
                category: int(subset["interpretation"].eq(category).sum())
                for category in INTERPRETATION_ORDER
            },
        }

    pairwise: dict[str, Any] = {}
    for phylum_x, phylum_y in [
        ("Rhodophyta", "Ochrophyta"),
        ("Rhodophyta", "Chlorophyta"),
        ("Ochrophyta", "Chlorophyta"),
    ]:
        col_x, col_y = f"{phylum_x}_r", f"{phylum_y}_r"
        subset = meta.loc[meta[col_x].notna() & meta[col_y].notna(), [col_x, col_y]]
        result = stats.pearsonr(subset[col_x], subset[col_y])
        pairwise[f"{phylum_x}_vs_{phylum_y}"] = {
            "pearson_r": float(result.statistic),
            "p": float(result.pvalue),
            "n_pairs": int(len(subset)),
        }

    specific_counts = (
        specific.groupby(["specific_to", "env_var"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=PHYLA, columns=list(ENVIRONMENTS), fill_value=0)
    )
    ers_summary: dict[str, Any] = {}
    for phylum in PHYLA:
        subset = ers.loc[ers["phylum"].eq(phylum)]
        ers_summary[phylum] = {
            "tested_pfams": int(len(subset)),
            "nominally_responsive_pfams": int(subset["ers"].gt(0).sum()),
            "nominally_responsive_percent": float(100 * subset["ers"].gt(0).mean()),
            "median_nonzero_ers": float(subset.loc[subset["ers"].gt(0), "ers"].median()),
            "maximum_ers": float(subset["ers"].max()),
            "maximum_ers_pfam": str(subset.loc[subset["ers"].idxmax(), "pfam"]),
        }

    top_columns = [
        "pfam",
        "env_var",
        "pooled_r",
        "stouffer_z",
        "stouffer_p",
        "stouffer_p_fdr",
        "I_squared",
        "Rhodophyta_r",
        "Ochrophyta_r",
        "Chlorophyta_r",
    ]
    return {
        "within_phylum": within,
        "meta_analysis": {
            "eligible_tests_at_least_2_phyla": int(len(meta)),
            "classification_basis": (
                "Nominal Stouffer p<0.05; the nominal shared-direction class additionally "
                "requires the same effect direction and I2<0.5. BH q is reported separately."
            ),
            "nominal_stouffer_p_lt_0_05": int(len(nominal_meta)),
            "BH_q_lt_0_05": int(meta["stouffer_p_fdr"].lt(0.05).sum()),
            "nominal_shared_direction_low_heterogeneity": int(len(shared_direction)),
            "fdr_supported_shared_direction": int(
                shared_direction["stouffer_p_fdr"].lt(0.05).sum()
            ),
            "interpretation_counts": interpretation_counts,
            "FDR_by_nominal_class_cross_tab": fdr_cross_tab,
            "median_absolute_pooled_r_among_nominal_p_lt_0_05": float(
                nominal_meta["pooled_r"].abs().median()
            ),
            "median_I_squared_among_nominal_p_lt_0_05": float(
                nominal_meta["I_squared"].median()
            ),
            "by_environment": by_environment,
            "pairwise_within_phylum_r_concordance": pairwise,
            "top_nominal_shared_direction_overall": records(shared_direction, top_columns),
            "top_nominal_shared_direction_temperature": records(
                shared_direction.loc[shared_direction["env_var"].eq("temperature")],
                top_columns,
            ),
            "top_nominal_shared_direction_latitude": records(
                shared_direction.loc[shared_direction["env_var"].eq("latitude")], top_columns
            ),
            "top_nominal_shared_direction_longitude": records(
                shared_direction.loc[shared_direction["env_var"].eq("longitude")], top_columns
            ),
        },
        "phylum_specific": {
            "classification_basis": (
                "Nominal p<0.05 detected in exactly one tested phylum; absence from another "
                "phylum may reflect either p>=0.05 or a non-variable, untested Pfam."
            ),
            "total": int(len(specific)),
            "counts": {
                phylum: {
                    env_name: int(specific_counts.loc[phylum, env_name])
                    for env_name in ENVIRONMENTS
                }
                for phylum in PHYLA
            },
            "totals_by_phylum": {
                phylum: int(specific["specific_to"].eq(phylum).sum()) for phylum in PHYLA
            },
        },
        "ers": ers_summary,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Run the complete full-data analysis and in-memory figure audit without writing outputs.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    run_started_at = iso_now()
    matrix_path = args.matrix.expanduser().resolve()
    metadata_path = args.metadata.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    print(f"Run started: {run_started_at}")
    print(f"Script: {SCRIPT_PATH}")
    print(f"Matrix: {matrix_path}")
    print(f"Metadata: {metadata_path}")
    print("Data policy: real full data only; no simulation, randomness, or subsampling")

    joined, pfam_columns, metadata, validation = load_and_validate(matrix_path, metadata_path)
    filtered_pfams, pfam_stats = filter_pfams(joined, pfam_columns)
    print(
        f"Validated join: n={len(joined)}; phyla={validation['phylum_counts']}; "
        f"strict Pfams={len(pfam_columns):,}; filtered Pfams={len(filtered_pfams):,}"
    )

    analysis_started = time.perf_counter()
    correlations, spearman_crosscheck = compute_correlations(joined, filtered_pfams)
    meta = meta_analyze(correlations)
    specific = identify_specific(correlations)
    ers = calculate_ers(correlations)
    analysis_seconds = time.perf_counter() - analysis_started

    figure_started = time.perf_counter()
    fig = build_figure(correlations, meta, ers, specific)
    figure_audit = audit_figure(fig)
    figure_draw_seconds = time.perf_counter() - figure_started
    numerical_summary = summarize_results(correlations, meta, ers, specific)

    print(
        "Computed full data: "
        f"within-phylum tests={len(correlations):,}; meta tests={len(meta):,}; "
        f"nominal Stouffer p<0.05="
        f"{numerical_summary['meta_analysis']['nominal_stouffer_p_lt_0_05']:,}; "
        f"nominal shared-direction="
        f"{numerical_summary['meta_analysis']['nominal_shared_direction_low_heterogeneity']:,}; "
        f"FDR-supported shared-direction="
        f"{numerical_summary['meta_analysis']['fdr_supported_shared_direction']:,}"
    )

    loaded_bytes = (
        validation["matrix_memory_bytes"]
        + validation["metadata_memory_bytes"]
        + validation["joined_memory_bytes"]
    )
    numeric_array_bytes = int(
        joined.loc[:, filtered_pfams].to_numpy(dtype=np.float64, copy=False).nbytes
    )
    result_memory_bytes = sum(
        dataframe_bytes(frame)
        for frame in [pfam_stats, correlations, meta, specific, ers]
    )
    estimated_peak_bytes = int(loaded_bytes + 4 * numeric_array_bytes + 2 * result_memory_bytes)
    elapsed_before_output = time.perf_counter() - started
    estimated_full_runtime_seconds = float(elapsed_before_output + 2 * figure_draw_seconds)
    resource_summary = {
        "input_dataframes_memory_bytes": int(loaded_bytes),
        "filtered_numeric_matrix_bytes": numeric_array_bytes,
        "result_dataframes_memory_bytes": int(result_memory_bytes),
        "estimated_peak_memory_bytes": estimated_peak_bytes,
        "observed_process_peak_rss_bytes": peak_rss_bytes(),
        "analysis_seconds": float(analysis_seconds),
        "figure_draw_and_audit_seconds": float(figure_draw_seconds),
        "elapsed_before_output_seconds": float(elapsed_before_output),
        "estimated_full_runtime_seconds": estimated_full_runtime_seconds,
        "runtime_estimate_formula": (
            "measured full-data validation elapsed + 2 * measured in-memory figure draw time "
            "for PDF and SVG serialization"
        ),
        "memory_estimate_formula": (
            "loaded input/join DataFrames + 4 filtered float matrices + 2 result DataFrame footprints"
        ),
    }

    if args.validate_only:
        plt.close(fig)
        print("VALIDATION PASSED (full 126-genome data; no output files written)")
        print(
            f"Measured analysis={analysis_seconds:.2f}s; draw/audit={figure_draw_seconds:.2f}s; "
            f"estimated full run={estimated_full_runtime_seconds:.2f}s"
        )
        print(
            f"Estimated peak memory={estimated_peak_bytes / 1024**2:.1f} MiB; "
            f"observed process peak RSS={resource_summary['observed_process_peak_rss_bytes'] / 1024**2:.1f} MiB"
        )
        return 0

    output_dir.mkdir(parents=True, exist_ok=True)
    run_id = datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")
    paths = output_paths(output_dir, run_id)
    assert_non_overwriting(paths)

    write_csv(pfam_stats, paths["pfam_filter_stats"])
    write_csv(correlations, paths["within_phylum_correlations"])
    write_csv(meta, paths["meta_analysis"])
    write_csv(ers, paths["ers"])
    write_csv(specific, paths["phylum_specific"])

    fig.savefig(
        paths["figure_pdf"],
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.03,
        transparent=True,
        edgecolor="none",
        metadata={
            "Title": "Figure 4: per-phylum recorded-environment analysis, authenticated 126 genomes",
            "Author": "Recomputed from authenticated project data",
            "Subject": f"Generated by {SCRIPT_PATH.name}",
            "CreationDate": datetime.now().astimezone(),
        },
    )
    fig.savefig(
        paths["figure_svg"],
        format="svg",
        bbox_inches="tight",
        pad_inches=0.03,
        transparent=True,
        edgecolor="none",
        metadata={
            "Title": "Figure 4: per-phylum recorded-environment analysis, authenticated 126 genomes",
            "Description": f"Generated by {SCRIPT_PATH.name}",
            "Date": iso_now(),
        },
    )
    plt.close(fig)
    vector_audit = inspect_vector_outputs(paths["figure_pdf"], paths["figure_svg"])

    formulas = {
        "global_filter": (
            "retain Pfam j when sample_variance(count_j, ddof=1) >= 0.1 and "
            "mean(count_j > 0) >= 0.1 across the 126 joined genomes"
        ),
        "spearman": "r_s = Pearson correlation of average ranks, computed separately within each phylum",
        "spearman_p": "t = r_s * sqrt((n-2)/((1+r_s)*(1-r_s))); two-sided p = 2*Pr(T[n-2] >= |t|)",
        "signed_z": "z_i = sign(r_i) * Phi^-1(1 - p_i/2), evaluated as sign(r_i)*norm.isf(p_i/2)",
        "weighted_stouffer": "Z = sum_i(sqrt(n_i)*z_i) / sqrt(sum_i(n_i)); p = 2*Phi(-|Z|)",
        "fisher_transform": "y_i = atanh(clip(r_i, -0.999, 0.999)); Var(y_i) = 1/(n_i-3)",
        "pooled_effect": "y_pool = sum_i((n_i-3)*y_i)/sum_i(n_i-3); r_pool = tanh(y_pool)",
        "cochran_Q": "Q = sum_i((n_i-3)*(y_i-y_pool)^2); df = k-1",
        "I_squared": "I2 = max(0, (Q-df)/Q), with I2=0 when Q=0",
        "nominal_shared_direction_class": (
            "Stouffer nominal p < 0.05, all nonzero r_i have the same direction, and I2 < 0.5; "
            "this is an evidence-screening label, not a claim of evolutionary convergence"
        ),
        "meta_analysis_eligibility": "a Pfam-environment pair must be tested in at least 2 phyla",
        "ers": "ERS_j,phylum = mean(r_s^2) over environment tests with nominal p < 0.05; otherwise 0",
        "multiple_testing": (
            "Benjamini-Hochberg FDR is computed within phylum over all its environment tests and "
            "separately over all eligible meta-analysis tests"
        ),
    }
    versions = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "statsmodels": statsmodels.__version__,
        "matplotlib": mpl.__version__,
    }

    output_hashes = {
        key: {
            "path": realpath(path),
            "bytes": int(path.stat().st_size),
            "sha256": sha256_file(path),
        }
        for key, path in paths.items()
        if key != "provenance"
    }
    provenance = {
        "run": {
            "started_at": run_started_at,
            "completed_at": iso_now(),
            "run_id": run_id,
            "validate_only": False,
            "randomness_used": False,
            "subsampling_used": False,
            "synthetic_data_used": False,
        },
        "script": {
            "path": realpath(SCRIPT_PATH),
            "sha256": sha256_file(SCRIPT_PATH),
            "created_timestamp_in_filename": "20260715_232158",
        },
        "inputs": {
            "raw_pfam_counts": {
                "path": realpath(matrix_path),
                "bytes": int(matrix_path.stat().st_size),
                "sha256": sha256_file(matrix_path),
            },
            "reconciled_metadata": {
                "path": realpath(metadata_path),
                "bytes": int(metadata_path.stat().st_size),
                "sha256": sha256_file(metadata_path),
            },
        },
        "parameters": {
            "pfam_column_regex": PFAM_PATTERN.pattern,
            "minimum_global_sample_variance": MIN_VARIANCE,
            "minimum_global_nonzero_prevalence": MIN_PREVALENCE,
            "minimum_within_phylum_samples": MIN_SAMPLES,
            "minimum_meta_analysis_phyla": MIN_PHYLA,
            "nominal_alpha": 0.05,
            "nominal_shared_direction_I_squared_threshold": 0.5,
            "environment_columns": dict(ENVIRONMENTS),
            "phylum_order": PHYLA,
            "color_palette": {
                "name": "Okabe-Ito",
                "phylum_colors": PHYLUM_COLORS,
                "environment_colors": ENV_COLORS,
                "classification_colors": INTERPRETATION_COLORS,
                "redundant_text_labels": True,
            },
        },
        "validation": {
            **validation,
            "filtered_pfam_columns": int(len(filtered_pfams)),
            "spearman_crosscheck": spearman_crosscheck,
        },
        "formulas": formulas,
        "software_versions": versions,
        "resource_accounting": resource_summary,
        "numerical_summary": numerical_summary,
        "figure_quality_audit": {**figure_audit, **vector_audit},
        "outputs": output_hashes,
        "provenance_note": (
            "The provenance JSON cannot contain its own SHA-256 digest without recursion; "
            "all scientific result and figure outputs are hashed here."
        ),
    }
    if paths["provenance"].exists():
        fail(f"Refusing to overwrite: {paths['provenance']}")
    with paths["provenance"].open("x", encoding="utf-8") as handle:
        json.dump(safe_float(provenance), handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")

    print("FULL RUN PASSED")
    for key, path in paths.items():
        print(f"{key}: {realpath(path)}")
    print(f"provenance_sha256: {sha256_file(paths['provenance'])}")
    print(f"Total wall time: {time.perf_counter() - started:.2f}s")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise

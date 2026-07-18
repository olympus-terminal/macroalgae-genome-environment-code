#!/usr/bin/env python3
"""Recompute the full-cohort Pfam--AEF screen from validated real inputs.

The inferential family is all variable, strict ``PF\d{5}`` count profiles
against A00--A63 in the exact 126-genome cohort. Spearman correlations use
average ranks, two-sided asymptotic P values, and one global Benjamini--
Hochberg correction. Inferential outputs and the Figure 5 display retain every
Pfam. A complete exact-duplicate membership table is written as an audit; it
does not alter inference or display weighting.

No synthetic, imputed, randomly generated, or hard-coded result values are
used. Every reported number is computed from the three supplied input files.

Created: 2026-07-18
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import resource
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist
from statsmodels.stats.multitest import multipletests
import statsmodels


SCRIPT_VERSION = "2026-07-18.6"
PFAM_RE = re.compile(r"^PF\d{5}$")
AEF_RE = re.compile(r"^A\d{2}$")
EXPECTED_AXES = [f"A{i:02d}" for i in range(64)]
EXPECTED_COHORT_SIZE = 126


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--pfam-counts", type=Path, required=True)
    parser.add_argument("--aef-embeddings", type=Path, required=True)
    parser.add_argument(
        "--output-parent",
        type=Path,
        required=True,
        help="A new timestamped run directory is created below this directory.",
    )
    parser.add_argument(
        "--skip-figure",
        action="store_true",
        help="Write inferential outputs without hierarchical clustering or Figure 5.",
    )
    return parser.parse_args()


def sha256(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def as_bool(series: pd.Series, name: str) -> pd.Series:
    mapping = {
        True: True,
        False: False,
        1: True,
        0: False,
        "true": True,
        "false": False,
        "True": True,
        "False": False,
        "TRUE": True,
        "FALSE": False,
    }
    parsed = series.map(mapping)
    if parsed.isna().any():
        bad = series.loc[parsed.isna()].astype(str).unique().tolist()
        raise ValueError(f"{name} contains unrecognized Boolean values: {bad}")
    return parsed.astype(bool)


def require_unique_nonblank(df: pd.DataFrame, label: str) -> None:
    if "Genome" not in df.columns:
        raise ValueError(f"{label} has no Genome column")
    values = df["Genome"].astype("string")
    if values.isna().any() or values.str.strip().eq("").any():
        raise ValueError(f"{label} contains a blank Genome identifier")
    if values.duplicated().any():
        duplicates = values.loc[values.duplicated(keep=False)].unique().tolist()
        raise ValueError(f"{label} contains duplicate Genome identifiers: {duplicates}")


def load_validated_inputs(
    manifest_path: Path, pfam_path: Path, aef_path: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    for path in (manifest_path, pfam_path, aef_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    manifest_all = pd.read_csv(manifest_path, low_memory=False)
    counts_all = pd.read_csv(pfam_path, low_memory=False)
    aef_all = pd.read_csv(aef_path, low_memory=False)
    for frame, label in (
        (manifest_all, "manifest"),
        (counts_all, "Pfam count table"),
        (aef_all, "AEF embedding table"),
    ):
        require_unique_nonblank(frame, label)

    required_manifest = {
        "master_row",
        "Genome",
        "Species",
        "Phylum",
        "DD latitude",
        "DD longitude",
        "aef_present",
        "aef_id_number",
        "safe_for_raw_pfam_analysis",
        "safe_for_aef_pfam_analysis",
    }
    missing = sorted(required_manifest - set(manifest_all.columns))
    if missing:
        raise ValueError(f"Manifest is missing columns: {missing}")

    safe = as_bool(
        manifest_all["safe_for_aef_pfam_analysis"], "safe_for_aef_pfam_analysis"
    )
    manifest = manifest_all.loc[safe].sort_values("master_row", kind="stable").copy()
    if len(manifest) != EXPECTED_COHORT_SIZE:
        raise ValueError(
            f"Expected {EXPECTED_COHORT_SIZE} safe AEF genomes, found {len(manifest)}"
        )
    if not as_bool(manifest["aef_present"], "aef_present").all():
        raise ValueError("A safe AEF record is not marked aef_present")
    if not as_bool(
        manifest["safe_for_raw_pfam_analysis"], "safe_for_raw_pfam_analysis"
    ).all():
        raise ValueError("A safe AEF record is not safe_for_raw_pfam_analysis")

    expected_ids = manifest["Genome"].astype(str).tolist()
    expected_set = set(expected_ids)
    for frame, label in ((counts_all, "Pfam counts"), (aef_all, "AEF embeddings")):
        observed = set(frame["Genome"].astype(str))
        if observed != expected_set:
            raise ValueError(
                f"{label} Genome set differs from safe manifest: "
                f"missing={sorted(expected_set-observed)}, extra={sorted(observed-expected_set)}"
            )

    counts = counts_all.set_index("Genome").loc[expected_ids].reset_index()
    aef = aef_all.set_index("Genome").loc[expected_ids].reset_index()
    manifest["Genome"] = manifest["Genome"].astype(str)
    counts["Genome"] = counts["Genome"].astype(str)
    aef["Genome"] = aef["Genome"].astype(str)
    if not (
        manifest["Genome"].tolist()
        == counts["Genome"].tolist()
        == aef["Genome"].tolist()
    ):
        raise AssertionError("Exact-ID ordering failed")

    axes = [column for column in aef.columns if AEF_RE.fullmatch(str(column))]
    if axes != EXPECTED_AXES:
        raise ValueError(f"Expected exactly A00--A63 in order; observed {axes}")
    pfams = [column for column in counts.columns if PFAM_RE.fullmatch(str(column))]
    if not pfams:
        raise ValueError("No strict PFxxxxx columns were found")
    if len(set(pfams)) != len(pfams):
        raise ValueError("Duplicate Pfam accession columns were found")

    count_values = counts[pfams].apply(pd.to_numeric, errors="raise").to_numpy(float)
    aef_values = aef[axes].apply(pd.to_numeric, errors="raise").to_numpy(float)
    if not np.isfinite(count_values).all() or not np.isfinite(aef_values).all():
        raise ValueError("Non-finite Pfam or AEF values were found")
    if (count_values < 0).any() or not np.equal(count_values, np.floor(count_values)).all():
        raise ValueError("Pfam counts must be nonnegative integers")
    if (count_values.sum(axis=1) == 0).any():
        bad = counts.loc[count_values.sum(axis=1) == 0, "Genome"].tolist()
        raise ValueError(f"All-zero genome Pfam profiles were found: {bad}")
    variable = np.ptp(count_values, axis=0) > 0
    if not variable.all():
        pfams = [pfam for pfam, keep in zip(pfams, variable) if keep]
        count_values = count_values[:, variable]
        counts = pd.concat(
            [counts[["Genome"]].reset_index(drop=True), pd.DataFrame(count_values, columns=pfams)],
            axis=1,
        )
    if not (np.ptp(aef_values, axis=0) > 0).all():
        raise ValueError("At least one AEF axis is invariant")

    if "ID number" not in aef.columns:
        raise ValueError("AEF table lacks ID number")
    if not np.array_equal(
        pd.to_numeric(aef["ID number"], errors="raise").to_numpy(float),
        pd.to_numeric(manifest["aef_id_number"], errors="raise").to_numpy(float),
    ):
        raise ValueError("AEF ID number differs from reconciled manifest")
    for column in ("Species", "DD latitude", "DD longitude"):
        if column not in aef.columns:
            raise ValueError(f"AEF table lacks {column}")
    if not np.array_equal(
        aef["Species"].fillna("").astype(str).to_numpy(),
        manifest["Species"].fillna("").astype(str).to_numpy(),
    ):
        raise ValueError("AEF Species labels differ from reconciled manifest")
    for column in ("DD latitude", "DD longitude"):
        if not np.allclose(
            pd.to_numeric(aef[column], errors="raise").to_numpy(float),
            pd.to_numeric(manifest[column], errors="raise").to_numpy(float),
            rtol=0,
            atol=1e-10,
        ):
            raise ValueError(f"AEF {column} differs from reconciled manifest")

    counts[pfams] = count_values.astype(np.int64)
    aef[axes] = aef_values
    return manifest.reset_index(drop=True), counts, aef, pfams


def compute_screen(
    count_values: np.ndarray,
    aef_values: np.ndarray,
    pfams: list[str],
    axes: list[str],
) -> dict[str, object]:
    if count_values.shape != (aef_values.shape[0], len(pfams)):
        raise ValueError("Count matrix shape does not match Pfam labels or AEF rows")
    if aef_values.shape[1] != len(axes):
        raise ValueError("AEF matrix shape does not match axis labels")
    n = count_values.shape[0]
    if n < 3:
        raise ValueError("At least three records are required")

    count_ranks = stats.rankdata(count_values, axis=0, method="average")
    aef_ranks = stats.rankdata(aef_values, axis=0, method="average")
    count_centered = count_ranks - count_ranks.mean(axis=0, keepdims=True)
    aef_centered = aef_ranks - aef_ranks.mean(axis=0, keepdims=True)
    count_ss = np.einsum("ij,ij->j", count_centered, count_centered)
    aef_ss = np.einsum("ij,ij->j", aef_centered, aef_centered)
    if (count_ss <= 0).any() or (aef_ss <= 0).any():
        raise ValueError("Invariant vectors reached the Spearman calculation")
    rho = (count_centered.T @ aef_centered) / np.sqrt(
        count_ss[:, None] * aef_ss[None, :]
    )
    rho = np.clip(rho, -1.0, 1.0)
    dof = n - 2
    denominator = np.maximum((1.0 - rho) * (1.0 + rho), np.finfo(float).tiny)
    t_stat = rho * np.sqrt(dof / denominator)
    p_value = 2.0 * stats.t.sf(np.abs(t_stat), df=dof)
    flat_p = p_value.ravel(order="C")
    if not np.isfinite(flat_p).all():
        raise ValueError("A non-finite P value was produced")
    q_value = multipletests(flat_p, method="fdr_bh")[1].reshape(rho.shape)

    row_mean = rho.mean(axis=1, keepdims=True)
    row_sd = rho.std(axis=1, ddof=1, keepdims=True)
    if (row_sd <= 0).any():
        raise ValueError("A constant 64-axis correlation profile cannot be z scored")
    row_z = (rho - row_mean) / row_sd
    return {
        "n": n,
        "rho": rho,
        "p": p_value,
        "q": q_value,
        "z": row_z,
    }


def long_table(result: dict[str, object], pfams: list[str], axes: list[str]) -> pd.DataFrame:
    rho = np.asarray(result["rho"])
    p_value = np.asarray(result["p"])
    q_value = np.asarray(result["q"])
    return pd.DataFrame(
        {
            "pfam": np.repeat(np.asarray(pfams, dtype=object), len(axes)),
            "embedding_dim": np.tile(np.asarray(axes, dtype=object), len(pfams)),
            "spearman_rho": rho.ravel(order="C"),
            "p_value_two_sided": p_value.ravel(order="C"),
            "q_value_bh_global": q_value.ravel(order="C"),
            "n_genomes": int(result["n"]),
        }
    )


def validate_sentinels(
    count_values: np.ndarray,
    aef_values: np.ndarray,
    pfams: list[str],
    axes: list[str],
    result: dict[str, object],
) -> pd.DataFrame:
    requested = [
        (pfams[0], axes[0]),
        (pfams[len(pfams) // 2], axes[31]),
        (pfams[-1], axes[-1]),
        ("PF01638", "A52"),
        ("PF01638", "A53"),
        ("PF10988", "A36"),
        ("PF13411", "A18"),
        ("PF00092", "A06"),
    ]
    rows: list[dict[str, object]] = []
    rho = np.asarray(result["rho"])
    p_value = np.asarray(result["p"])
    for pfam, axis in requested:
        if pfam not in pfams or axis not in axes:
            continue
        i, j = pfams.index(pfam), axes.index(axis)
        scipy_result = stats.spearmanr(count_values[:, i], aef_values[:, j])
        delta_rho = float(rho[i, j] - scipy_result.statistic)
        delta_p = float(p_value[i, j] - scipy_result.pvalue)
        if abs(delta_rho) > 1e-12 or abs(delta_p) > 1e-12:
            raise AssertionError(
                f"Vectorized result differs from scipy.stats.spearmanr for {pfam}/{axis}"
            )
        rows.append(
            {
                "pfam": pfam,
                "embedding_dim": axis,
                "vectorized_rho": rho[i, j],
                "scipy_rho": scipy_result.statistic,
                "absolute_rho_difference": abs(delta_rho),
                "vectorized_p": p_value[i, j],
                "scipy_p": scipy_result.pvalue,
                "absolute_p_difference": abs(delta_p),
            }
        )
    return pd.DataFrame(rows)


def write_screen_outputs(
    prefix: str,
    out_dir: Path,
    result: dict[str, object],
    pfams: list[str],
    axes: list[str],
    prevalence: np.ndarray,
) -> dict[str, Path]:
    rho = np.asarray(result["rho"])
    p_value = np.asarray(result["p"])
    q_value = np.asarray(result["q"])
    row_z = np.asarray(result["z"])
    full = long_table(result, pfams, axes)
    paths = {
        "full": out_dir / f"{prefix}_full_correlations.csv.gz",
        "significant": out_dir / f"{prefix}_significant_q_lt_0.05.csv",
        "rho": out_dir / f"{prefix}_correlation_matrix.csv.gz",
        "z": out_dir / f"{prefix}_row_z_matrix.csv.gz",
        "dimensions": out_dir / f"{prefix}_dimension_summary.csv",
        "pfams": out_dir / f"{prefix}_pfam_summary.csv.gz",
    }
    full.to_csv(paths["full"], index=False, compression="gzip", float_format="%.17g")
    full.loc[full["q_value_bh_global"] < 0.05].to_csv(
        paths["significant"], index=False, float_format="%.17g"
    )
    pd.DataFrame(rho, index=pfams, columns=axes).rename_axis("pfam").to_csv(
        paths["rho"], compression="gzip", float_format="%.17g"
    )
    pd.DataFrame(row_z, index=pfams, columns=axes).rename_axis("pfam").to_csv(
        paths["z"], compression="gzip", float_format="%.17g"
    )
    dimension_summary = pd.DataFrame(
        {
            "embedding_dim": axes,
            "n_tests": len(pfams),
            "n_p_lt_0.05": (p_value < 0.05).sum(axis=0),
            "n_p_lt_0.001": (p_value < 0.001).sum(axis=0),
            "n_q_lt_0.05": (q_value < 0.05).sum(axis=0),
            "minimum_q": q_value.min(axis=0),
            "minimum_rho": rho.min(axis=0),
            "maximum_rho": rho.max(axis=0),
        }
    )
    dimension_summary.to_csv(paths["dimensions"], index=False, float_format="%.17g")
    pfam_summary = pd.DataFrame(
        {
            "pfam": pfams,
            "n_nonzero_genomes": prevalence,
            "n_q_lt_0.05": (q_value < 0.05).sum(axis=1),
            "minimum_q": q_value.min(axis=1),
            "minimum_rho": rho.min(axis=1),
            "maximum_rho": rho.max(axis=1),
            "maximum_absolute_rho": np.abs(rho).max(axis=1),
        }
    )
    pfam_summary.to_csv(paths["pfams"], index=False, compression="gzip", float_format="%.17g")
    return paths


def duplicate_profile_map(
    pfams: list[str], row_z: np.ndarray
) -> tuple[np.ndarray, pd.DataFrame]:
    representative_by_bytes: dict[bytes, int] = {}
    unique_indices: list[int] = []
    rows: list[dict[str, object]] = []
    for index, pfam in enumerate(pfams):
        key = np.ascontiguousarray(row_z[index], dtype=np.float64).tobytes()
        representative_index = representative_by_bytes.get(key)
        if representative_index is None:
            representative_by_bytes[key] = index
            representative_index = index
            unique_indices.append(index)
        rows.append(
            {
                "pfam": pfam,
                "representative_pfam": pfams[representative_index],
                "is_representative": index == representative_index,
            }
        )
    return np.asarray(unique_indices, dtype=int), pd.DataFrame(rows)


def configure_figure_style() -> LinearSegmentedColormap:
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
        }
    )
    colors = [
        (0.0, 0.3, 0.7),
        (0.3, 0.5, 0.9),
        (0.7, 0.8, 0.95),
        (0.95, 0.95, 0.95),
        (0.95, 0.8, 0.7),
        (0.9, 0.4, 0.3),
        (0.7, 0.1, 0.1),
    ]
    cmap = LinearSegmentedColormap.from_list("protocol_diverging", colors)
    cmap.set_bad("white")
    return cmap


def cluster_and_plot(
    out_dir: Path,
    pfams: list[str],
    axes: list[str],
    result: dict[str, object],
) -> dict[str, Path]:
    rho = np.asarray(result["rho"])
    q_value = np.asarray(result["q"])
    row_z = np.asarray(result["z"])
    unique_indices, membership = duplicate_profile_map(pfams, row_z)
    if row_z.shape[0] < 2:
        raise ValueError("At least two profiles are required for clustering")

    row_dist = pdist(row_z, metric="correlation")
    col_dist = pdist(row_z.T, metric="correlation")
    if not np.isfinite(row_dist).all() or not np.isfinite(col_dist).all():
        raise ValueError("Non-finite correlation distance was produced")
    row_linkage = linkage(row_dist, method="average")
    col_linkage = linkage(col_dist, method="average")
    row_order = leaves_list(row_linkage)
    col_order = leaves_list(col_linkage)

    paths = {
        "membership": out_dir / "pooled_duplicate_profile_membership.csv.gz",
        "row_order": out_dir / "figure5_display_row_order.csv",
        "column_order": out_dir / "figure5_display_column_order.csv",
        "display_z": out_dir / "figure5_all_profiles_row_z_matrix.csv.gz",
        "pdf": out_dir / "Figure5_corrected_AEF_Pfam_landscape.pdf",
        "svg": out_dir / "Figure5_corrected_AEF_Pfam_landscape.svg",
    }
    membership.to_csv(paths["membership"], index=False, compression="gzip")
    pd.DataFrame(
        {
            "display_position": np.arange(len(row_order)),
            "pfam": np.asarray(pfams, dtype=object)[row_order],
            "original_pfam_index": row_order,
        }
    ).to_csv(paths["row_order"], index=False)
    pd.DataFrame(
        {
            "display_position": np.arange(len(col_order)),
            "embedding_dim": np.asarray(axes, dtype=object)[col_order],
        }
    ).to_csv(paths["column_order"], index=False)
    pd.DataFrame(row_z, index=pfams, columns=axes).rename_axis("pfam").to_csv(
        paths["display_z"], compression="gzip", float_format="%.17g"
    )

    cmap = configure_figure_style()
    fig = plt.figure(figsize=(8.27, 10.2))
    layout = gridspec.GridSpec(
        3,
        4,
        figure=fig,
        height_ratios=[0.62, 0.36, 8.2],
        width_ratios=[0.65, 7.2, 0.30, 0.16],
        hspace=0.02,
        wspace=0.03,
    )
    ax_col_tree = fig.add_subplot(layout[0, 1])
    dendrogram(
        col_linkage,
        ax=ax_col_tree,
        no_labels=True,
        color_threshold=0,
        above_threshold_color="black",
    )
    for collection in ax_col_tree.collections:
        collection.set_linewidth(0.25)
    ax_col_tree.set_axis_off()

    ax_counts = fig.add_subplot(layout[1, 1])
    axis_sig = (q_value < 0.05).sum(axis=0)[col_order]
    ax_counts.bar(
        np.arange(len(axes)), axis_sig, width=0.9, color="#315a8a", edgecolor="none"
    )
    ax_counts.set_xlim(-0.5, len(axes) - 0.5)
    ax_counts.set_ylabel("BH q<0.05")
    ax_counts.set_xticks([])
    ax_counts.spines[["top", "right"]].set_visible(False)

    ax_row_tree = fig.add_subplot(layout[2, 0])
    dendrogram(
        row_linkage,
        ax=ax_row_tree,
        orientation="left",
        no_labels=True,
        color_threshold=0,
        above_threshold_color="black",
    )
    for collection in ax_row_tree.collections:
        collection.set_linewidth(0.25)
    ax_row_tree.set_axis_off()

    ax_heat = fig.add_subplot(layout[2, 1])
    display_z = row_z[row_order][:, col_order]
    image = ax_heat.imshow(
        display_z,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=-3,
        vmax=3,
        origin="upper",
    )
    ax_heat.set_xticks(np.arange(len(axes)))
    ax_heat.set_xticklabels(np.asarray(axes, dtype=object)[col_order], rotation=90)
    ax_heat.xaxis.tick_bottom()
    ax_heat.set_yticks([])
    for spine in ax_heat.spines.values():
        spine.set_linewidth(0.25)

    ax_strength = fig.add_subplot(layout[2, 2], sharey=ax_heat)
    max_abs = np.abs(rho).max(axis=1)[row_order][:, None]
    ax_strength.imshow(
        max_abs,
        aspect="auto",
        interpolation="nearest",
        cmap="Greys",
        vmin=0,
        vmax=max(0.5, float(max_abs.max())),
        origin="upper",
    )
    ax_strength.set_xticks([0])
    ax_strength.set_xticklabels(["max |rho|"], rotation=90)
    ax_strength.xaxis.tick_bottom()
    ax_strength.set_yticks([])
    for spine in ax_strength.spines.values():
        spine.set_linewidth(0.25)

    cax = fig.add_subplot(layout[2, 3])
    colorbar = fig.colorbar(image, cax=cax)
    colorbar.set_label("Row-wise z score of Spearman rho")
    colorbar.outline.set_linewidth(0.25)
    colorbar.ax.tick_params(width=0.25, length=2)

    fig.suptitle(
        "Exploratory AEF--Pfam correlation-profile landscape; "
        "average linkage, correlation distance",
        x=0.52,
        y=0.998,
        fontsize=6,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.06, right=0.985, top=0.975, bottom=0.06)
    fig.savefig(paths["pdf"], transparent=True, edgecolor="none", bbox_inches="tight")
    fig.savefig(paths["svg"], transparent=True, edgecolor="none", bbox_inches="tight")
    plt.close(fig)
    return paths


def output_hashes(out_dir: Path, exclude: set[Path] | None = None) -> dict[str, str]:
    excluded = {path.resolve() for path in (exclude or set())}
    hashes: dict[str, str] = {}
    for path in sorted(out_dir.iterdir()):
        if path.is_file() and path.resolve() not in excluded:
            hashes[path.name] = sha256(path)
    return hashes


def main() -> None:
    args = parse_args()
    started_epoch = time.time()
    started_utc = datetime.now(timezone.utc)
    run_id = started_utc.strftime("%Y%m%d_%H%M%S")
    out_dir = args.output_parent.resolve() / f"full_aef_corrected_run_{run_id}"
    if out_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing run directory: {out_dir}")
    out_dir.mkdir(parents=True)

    input_paths = {
        "reconciled_manifest": args.manifest.resolve(),
        "reconstructed_raw_pfam_counts": args.pfam_counts.resolve(),
        "aef_embeddings": args.aef_embeddings.resolve(),
    }
    manifest, counts, aef, pfams = load_validated_inputs(
        input_paths["reconciled_manifest"],
        input_paths["reconstructed_raw_pfam_counts"],
        input_paths["aef_embeddings"],
    )
    axes = EXPECTED_AXES
    count_values = counts[pfams].to_numpy(float)
    aef_values = aef[axes].to_numpy(float)
    prevalence = (count_values > 0).sum(axis=0)

    join_audit = manifest[
        [
            "master_row",
            "Genome",
            "Species",
            "Phylum",
            "DD latitude",
            "DD longitude",
            "aef_id_number",
            "safe_for_raw_pfam_analysis",
            "safe_for_aef_pfam_analysis",
        ]
    ].copy()
    join_audit["pfam_row_exact_id"] = counts["Genome"]
    join_audit["aef_row_exact_id"] = aef["Genome"]
    join_path = out_dir / "exact_id_join_audit.csv"
    join_audit.to_csv(join_path, index=False)

    pooled = compute_screen(count_values, aef_values, pfams, axes)
    pooled_paths = write_screen_outputs(
        "pooled_126", out_dir, pooled, pfams, axes, prevalence
    )
    sentinel = validate_sentinels(count_values, aef_values, pfams, axes, pooled)
    sentinel_path = out_dir / "scipy_sentinel_validation.csv"
    sentinel.to_csv(sentinel_path, index=False, float_format="%.17g")

    subgroup_summaries: list[dict[str, object]] = []
    phylum = manifest["Phylum"].astype(str).to_numpy()
    phylum_centered_values = count_values.copy()
    for phylum_name in sorted(set(phylum)):
        group_mask = phylum == phylum_name
        phylum_centered_values[group_mask] -= phylum_centered_values[group_mask].mean(
            axis=0, keepdims=True
        )
    centered_eligible = np.ptp(phylum_centered_values, axis=0) > 0
    centered_pfams = [pfam for pfam, keep in zip(pfams, centered_eligible) if keep]
    centered_result = compute_screen(
        phylum_centered_values[:, centered_eligible],
        aef_values,
        centered_pfams,
        axes,
    )
    write_screen_outputs(
        "phylum_centered_126",
        out_dir,
        centered_result,
        centered_pfams,
        axes,
        prevalence[centered_eligible],
    )
    centered_summary = {
        "n_genomes": len(manifest),
        "n_pfams": len(centered_pfams),
        "n_tests": int(np.asarray(centered_result["p"]).size),
        "n_p_lt_0.001": int((np.asarray(centered_result["p"]) < 0.001).sum()),
        "n_q_lt_0.05": int((np.asarray(centered_result["q"]) < 0.05).sum()),
        "minimum_q": float(np.asarray(centered_result["q"]).min()),
    }

    subgroup_specs = [(name, 3) for name in sorted(set(phylum))]
    subgroup_specs.append(("Rhodophyta", 10))
    for phylum_name, minimum_nonzero in subgroup_specs:
        subgroup_mask = phylum == phylum_name
        subgroup_counts = count_values[subgroup_mask]
        subgroup_aef = aef_values[subgroup_mask]
        eligible = ((subgroup_counts > 0).sum(axis=0) >= minimum_nonzero) & (
            np.ptp(subgroup_counts, axis=0) > 0
        )
        subgroup_pfams = [pfam for pfam, keep in zip(pfams, eligible) if keep]
        subgroup_result = compute_screen(
            subgroup_counts[:, eligible], subgroup_aef, subgroup_pfams, axes
        )
        subgroup_prevalence = (subgroup_counts[:, eligible] > 0).sum(axis=0)
        prefix = (
            f"{phylum_name.lower()}_n{int(subgroup_mask.sum())}_"
            f"min{minimum_nonzero}_nonzero"
        )
        write_screen_outputs(
            prefix,
            out_dir,
            subgroup_result,
            subgroup_pfams,
            axes,
            subgroup_prevalence,
        )
        subgroup_summaries.append(
            {
                "phylum": phylum_name,
                "n_genomes": int(subgroup_mask.sum()),
                "minimum_nonzero_genomes": minimum_nonzero,
                "n_pfams": len(subgroup_pfams),
                "n_tests": len(subgroup_pfams) * len(axes),
                "n_p_lt_0.001": int((np.asarray(subgroup_result["p"]) < 0.001).sum()),
                "n_q_lt_0.05": int((np.asarray(subgroup_result["q"]) < 0.05).sum()),
                "minimum_q": float(np.asarray(subgroup_result["q"]).min()),
            }
        )
    subgroup_path = out_dir / "within_phylum_specification_summary.csv"
    pd.DataFrame(subgroup_summaries).to_csv(
        subgroup_path, index=False, float_format="%.17g"
    )

    figure_paths: dict[str, Path] = {}
    if not args.skip_figure:
        figure_paths = cluster_and_plot(out_dir, pfams, axes, pooled)

    pooled_rho = np.asarray(pooled["rho"])
    pooled_p = np.asarray(pooled["p"])
    pooled_q = np.asarray(pooled["q"])
    if figure_paths:
        membership = pd.read_csv(figure_paths["membership"])
        unique_profile_count = int(membership["is_representative"].sum())
    else:
        unique_profile_count = None

    summary = {
        "cohort_size": len(manifest),
        "phylum_counts": manifest["Phylum"].value_counts().sort_index().to_dict(),
        "strict_variable_pfam_profiles": len(pfams),
        "aef_axes": len(axes),
        "global_tests": int(pooled_p.size),
        "p_lt_0.05": int((pooled_p < 0.05).sum()),
        "p_lt_0.001": int((pooled_p < 0.001).sum()),
        "bonferroni_p_lt_0.05": int((pooled_p < (0.05 / pooled_p.size)).sum()),
        "bh_q_lt_0.05": int((pooled_q < 0.05).sum()),
        "minimum_q": float(pooled_q.min()),
        "rho_minimum": float(pooled_rho.min()),
        "rho_maximum": float(pooled_rho.max()),
        "displayed_pfam_profiles": len(pfams) if figure_paths else None,
        "unique_correlation_profiles": unique_profile_count,
        "exact_duplicate_profiles": (
            len(pfams) - unique_profile_count if unique_profile_count is not None else None
        ),
        "phylum_centered_screen": centered_summary,
        "within_phylum_screens": subgroup_summaries,
    }
    summary_path = out_dir / "computed_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    completed_utc = datetime.now(timezone.utc)
    manifest_path = out_dir / "run_manifest.json"
    software = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "statsmodels": statsmodels.__version__,
        "matplotlib": mpl.__version__,
    }
    run_manifest = {
        "script_version": SCRIPT_VERSION,
        "script": {
            "path_at_runtime": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "command": [str(value) for value in sys.argv],
        "started_utc": started_utc.isoformat(),
        "completed_utc": completed_utc.isoformat(),
        "elapsed_seconds": completed_utc.timestamp() - started_epoch,
        "platform": platform.platform(),
        "software": software,
        "inputs": {
            name: {"path_at_runtime": str(path), "sha256": sha256(path)}
            for name, path in input_paths.items()
        },
        "methods": {
            "join": "exact Genome ID; safe_for_aef_pfam_analysis records; manifest order",
            "pfam_columns": "strict regular expression ^PF\\d{5}$; nonnegative integer counts; variable across cohort",
            "aef_columns": "exact A00--A63",
            "association": "Spearman rank correlation with average tie ranks",
            "p_values": "two-sided asymptotic Student t transformation; df=n-2",
            "multiple_testing": "one global Benjamini--Hochberg family for each stated screen",
            "display_transform": "row-wise z score of raw rho across 64 axes; sample standard deviation ddof=1",
            "clustering": "average-linkage hierarchical clustering with correlation distance; rows and columns clustered once",
            "duplicate_handling": "all Pfams retained for inference and display; exact duplicate profiles documented in a membership audit",
        },
        "computed_summary": summary,
        "peak_rss_raw": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "output_files_sha256": output_hashes(out_dir, exclude={manifest_path}),
    }
    manifest_path.write_text(json.dumps(run_manifest, indent=2, sort_keys=True) + "\n")

    print(json.dumps({"output_directory": str(out_dir), **summary}, indent=2))


if __name__ == "__main__":
    main()

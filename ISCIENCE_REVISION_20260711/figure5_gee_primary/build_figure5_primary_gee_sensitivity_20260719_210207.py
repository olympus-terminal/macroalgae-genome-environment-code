#!/usr/bin/env python3
"""Build the primary GEE selected-family sensitivity figure from real outputs.

The figure displays every one of the 84 exact-ID GEE pairs selected by the
global discovery-wide Benjamini--Hochberg threshold.  Scientific values are
read directly from the authenticated sensitivity and refined structured-null
tables.  No values are simulated, reconstructed from summaries, or hardcoded.

Created: 2026-07-19.
"""

from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


SCRIPT = Path(__file__).resolve()
ROOT = SCRIPT.parents[2]
REVISION = ROOT / "ISCIENCE_REVISION_20260711"
GEE = REVISION / "gee_validation"
OUTDIR = SCRIPT.parent
RUN_ID = "20260719_214228"

SENSITIVITY = GEE / "GEE_primary_selected84_sensitivity_20260719_205317.csv"
REFINED_NULL = GEE / "GEE_primary_selected84_structured_null_refined99999_20260719_210049.csv"
REFINED_CANDIDATES = GEE / "GEE_primary_selected84_candidate_summary_refined99999_20260719_210049.csv"
REFINED_MANIFEST = GEE / "GEE_primary_selected84_manifest_refined99999_20260719_210049.json"

PDF = OUTDIR / f"Figure5_primary_GEE_sensitivity_{RUN_ID}.pdf"
SVG = OUTDIR / f"Figure5_primary_GEE_sensitivity_{RUN_ID}.svg"
INTEGRITY = OUTDIR / f"Figure5_primary_GEE_sensitivity_integrity_{RUN_ID}.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_new(paths: list[Path]) -> None:
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError(f"Refusing to overwrite existing outputs: {existing}")


ENVIRONMENT_ORDER = [
    "sst_mean_c",
    "sst_max_c",
    "sst_min_c",
    "sst_summer_c",
    "sst_winter_c",
    "depth_meters",
    "distance_coast_km",
]

ENVIRONMENT_LABELS = {
    "sst_mean_c": "SST mean",
    "sst_max_c": "SST maximum",
    "sst_min_c": "SST minimum",
    "sst_summer_c": "SST summer",
    "sst_winter_c": "SST winter",
    "depth_meters": "Depth",
    "distance_coast_km": "Distance to coast",
}

CHECKS = [
    ("raw_site", "raw_count", "site_mean_spearman", "Raw\nsite"),
    ("total_site", "per_total_pfam_hit", "site_mean_spearman", "Total-hit\nsite"),
    ("peptide_site", "per_final_peptide_record", "site_mean_spearman", "Peptide\nsite"),
    ("quality", "per_total_pfam_hit", "quality_phylum_sitecluster", "Quality +\nphylum"),
    ("busco50", "per_total_pfam_hit", "site_mean_busco_ge50_spearman", "BUSCO >=50%\nsite"),
    ("topology", "per_total_pfam_hit", "phylo_pev3_quality_sitecluster", "Topology\nPEV3"),
]


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    for path in (SENSITIVITY, REFINED_NULL, REFINED_CANDIDATES, REFINED_MANIFEST):
        if not path.is_file():
            raise FileNotFoundError(path)
    sensitivity = pd.read_csv(SENSITIVITY, low_memory=False)
    null = pd.read_csv(REFINED_NULL, low_memory=False)
    candidates = pd.read_csv(REFINED_CANDIDATES, low_memory=False)
    manifest = json.loads(REFINED_MANIFEST.read_text(encoding="utf-8"))
    if len(candidates) != 84 or len(null) != 84:
        raise RuntimeError("Expected 84 candidate and 84 structured-null rows")
    if candidates[["pfam", "environment"]].duplicated().any():
        raise RuntimeError("Candidate table contains duplicate Pfam--environment pairs")
    if int(candidates["robust_candidate_all_required_checks"].sum()) != 49:
        raise RuntimeError("Refined candidate result is not the audited 49/84 outcome")
    if int(candidates["discovery_bonferroni_significant"].sum()) != 13:
        raise RuntimeError("Expected 13 discovery Bonferroni pairs")
    bonf = candidates[candidates["discovery_bonferroni_significant"]]
    if not bonf["robust_candidate_all_required_checks"].all():
        raise RuntimeError("Not every discovery Bonferroni pair retained all seven checks")
    return sensitivity, null, candidates, manifest


def build_matrix(
    sensitivity: pd.DataFrame,
    null: pd.DataFrame,
    candidates: pd.DataFrame,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    order_map = {name: index for index, name in enumerate(ENVIRONMENT_ORDER)}
    ordered = candidates.assign(
        _environment_order=candidates["environment"].map(order_map),
    ).sort_values(
        ["_environment_order", "discovery_global_bh_q", "pfam"],
        kind="mergesort",
    ).reset_index(drop=True)
    if ordered["_environment_order"].isna().any():
        unknown = ordered.loc[ordered["_environment_order"].isna(), "environment"].unique()
        raise RuntimeError(f"Unexpected environment variables: {unknown.tolist()}")

    effects = np.full((len(ordered), 7), np.nan, dtype=float)
    supported = np.zeros((len(ordered), 7), dtype=bool)
    for row_index, pair in ordered.iterrows():
        pair_rows = sensitivity[
            sensitivity["pfam"].eq(pair["pfam"])
            & sensitivity["environment"].eq(pair["environment"])
        ]
        for column_index, (_, representation, method, _) in enumerate(CHECKS):
            match = pair_rows[
                pair_rows["representation"].eq(representation)
                & pair_rows["method"].eq(method)
            ]
            if len(match) != 1:
                raise RuntimeError(
                    f"Expected one sensitivity row for {pair['pfam']}--{pair['environment']} "
                    f"{representation}/{method}; found {len(match)}"
                )
            record = match.iloc[0]
            effects[row_index, column_index] = float(record["effect"])
            supported[row_index, column_index] = bool(
                np.isfinite(record["selected_GEE_set_bh_q"])
                and float(record["selected_GEE_set_bh_q"]) < 0.05
                and np.sign(float(record["effect"])) == np.sign(float(pair["discovery_spearman_r"]))
            )

        null_match = null[
            null["pfam"].eq(pair["pfam"])
            & null["environment"].eq(pair["environment"])
        ]
        if len(null_match) != 1:
            raise RuntimeError(
                f"Expected one refined-null row for {pair['pfam']}--{pair['environment']}"
            )
        record = null_match.iloc[0]
        effects[row_index, 6] = float(record["observed_spearman_rho"])
        supported[row_index, 6] = bool(
            float(record["selected_GEE_empirical_bh_q"]) < 0.05
            and np.sign(float(record["observed_spearman_rho"]))
            == np.sign(float(pair["discovery_spearman_r"]))
        )

    expected = ordered["robust_candidate_all_required_checks"].to_numpy(bool)
    observed = supported.all(axis=1)
    if not np.array_equal(expected, observed):
        bad = ordered.loc[expected != observed, ["pfam", "environment"]]
        raise RuntimeError(f"Figure check matrix disagrees with candidate summary: {bad.to_dict('records')}")
    return ordered, effects, supported


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
        }
    )


def make_figure(ordered: pd.DataFrame, effects: np.ndarray, supported: np.ndarray) -> None:
    configure_matplotlib()
    colors = [
        (0.0, 0.3, 0.7),
        (0.3, 0.5, 0.9),
        (0.7, 0.8, 0.95),
        (0.95, 0.95, 0.95),
        (0.95, 0.8, 0.7),
        (0.9, 0.4, 0.3),
        (0.7, 0.1, 0.1),
    ]
    cmap = LinearSegmentedColormap.from_list("protocol_diverging", colors, N=256)

    fig = plt.figure(figsize=(7.0, 11.5), constrained_layout=False)
    grid = fig.add_gridspec(
        3,
        4,
        height_ratios=[1.05, 0.18, 9.15],
        width_ratios=[5.45, 0.22, 0.32, 0.16],
        hspace=0.13,
        wspace=0.07,
        left=0.205,
        right=0.965,
        top=0.985,
        bottom=0.052,
    )

    # Panel A: selected and retained counts by named GEE variable.
    ax_counts = fig.add_subplot(grid[0, 0])
    selected_counts = ordered.groupby("environment", sort=False).size().reindex(ENVIRONMENT_ORDER)
    retained_counts = (
        ordered.groupby("environment", sort=False)["robust_candidate_all_required_checks"]
        .sum()
        .reindex(ENVIRONMENT_ORDER)
    )
    x = np.arange(len(ENVIRONMENT_ORDER))
    ax_counts.bar(x, selected_counts, width=0.72, color="#d8e1ea", edgecolor="#5b7083", linewidth=0.25, label="Global-BH discovery pairs")
    ax_counts.bar(x, retained_counts, width=0.46, color="#244f73", edgecolor="#244f73", linewidth=0.25, label="Retained in all seven specifications")
    for index, (selected, retained) in enumerate(zip(selected_counts, retained_counts)):
        ax_counts.text(index, selected + 0.35, f"{int(retained)}/{int(selected)}", ha="center", va="bottom")
    ax_counts.set_xticks(x)
    ax_counts.set_xticklabels([ENVIRONMENT_LABELS[name] for name in ENVIRONMENT_ORDER], rotation=24, ha="right")
    ax_counts.set_ylabel("Pairs")
    ax_counts.set_ylim(0, max(selected_counts) + 4)
    ax_counts.spines[["top", "right"]].set_visible(False)
    ax_counts.grid(axis="y", color="#b7b7b7", linewidth=0.25, alpha=0.7)
    ax_counts.set_axisbelow(True)
    ax_counts.legend(frameon=False, loc="upper right", ncol=2, handlelength=1.1, columnspacing=1.2)
    ax_counts.set_title("A  Primary exact-ID GEE associations after seven sensitivity specifications", loc="left", fontweight="bold")

    colorbar_axis = fig.add_subplot(grid[1, 0])

    # Panel B: all 84 pairs and all seven required checks.
    ax_heat = fig.add_subplot(grid[2, 0])
    image = ax_heat.imshow(effects, aspect="auto", interpolation="nearest", cmap=cmap, vmin=-0.65, vmax=0.65)
    row_labels = [
        f"{row.pfam} | {ENVIRONMENT_LABELS[row.environment]}"
        for row in ordered.itertuples(index=False)
    ]
    ax_heat.set_yticks(np.arange(len(ordered)))
    ax_heat.set_yticklabels(row_labels)
    ax_heat.set_xticks(np.arange(7))
    ax_heat.set_xticklabels([entry[3] for entry in CHECKS] + ["Structured\nnull"])
    ax_heat.xaxis.tick_top()
    ax_heat.tick_params(axis="x", pad=2)
    ax_heat.tick_params(axis="y", length=0)
    ax_heat.set_title("B  Effect estimates; dots mark selected-family q < 0.05 with discovery-consistent direction", loc="left", fontweight="bold", pad=19)
    dot_y, dot_x = np.where(supported)
    ax_heat.scatter(dot_x, dot_y, s=3.2, c="black", marker="o", linewidths=0)

    # Separators and group labels are derived from the ordered real table.
    starts: list[tuple[int, str]] = []
    cursor = 0
    for environment in ENVIRONMENT_ORDER:
        count = int((ordered["environment"] == environment).sum())
        starts.append((cursor, environment))
        cursor += count
        if cursor < len(ordered):
            ax_heat.axhline(cursor - 0.5, color="black", linewidth=0.25)

    # Thin discovery-effect track.
    ax_discovery = fig.add_subplot(grid[2, 1], sharey=ax_heat)
    discovery = ordered["discovery_spearman_r"].to_numpy(float)[:, None]
    ax_discovery.imshow(discovery, aspect="auto", interpolation="nearest", cmap=cmap, vmin=-0.65, vmax=0.65)
    ax_discovery.set_xticks([0])
    ax_discovery.set_xticklabels(["Discovery\nr"])
    ax_discovery.xaxis.tick_top()
    ax_discovery.tick_params(axis="y", left=False, labelleft=False)
    ax_discovery.tick_params(axis="x", length=0, pad=2)
    for start, environment in starts[1:]:
        ax_discovery.axhline(start - 0.5, color="black", linewidth=0.25)

    # Check-count track, with all-seven pairs visually distinct.
    ax_checks = fig.add_subplot(grid[2, 2], sharey=ax_heat)
    check_counts = ordered["checks_passed"].to_numpy(int)
    check_color = np.where(check_counts == 7, "#1b7837", "#b8b8b8")
    ax_checks.barh(np.arange(len(ordered)), check_counts, height=0.82, color=check_color, edgecolor="none")
    ax_checks.set_xlim(0, 7.2)
    ax_checks.set_xticks([0, 7])
    ax_checks.set_xticklabels(["0", "7"])
    ax_checks.xaxis.tick_top()
    ax_checks.tick_params(axis="y", left=False, labelleft=False)
    ax_checks.tick_params(axis="x", length=0, pad=2)
    ax_checks.set_title("Checks", pad=19)
    ax_checks.spines[["left", "right", "bottom"]].set_visible(False)
    for start, environment in starts[1:]:
        ax_checks.axhline(start - 0.5, color="black", linewidth=0.25)

    # Bonferroni discovery marker track.
    ax_bonf = fig.add_subplot(grid[2, 3], sharey=ax_heat)
    bonf_y = np.flatnonzero(ordered["discovery_bonferroni_significant"].to_numpy(bool))
    ax_bonf.scatter(np.zeros(len(bonf_y)), bonf_y, s=7, marker="D", color="#54278f", linewidths=0)
    ax_bonf.set_xlim(-0.6, 0.6)
    ax_bonf.set_xticks([0])
    ax_bonf.set_xticklabels(["Bonf."])
    ax_bonf.xaxis.tick_top()
    ax_bonf.tick_params(axis="x", length=0, pad=2)
    ax_bonf.tick_params(axis="y", left=False, labelleft=False)
    ax_bonf.spines[["left", "right", "bottom"]].set_visible(False)
    for start, environment in starts[1:]:
        ax_bonf.axhline(start - 0.5, color="black", linewidth=0.25)

    colorbar = fig.colorbar(image, cax=colorbar_axis, orientation="horizontal")
    colorbar.set_label("Rank-association effect estimate")
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.outline.set_linewidth(0.25)
    colorbar.ax.tick_params(width=0.25, length=2)

    fig.text(
        0.205,
        0.009,
        "All 84 global-BH discovery pairs are shown. Required checks: unique-site raw, total-hit-normalized, and peptide-normalized Spearman;\n"
        "total-hit quality/phylum site-cluster model; BUSCO >=50% unique-site Spearman; total-hit PEV3 topology/site-cluster model;\n"
        "and total-hit site-phylum-composition structured null (99,999 permutations).",
        ha="left",
        va="bottom",
    )

    fig.savefig(PDF, format="pdf", transparent=True, edgecolor="none")
    fig.savefig(SVG, format="svg", transparent=True, edgecolor="none")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    require_new([PDF, SVG, INTEGRITY])
    sensitivity, null, candidates, manifest = load_inputs()
    ordered, effects, supported = build_matrix(sensitivity, null, candidates)
    make_figure(ordered, effects, supported)
    if not PDF.is_file() or not SVG.is_file():
        raise RuntimeError("Figure export failed")

    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "generator": {"path": str(SCRIPT), "sha256": sha256(SCRIPT)},
        "inputs": {
            str(path): {"sha256": sha256(path), "bytes": path.stat().st_size}
            for path in (SENSITIVITY, REFINED_NULL, REFINED_CANDIDATES, REFINED_MANIFEST)
        },
        "results": {
            "selected_pairs": len(ordered),
            "pairs_retaining_all_seven": int(ordered["robust_candidate_all_required_checks"].sum()),
            "discovery_bonferroni_pairs": int(ordered["discovery_bonferroni_significant"].sum()),
            "bonferroni_pairs_retaining_all_seven": int(
                ordered.loc[ordered["discovery_bonferroni_significant"], "robust_candidate_all_required_checks"].sum()
            ),
            "matrix_shape": list(effects.shape),
            "supported_cells": int(supported.sum()),
            "variable_counts": {
                name: {
                    "selected": int((ordered["environment"] == name).sum()),
                    "all_seven": int(
                        ordered.loc[ordered["environment"] == name, "robust_candidate_all_required_checks"].sum()
                    ),
                }
                for name in ENVIRONMENT_ORDER
            },
        },
        "outputs": {
            str(path): {"sha256": sha256(path), "bytes": path.stat().st_size}
            for path in (PDF, SVG)
        },
        "rendering": {
            "font": "Arial",
            "font_size_pt": 6,
            "line_width_pt": 0.25,
            "pdf_fonttype": 42,
            "svg_fonttype": "none",
            "transparent_background": True,
            "vector_only": True,
        },
        "data_integrity": {
            "real_data_only": True,
            "synthetic_scientific_values": False,
            "randomness_used_in_figure_build": False,
            "refined_null_permutations": manifest["parameters"]["permutations"],
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": mpl.__version__,
        },
        "result": "PASS",
    }
    INTEGRITY.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"pdf": str(PDF), "svg": str(SVG), "integrity": str(INTEGRITY), "all_seven": 49}, indent=2))


if __name__ == "__main__":
    main()

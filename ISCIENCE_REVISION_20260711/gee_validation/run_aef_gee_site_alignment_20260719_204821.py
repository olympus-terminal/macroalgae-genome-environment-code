#!/usr/bin/env python3
"""Quantify site-level alignment between AEF axes and named GEE variables.

This descriptive cross-representation alignment compares a fixed latent
representation with named environmental variables at the same recorded sites;
axis meanings remain distributed and nonunique. All 64 AEF axes are tested
against all 13 named GEE variables after exact-ID joining and aggregation to the
90 unique observed coordinate sites. Missing GEE source pixels are retained as
missing; no values are imputed or simulated.

Created: 2026-07-19.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats
import statsmodels
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "ISCIENCE_REVISION_20260711" / "gee_validation"
SCRIPT = Path(__file__).resolve()
AEF = ROOT / "AlphaEarth" / "CSV" / "alphaearth_embeddings_20251019_122918.csv"
GEE = OUTDIR / "exact_id_gee_environmental_extraction_20260712_071838.csv"
GEE_MANIFEST = OUTDIR / "exact_id_gee_environmental_extraction_manifest_20260712_071838.json"

EXPECTED = {
    AEF: "e0b05e727aec4a5b45565c9026de21c46561a18c06786f235eae371771a4cf87",
    GEE: "6d7c1e41eb5651c34464c41c6fadd637db241063d3dfe3c1cc1e614bd6e40418",
}

AXES = [f"A{i:02d}" for i in range(64)]
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
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def authenticate() -> dict[str, dict[str, object]]:
    records = {}
    for path, expected in EXPECTED.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        observed = sha256(path)
        if observed != expected:
            raise RuntimeError(f"Authenticated input changed: {path}: {observed}")
        records[path.name] = {
            "path": str(path.relative_to(ROOT)),
            "realpath": str(path.resolve()),
            "sha256": observed,
            "bytes": path.stat().st_size,
        }
    if not GEE_MANIFEST.is_file():
        raise FileNotFoundError(GEE_MANIFEST)
    gee_manifest = json.loads(GEE_MANIFEST.read_text(encoding="utf-8"))
    if gee_manifest["output"]["sha256"] != EXPECTED[GEE]:
        raise RuntimeError("GEE provenance manifest does not authenticate extraction")
    records[GEE_MANIFEST.name] = {
        "path": str(GEE_MANIFEST.relative_to(ROOT)),
        "realpath": str(GEE_MANIFEST.resolve()),
        "sha256": sha256(GEE_MANIFEST),
        "bytes": GEE_MANIFEST.stat().st_size,
    }
    return records


def load_site_table() -> tuple[pd.DataFrame, dict[str, object]]:
    aef = pd.read_csv(AEF, low_memory=False)
    gee = pd.read_csv(GEE, low_memory=False)
    if len(aef) != 126 or aef["Genome"].nunique() != 126:
        raise RuntimeError("AEF input is not 126 unique genome IDs")
    if len(gee) != 126 or gee["genome_id"].nunique() != 126:
        raise RuntimeError("GEE input is not 126 unique genome IDs")
    if set(aef["Genome"]) != set(gee["genome_id"]):
        raise RuntimeError("AEF and GEE exact genome-ID sets differ")
    observed_axes = [column for column in aef.columns if column.startswith("A") and len(column) == 3]
    if observed_axes != AXES:
        raise RuntimeError(f"Expected A00--A63, observed {observed_axes}")

    merged = gee.merge(aef, left_on="genome_id", right_on="Genome", validate="one_to_one")
    for left, right in (("latitude", "DD latitude"), ("longitude", "DD longitude")):
        difference = (
            pd.to_numeric(merged[left], errors="raise")
            - pd.to_numeric(merged[right], errors="raise")
        ).abs()
        if float(difference.max()) > 1e-10:
            raise RuntimeError(f"Exact-ID AEF and GEE coordinates differ for {left}")
    merged["site_id"] = (
        merged["latitude"].round(10).map(lambda value: f"{value:.10f}")
        + "|"
        + merged["longitude"].round(10).map(lambda value: f"{value:.10f}")
    )
    grouped = merged.groupby("site_id", sort=True)
    if grouped.ngroups != 90:
        raise RuntimeError(f"Expected 90 unique coordinate sites, found {grouped.ngroups}")
    aef_variation = grouped[AXES].nunique(dropna=False).gt(1).any(axis=1)
    gee_variation = grouped[ENVIRONMENTAL_COLUMNS].nunique(dropna=False).gt(1).any(axis=1)
    if aef_variation.any() or gee_variation.any():
        raise RuntimeError("AEF or GEE values differ within an exact coordinate site")

    site = grouped[AXES].first().join(grouped[ENVIRONMENTAL_COLUMNS].first())
    site["n_genomes"] = grouped.size()
    audit = {
        "exact_id_genomes": len(merged),
        "unique_coordinate_sites": len(site),
        "repeated_coordinate_sites": int((site["n_genomes"] > 1).sum()),
        "maximum_genomes_per_site": int(site["n_genomes"].max()),
        "within_site_aef_disagreements": int(aef_variation.sum()),
        "within_site_gee_disagreements": int(gee_variation.sum()),
    }
    return site, audit


def calculate(site: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    for axis in AXES:
        for environment in ENVIRONMENTAL_COLUMNS:
            pair = site[[axis, environment]].dropna()
            if len(pair) < 10 or pair[axis].nunique() < 2 or pair[environment].nunique() < 2:
                rho = np.nan
                p_value = np.nan
            else:
                result = stats.spearmanr(pair[axis].to_numpy(float), pair[environment].to_numpy(float))
                rho = float(result.statistic)
                p_value = float(result.pvalue)
            rows.append({
                "aef_axis": axis,
                "gee_variable": environment,
                "n_unique_sites": len(pair),
                "spearman_rho": rho,
                "p_value_two_sided": p_value,
            })
    correlations = pd.DataFrame(rows)
    valid = correlations["p_value_two_sided"].notna()
    if int(valid.sum()) != 64 * 13:
        raise RuntimeError("Not all 832 site-level AEF--GEE pairs were estimable")
    correlations.loc[valid, "global_bh_q_832"] = multipletests(
        correlations.loc[valid, "p_value_two_sided"].to_numpy(float), method="fdr_bh"
    )[1]
    correlations.loc[valid, "global_bonferroni_p_832"] = np.minimum(
        correlations.loc[valid, "p_value_two_sided"].to_numpy(float) * int(valid.sum()), 1.0
    )
    correlations["global_bh_q_lt_0.05"] = correlations["global_bh_q_832"] < 0.05
    correlations["global_bonferroni_p_lt_0.05"] = correlations["global_bonferroni_p_832"] < 0.05
    correlations["interpretation_boundary"] = (
        "site-level external alignment; not a semantic decoder or a unique physical label for an AEF axis"
    )
    correlations = correlations.sort_values(
        ["p_value_two_sided", "aef_axis", "gee_variable"], kind="stable"
    ).reset_index(drop=True)

    variable_summary = (
        correlations.groupby("gee_variable", sort=False)
        .agg(
            n_sites_min=("n_unique_sites", "min"),
            n_sites_max=("n_unique_sites", "max"),
            n_axes_global_bh_q_lt_0_05=("global_bh_q_lt_0.05", "sum"),
            n_axes_global_bonferroni_p_lt_0_05=("global_bonferroni_p_lt_0.05", "sum"),
        )
        .reset_index()
    )
    strongest_by_variable = correlations.loc[
        correlations.groupby("gee_variable")["spearman_rho"].apply(lambda values: values.abs().idxmax())
    ][["gee_variable", "aef_axis", "spearman_rho", "global_bh_q_832"]]
    variable_summary = variable_summary.merge(
        strongest_by_variable.rename(columns={
            "aef_axis": "strongest_axis",
            "spearman_rho": "strongest_axis_rho",
            "global_bh_q_832": "strongest_axis_global_bh_q",
        }),
        on="gee_variable", how="left", validate="one_to_one",
    )

    strongest_by_axis = correlations.loc[
        correlations.groupby("aef_axis")["spearman_rho"].apply(lambda values: values.abs().idxmax())
    ][["aef_axis", "gee_variable", "n_unique_sites", "spearman_rho", "global_bh_q_832"]]
    axis_summary = strongest_by_axis.rename(columns={
        "gee_variable": "strongest_aligned_gee_variable",
        "spearman_rho": "strongest_alignment_rho",
        "global_bh_q_832": "strongest_alignment_global_bh_q",
    }).sort_values("strongest_alignment_rho", key=lambda values: values.abs(), ascending=False)
    return correlations, variable_summary, axis_summary


def main() -> int:
    started = datetime.now(timezone.utc)
    clock = time.monotonic()
    inputs = authenticate()
    site, site_audit = load_site_table()
    correlations, variable_summary, axis_summary = calculate(site)
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    outputs = {
        "correlations": OUTDIR / f"AEF_GEE_site_alignment_832pairs_{run_id}.csv",
        "variable_summary": OUTDIR / f"AEF_GEE_site_alignment_variable_summary_{run_id}.csv",
        "axis_summary": OUTDIR / f"AEF_GEE_site_alignment_axis_summary_{run_id}.csv",
        "manifest": OUTDIR / f"AEF_GEE_site_alignment_manifest_{run_id}.json",
    }
    for path in outputs.values():
        if path.exists():
            raise FileExistsError(path)
    correlations.to_csv(outputs["correlations"], index=False, lineterminator="\n")
    variable_summary.to_csv(outputs["variable_summary"], index=False, lineterminator="\n")
    axis_summary.to_csv(outputs["axis_summary"], index=False, lineterminator="\n")

    for label in ("correlations", "variable_summary", "axis_summary"):
        if not outputs[label].is_file() or outputs[label].stat().st_size == 0:
            raise RuntimeError(f"Output was not materialized: {outputs[label]}")
    reread = pd.read_csv(outputs["correlations"])
    if len(reread) != 832 or int(reread["global_bh_q_lt_0.05"].sum()) != int(correlations["global_bh_q_lt_0.05"].sum()):
        raise RuntimeError("Materialized correlation output failed row/count audit")

    manifest = {
        "purpose": "site-level external alignment of 64 unitless AEF axes with 13 named GEE variables",
        "interpretation_boundary": "not a semantic decoder; no AEF axis is assigned a unique physical meaning",
        "generated_by": str(SCRIPT),
        "generator_sha256": sha256(SCRIPT),
        "started_utc": started.isoformat(),
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "runtime_seconds": time.monotonic() - clock,
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "statsmodels": statsmodels.__version__,
        },
        "inputs": inputs,
        "site_audit": site_audit,
        "results": {
            "tests": len(correlations),
            "global_bh_q_lt_0.05": int(correlations["global_bh_q_lt_0.05"].sum()),
            "global_bonferroni_p_lt_0.05": int(correlations["global_bonferroni_p_lt_0.05"].sum()),
            "a52": correlations[correlations["aef_axis"].eq("A52")].to_dict("records"),
            "a53": correlations[correlations["aef_axis"].eq("A53")].to_dict("records"),
        },
        "outputs": {
            label: {
                "path": str(path.relative_to(ROOT)),
                "realpath": str(path.resolve()),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
            for label, path in outputs.items()
            if label != "manifest"
        },
    }
    outputs["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "run_id": run_id,
        "unique_sites": site_audit["unique_coordinate_sites"],
        "tests": len(correlations),
        "global_bh_q_lt_0.05": manifest["results"]["global_bh_q_lt_0.05"],
        "global_bonferroni_p_lt_0.05": manifest["results"]["global_bonferroni_p_lt_0.05"],
        "outputs": {label: str(path.resolve()) for label, path in outputs.items()},
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

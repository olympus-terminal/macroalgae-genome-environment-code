#!/usr/bin/env python3
"""Run selected-family sensitivity analyses for the primary exact-ID GEE results.

The selected family is defined without manual candidate picking: every one of
the 84 Pfam--GEE pairs that met global Benjamini--Hochberg q < 0.05 in the
authenticated 87,123-test exact-ID discovery screen is carried forward.

All calculations use the observed 126-genome data. Random-number generation is
used only for paired bootstrap resampling of observed coordinate sites and for
structured permutation of intact observed site labels within site
phylum-composition strata. No synthetic scientific values are generated.

Required seven-check rule, matched to the existing selected AEF workflow:
1. raw-count unique-site Spearman;
2. total-Pfam-hit-normalized unique-site Spearman;
3. peptide-record-normalized unique-site Spearman;
4. total-hit coordinate-clustered quality/phylum model;
5. total-hit BUSCO >=50% unique-site Spearman;
6. total-hit primary three-axis topology/coordinate-clustered model; and
7. total-hit structured site-label empirical null.

BUSCO >=70% and five-/ten-topology-axis fits are reported as additional
sensitivities and are not part of the seven required checks.

Created: 2026-07-19.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata
import statsmodels


SCRIPT_VERSION = "2026-07-19.1"
ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "ISCIENCE_REVISION_20260711" / "gee_validation"
ANALYSIS_DIR = ROOT / "ISCIENCE_REVISION_20260711" / "analysis_stats"

ROBUST_SCRIPT = ANALYSIS_DIR / "run_robustness_20260711_085930.py"
ENV = OUTDIR / "exact_id_gee_environmental_extraction_20260712_071838.csv"
ENV_MANIFEST = OUTDIR / "exact_id_gee_environmental_extraction_manifest_20260712_071838.json"
DISCOVERY = OUTDIR / "exact_id_gee_raw_pfam_fdr05_20260712_072151.csv"
DISCOVERY_MANIFEST = OUTDIR / "exact_id_gee_correlation_validation_manifest_20260712_072151.json"
RAW = ANALYSIS_DIR / "reconstructed_raw_pfam_counts_20260711_131706.csv.gz"
PEPTIDE = ANALYSIS_DIR / "final_peptide_denominator_manifest_20260711_131706.csv"
PEV = ANALYSIS_DIR / "topology_phylogenetic_eigenvectors_20260711_131706.csv"
MASTER = ROOT / "AlphaEarth" / "CSV" / "Metadata_Table_macroalgae-published.csv"

EXPECTED_SHA256 = {
    ENV: "6d7c1e41eb5651c34464c41c6fadd637db241063d3dfe3c1cc1e614bd6e40418",
    DISCOVERY: "3eb2bbda9d8d9de94f6b13794b664dd3b86dff9788d425fd038805418c1f1b30",
    RAW: "683228342a90d2ecf2897930cd7c147f23973de0444155bfc162d50b09dd22bb",
    PEPTIDE: "5e69ed8cfde0cbe0fa3009d27d6d162e2316d4c4633db85aab45b87ff5d93cb4",
    PEV: "533501f4c48284e23c0525a0e554464284d126fcf69e54a964bd98829e84d056",
    MASTER: "0bff6bf58efae0d22edc1eabefb093c58e8bc5c3280bc5bd4be8466972713921",
    ROBUST_SCRIPT: "ac5c3dc9676dfc69bddca373679445b8ecbde095e88489f17e34b34e0ba226ca",
}

ENV_COLUMNS = [
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
REPRESENTATIONS = ["raw_count", "per_total_pfam_hit", "per_final_peptide_record"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bootstraps", type=int, default=999)
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260719)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def authenticate_inputs() -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    for path, expected in EXPECTED_SHA256.items():
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

    for path in (ENV_MANIFEST, DISCOVERY_MANIFEST):
        if not path.is_file():
            raise FileNotFoundError(path)
    env_manifest = json.loads(ENV_MANIFEST.read_text(encoding="utf-8"))
    if env_manifest["output"]["sha256"] != EXPECTED_SHA256[ENV]:
        raise RuntimeError("Environmental extraction manifest hash does not authenticate CSV")
    discovery_manifest = json.loads(DISCOVERY_MANIFEST.read_text(encoding="utf-8"))
    discovery_key = str(DISCOVERY.relative_to(ROOT))
    if discovery_manifest["outputs"].get(discovery_key) != EXPECTED_SHA256[DISCOVERY]:
        raise RuntimeError("Discovery manifest hash does not authenticate FDR CSV")
    records[ENV_MANIFEST.name] = {
        "path": str(ENV_MANIFEST.relative_to(ROOT)),
        "realpath": str(ENV_MANIFEST.resolve()),
        "sha256": sha256(ENV_MANIFEST),
        "bytes": ENV_MANIFEST.stat().st_size,
    }
    records[DISCOVERY_MANIFEST.name] = {
        "path": str(DISCOVERY_MANIFEST.relative_to(ROOT)),
        "realpath": str(DISCOVERY_MANIFEST.resolve()),
        "sha256": sha256(DISCOVERY_MANIFEST),
        "bytes": DISCOVERY_MANIFEST.stat().st_size,
    }
    return records


def import_robustness_module():
    specification = importlib.util.spec_from_file_location("existing_robustness", ROBUST_SCRIPT)
    if specification is None or specification.loader is None:
        raise ImportError(ROBUST_SCRIPT)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def load_data(robust):
    env = pd.read_csv(ENV, low_memory=False)
    discovery = pd.read_csv(DISCOVERY, low_memory=False)
    raw = pd.read_csv(RAW, compression="gzip", low_memory=False)
    peptide = pd.read_csv(PEPTIDE, low_memory=False)
    pev = pd.read_csv(PEV, low_memory=False)
    master = pd.read_csv(MASTER, encoding="utf-8-sig", low_memory=False)

    if len(env) != 126 or env["genome_id"].nunique() != 126:
        raise RuntimeError("GEE input is not 126 unique exact genome IDs")
    if len(raw) != 126 or raw["Genome"].nunique() != 126:
        raise RuntimeError("Raw Pfam matrix is not 126 unique exact genome IDs")
    master = master[master["Genome"].notna()].copy()
    master = master[master["Genome"].isin(set(env["genome_id"]))].copy()
    if len(master) != 126 or master["Genome"].nunique() != 126:
        raise RuntimeError("Master metadata does not resolve the 126 exact GEE genome IDs once each")
    if len(discovery) != 84 or discovery["sig_fdr05"].fillna(False).astype(bool).sum() != 84:
        raise RuntimeError("Selected GEE family is not the expected 84 global-BH pairs")
    if discovery[["pfam", "env_var"]].duplicated().any():
        raise RuntimeError("Selected GEE family contains duplicate pairs")
    if discovery["pfam"].nunique() != 51:
        raise RuntimeError("Selected GEE family does not contain 51 unique Pfams")
    if int(discovery["sig_bonferroni05"].fillna(False).astype(bool).sum()) != 13:
        raise RuntimeError("Selected GEE family does not contain 13 Bonferroni pairs")
    if not set(discovery["env_var"]).issubset(ENV_COLUMNS):
        raise RuntimeError("Selected GEE family contains an unexpected environmental variable")

    id_sets = [set(env["genome_id"]), set(raw["Genome"]), set(master["Genome"])]
    if id_sets[0] != id_sets[1] or id_sets[0] != id_sets[2]:
        raise RuntimeError("Exact genome-ID sets differ among GEE, Pfam, and metadata inputs")

    raw = raw.set_index("Genome").loc[env["genome_id"]]
    pfam_columns = [column for column in raw.columns if str(column).startswith("PF")]
    raw = raw[pfam_columns].apply(pd.to_numeric, errors="raise")
    if raw.isna().any().any() or (raw < 0).any().any():
        raise RuntimeError("Raw Pfam matrix contains missing or negative values")

    meta_columns = ["Genome", "Phylum", "Nucleotides", "BUSCOs-%present", "DD latitude", "DD longitude"]
    metadata = env.rename(columns={"genome_id": "Genome"}).merge(
        master[meta_columns], on="Genome", how="left", validate="one_to_one"
    ).merge(
        peptide[["Genome", "final_peptide_records"]],
        on="Genome", how="left", validate="one_to_one",
    )
    if metadata["Phylum"].isna().any():
        raise RuntimeError("Master metadata does not cover exact-ID GEE cohort")
    for left, right in (("latitude", "DD latitude"), ("longitude", "DD longitude")):
        difference = (
            pd.to_numeric(metadata[left], errors="raise")
            - pd.to_numeric(metadata[right], errors="raise")
        ).abs()
        if float(difference.max()) > 1e-10:
            raise RuntimeError(f"GEE and master coordinates differ for {left}")
    metadata["Nucleotides"] = pd.to_numeric(metadata["Nucleotides"], errors="raise")
    metadata["BUSCOs-%present"] = pd.to_numeric(metadata["BUSCOs-%present"], errors="raise")
    metadata["total_pfam_hits"] = metadata["Genome"].map(raw.sum(axis=1))
    metadata["site_id"] = (
        metadata["latitude"].round(10).map(lambda value: f"{value:.10f}")
        + "|"
        + metadata["longitude"].round(10).map(lambda value: f"{value:.10f}")
    )
    if metadata["site_id"].nunique() != 90:
        raise RuntimeError("Expected 90 unique coordinate sites")
    metadata = metadata.set_index("Genome", drop=False).loc[raw.index]
    representations = robust.make_representations(raw, metadata)

    pev = pev.set_index("Genome")
    pev_columns = [f"PEV{i}" for i in range(1, 11)]
    if list(pev.columns) != pev_columns or len(pev) != 112:
        raise RuntimeError("Unexpected topology eigenvector matrix")

    missing_pfams = sorted(set(discovery["pfam"]) - set(raw.columns))
    if missing_pfams:
        raise RuntimeError(f"Selected Pfams absent from reconstructed matrix: {missing_pfams}")

    # Recompute every discovery coefficient from the exact joined input before
    # any post hoc analysis. This catches ordering or ID-join errors.
    maximum_difference = 0.0
    for row in discovery.itertuples(index=False):
        exposure = pd.to_numeric(metadata[row.env_var], errors="coerce")
        effect, _ = robust.spearman_effect(
            raw[row.pfam].to_numpy(float), exposure.to_numpy(float)
        )
        maximum_difference = max(maximum_difference, abs(effect - float(row.spearman_r)))
    if maximum_difference > 1e-12:
        raise RuntimeError(f"Discovery coefficient reproduction failed: {maximum_difference}")

    return metadata, discovery, raw, representations, pev, maximum_difference


def result_row(pair, representation: str, method: str, payload: dict[str, object]) -> dict[str, object]:
    return {
        "pfam": pair.pfam,
        "environment": pair.env_var,
        "representation": representation,
        "method": method,
        "discovery_spearman_r": float(pair.spearman_r),
        "discovery_global_bh_q": float(pair.p_fdr),
        "discovery_bonferroni_p": float(pair.p_bonferroni),
        "discovery_bonferroni_significant": bool(pair.sig_bonferroni05),
        **payload,
    }


def analyze_pairs(robust, metadata, discovery, representations, pev, n_boot, rng):
    rows: list[dict[str, object]] = []
    for pair in discovery.itertuples(index=False):
        for representation_name, representation in representations.items():
            outcome = representation.loc[metadata.index, pair.pfam]
            exposure = pd.to_numeric(metadata[pair.env_var], errors="coerce")
            valid = outcome.notna() & exposure.notna()
            m = metadata.loc[valid]

            effect, p_value = robust.spearman_effect(
                outcome.loc[valid].to_numpy(float), exposure.loc[valid].to_numpy(float)
            )
            rows.append(result_row(pair, representation_name, "genome_level_spearman_descriptive", {
                "n": int(valid.sum()),
                "n_sites": int(m["site_id"].nunique()),
                "effect": effect,
                "p_value": p_value,
                "status": "descriptive_only_repeated_site_exposures",
            }))

            site = robust.site_level_frame(outcome.loc[valid], exposure.loc[valid], m)
            site_effect, site_p = robust.spearman_effect(
                site["outcome"].to_numpy(float), site["exposure"].to_numpy(float)
            )
            ci_low, ci_high, kept = robust.bootstrap_spearman_ci(
                site["outcome"].to_numpy(float), site["exposure"].to_numpy(float), n_boot, rng
            )
            rows.append(result_row(pair, representation_name, "site_mean_spearman", {
                "n": len(site),
                "n_sites": len(site),
                "effect": site_effect,
                "p_value": site_p,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "bootstrap_replicates_retained": kept,
                "status": "unique_coordinate_site_inference",
            }))

            quality = robust.partial_rank_cluster(
                outcome.loc[valid], exposure.loc[valid], m,
                cluster=m["site_id"], include_phylum=True, include_peptide=True,
            )
            rows.append(result_row(pair, representation_name, "quality_phylum_sitecluster", {
                "n": quality["n"],
                "n_sites": quality["n_clusters"],
                "effect": quality["effect"],
                "se": quality["se"],
                "p_value": quality["p_value"],
                "ci_low": quality["ci_low"],
                "ci_high": quality["ci_high"],
                "df_resid": quality["df_resid"],
                "condition_number": quality["condition_number"],
                "status": quality["status"],
            }))

            for threshold, suffix in ((50.0, "busco_ge50"), (70.0, "busco_ge70")):
                ids = metadata.index[valid & metadata["BUSCOs-%present"].ge(threshold)]
                filtered_site = robust.site_level_frame(
                    outcome.loc[ids], exposure.loc[ids], metadata.loc[ids]
                )
                filtered_effect, filtered_p = robust.spearman_effect(
                    filtered_site["outcome"].to_numpy(float),
                    filtered_site["exposure"].to_numpy(float),
                )
                low, high, retained = robust.bootstrap_spearman_ci(
                    filtered_site["outcome"].to_numpy(float),
                    filtered_site["exposure"].to_numpy(float), n_boot, rng,
                )
                rows.append(result_row(pair, representation_name, f"site_mean_{suffix}_spearman", {
                    "n": len(filtered_site),
                    "n_sites": len(filtered_site),
                    "effect": filtered_effect,
                    "p_value": filtered_p,
                    "ci_low": low,
                    "ci_high": high,
                    "bootstrap_replicates_retained": retained,
                    "status": "unique_coordinate_site_inference",
                }))

            tree_ids = m.index.intersection(pev.index)
            tree_meta = m.loc[tree_ids]
            for axes, method in (
                (3, "phylo_pev3_quality_sitecluster"),
                (5, "phylo_pev5_quality_sitecluster"),
                (10, "phylo_pev10_quality_sitecluster"),
            ):
                topology = robust.partial_rank_cluster(
                    outcome.loc[tree_ids], exposure.loc[tree_ids], tree_meta,
                    cluster=tree_meta["site_id"], include_phylum=False,
                    include_peptide=True, extra_covariates=pev.loc[tree_ids, pev.columns[:axes]],
                )
                rows.append(result_row(pair, representation_name, method, {
                    "n": topology["n"],
                    "n_sites": topology["n_clusters"],
                    "effect": topology["effect"],
                    "se": topology["se"],
                    "p_value": topology["p_value"],
                    "ci_low": topology["ci_low"],
                    "ci_high": topology["ci_high"],
                    "df_resid": topology["df_resid"],
                    "condition_number": topology["condition_number"],
                    "status": topology["status"],
                }))

    result = pd.DataFrame(rows)
    result = robust.bh_adjust_by_group(
        result,
        ["representation", "method"],
        p_column="p_value",
        output_column="selected_GEE_set_bh_q",
    )
    return result


def structured_null(robust, metadata, discovery, representations, n_perm, rng):
    representation = representations["per_total_pfam_hit"]
    rows: list[dict[str, object]] = []
    for pair in discovery.itertuples(index=False):
        outcome = representation.loc[metadata.index, pair.pfam]
        exposure = pd.to_numeric(metadata[pair.env_var], errors="coerce")
        valid = outcome.notna() & exposure.notna()
        site = robust.site_level_frame(outcome.loc[valid], exposure.loc[valid], metadata.loc[valid])
        x_rank = rankdata(site["outcome"].to_numpy(float), method="average")
        y_rank = rankdata(site["exposure"].to_numpy(float), method="average")
        permuted_y, audit = robust.grouped_permutation_matrix(
            y_rank, site["site_phyla"].to_numpy(), n_perm, rng
        )
        x0 = x_rank - x_rank.mean()
        y0 = y_rank - y_rank.mean()
        observed = float(np.dot(x0, y0) / (np.linalg.norm(x0) * np.linalg.norm(y0)))
        permuted_y -= permuted_y.mean(axis=1, keepdims=True)
        denominator = np.linalg.norm(permuted_y, axis=1) * np.linalg.norm(x0)
        null = np.divide(
            permuted_y @ x0,
            denominator,
            out=np.full(n_perm, np.nan),
            where=denominator > 0,
        )
        finite = null[np.isfinite(null)]
        # Restricted permutations preserve site phylum-composition strata and
        # need not produce a null distribution centered on zero. A two-sided
        # test based on absolute rho would therefore be anti-conservative for
        # some pairs. Use the doubled smaller randomization tail instead.
        lower_tail_p = float((1 + np.sum(finite <= observed)) / (1 + len(finite)))
        upper_tail_p = float((1 + np.sum(finite >= observed)) / (1 + len(finite)))
        empirical = float(min(1.0, 2.0 * min(lower_tail_p, upper_tail_p)))
        rows.append({
            "pfam": pair.pfam,
            "environment": pair.env_var,
            "representation": "per_total_pfam_hit",
            "analysis_level": "unique_coordinate_site_mean",
            "n": len(site),
            "n_sites": len(site),
            "n_genomes_contributing": int(site["n_genomes"].sum()),
            "observed_spearman_rho": observed,
            "empirical_two_sided_p": empirical,
            "empirical_lower_tail_p": lower_tail_p,
            "empirical_upper_tail_p": upper_tail_p,
            "two_sided_method": "two_times_smaller_randomization_tail",
            "null_mean": float(finite.mean()),
            "null_median": float(np.median(finite)),
            "null_sd": float(finite.std(ddof=1)),
            "null_q025": float(np.percentile(finite, 2.5)),
            "null_q975": float(np.percentile(finite, 97.5)),
            "permutations": len(finite),
            **audit,
        })
    result = pd.DataFrame(rows)
    result = robust.bh_adjust_by_group(
        result,
        ["representation"],
        p_column="empirical_two_sided_p",
        output_column="selected_GEE_empirical_bh_q",
    )
    return result


def candidate_summary(results: pd.DataFrame, null: pd.DataFrame, discovery: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    discovery_lookup = discovery.set_index(["pfam", "env_var"])
    for (pfam, environment), group in results.groupby(["pfam", "environment"], sort=False):
        definitions = {
            "raw_site": ("raw_count", "site_mean_spearman"),
            "total_site": ("per_total_pfam_hit", "site_mean_spearman"),
            "peptide_site": ("per_final_peptide_record", "site_mean_spearman"),
            "total_quality": ("per_total_pfam_hit", "quality_phylum_sitecluster"),
            "total_busco50": ("per_total_pfam_hit", "site_mean_busco_ge50_spearman"),
            "total_tree": ("per_total_pfam_hit", "phylo_pev3_quality_sitecluster"),
        }
        checks: dict[str, bool] = {}
        effects: dict[str, float] = {}
        for label, (representation, method) in definitions.items():
            frame = group[
                group["representation"].eq(representation)
                & group["method"].eq(method)
            ]
            if len(frame) != 1:
                raise RuntimeError(f"Missing or duplicate required check: {pfam}, {environment}, {label}")
            checks[label] = bool(frame.iloc[0]["selected_GEE_set_bh_q"] < 0.05)
            effects[label] = float(frame.iloc[0]["effect"])
        empirical = null[
            null["pfam"].eq(pfam) & null["environment"].eq(environment)
        ]
        if len(empirical) != 1:
            raise RuntimeError(f"Missing structured null: {pfam}, {environment}")
        checks["total_structured_null"] = bool(
            empirical.iloc[0]["selected_GEE_empirical_bh_q"] < 0.05
        )
        effects["total_structured_null"] = float(empirical.iloc[0]["observed_spearman_rho"])

        reference_sign = np.sign(effects["raw_site"])
        direction_consistent = bool(
            reference_sign != 0
            and all(np.isfinite(value) and np.sign(value) == reference_sign for value in effects.values())
        )
        source = discovery_lookup.loc[(pfam, environment)]
        rows.append({
            "pfam": pfam,
            "environment": environment,
            "discovery_spearman_r": float(source["spearman_r"]),
            "discovery_global_bh_q": float(source["p_fdr"]),
            "discovery_bonferroni_p": float(source["p_bonferroni"]),
            "discovery_bonferroni_significant": bool(source["sig_bonferroni05"]),
            "reference_raw_site_effect": effects["raw_site"],
            "direction_consistent_required_checks": direction_consistent,
            **{f"selected_set_q_lt_0.05_{label}": value for label, value in checks.items()},
            "checks_passed": int(sum(checks.values())),
            "checks_required": 7,
            "robust_candidate_all_required_checks": bool(direction_consistent and all(checks.values())),
            "interpretation_boundary": "named-variable association; observational covariance does not establish adaptation or mechanism",
        })
    return pd.DataFrame(rows).sort_values(
        ["robust_candidate_all_required_checks", "checks_passed", "discovery_bonferroni_significant", "discovery_global_bh_q"],
        ascending=[False, False, False, True],
        kind="stable",
    )


def write_summary(path: Path, discovery: pd.DataFrame, results: pd.DataFrame, null: pd.DataFrame, candidates: pd.DataFrame) -> None:
    counts = candidates["checks_passed"].value_counts().sort_index(ascending=False)
    lines = [
        "# Primary exact-ID GEE selected-family sensitivity analysis",
        "",
        "## Selection and family",
        "",
        "- Primary discovery family: 87,123 exact-ID raw-count Pfam--GEE tests.",
        "- Selected post hoc family: all 84 pairs meeting global Benjamini--Hochberg q < 0.05, spanning 51 Pfams.",
        "- Thirteen of the 84 selected pairs also met global Bonferroni p < 0.05 in discovery.",
        "- Required post hoc q-values are adjusted across the full 84-pair selected family separately for each check.",
        "",
        "## Seven-check outcome",
        "",
        f"- Pairs passing all seven required checks with consistent direction: {int(candidates['robust_candidate_all_required_checks'].sum())}.",
    ]
    for number, count in counts.items():
        lines.append(f"- Pairs passing {int(number)}/7 checks: {int(count)}.")

    lines.extend(["", "## Bonferroni discovery pairs after sensitivity analysis", ""])
    bonf = candidates[candidates["discovery_bonferroni_significant"]]
    for row in bonf.itertuples(index=False):
        lines.append(
            f"- {row.pfam}--{row.environment}: discovery r={row.discovery_spearman_r:.3f}, "
            f"global q={row.discovery_global_bh_q:.3g}; raw-site r={row.reference_raw_site_effect:.3f}; "
            f"{row.checks_passed}/7 checks; all-required={row.robust_candidate_all_required_checks}."
        )

    lines.extend([
        "",
        "## Interpretation boundary",
        "",
        "These checks evaluate whether selected named-variable associations retain direction and selected-family support across repeated-site, abundance, measured-quality, completeness, topology, and structured-null specifications. They do not establish selection, adaptation, causality, or molecular mechanism.",
        "",
        f"Detailed model rows: `{results.name if getattr(results, 'name', None) else 'see timestamped sensitivity CSV'}`.",
        f"Detailed structured null rows: `{null.name if getattr(null, 'name', None) else 'see timestamped structured-null CSV'}`.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    if args.bootstraps < 99 or args.permutations < 99:
        raise ValueError("At least 99 real-data bootstrap/permutation iterations are required")
    started = datetime.now(timezone.utc)
    start_clock = time.monotonic()
    inputs = authenticate_inputs()
    robust = import_robustness_module()
    metadata, discovery, raw, representations, pev, reproduction_delta = load_data(robust)

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    outputs = {
        "sensitivity": OUTDIR / f"GEE_primary_selected84_sensitivity_{run_id}.csv",
        "structured_null": OUTDIR / f"GEE_primary_selected84_structured_null_{run_id}.csv",
        "candidates": OUTDIR / f"GEE_primary_selected84_candidate_summary_{run_id}.csv",
        "summary": OUTDIR / f"GEE_primary_selected84_summary_{run_id}.md",
        "manifest": OUTDIR / f"GEE_primary_selected84_manifest_{run_id}.json",
    }
    for path in outputs.values():
        if path.exists():
            raise FileExistsError(path)

    rng = np.random.default_rng(args.seed)  # observed-data resampling only
    results = analyze_pairs(
        robust, metadata, discovery, representations, pev, args.bootstraps, rng
    )
    null = structured_null(
        robust, metadata, discovery, representations, args.permutations, rng
    )
    candidates = candidate_summary(results, null, discovery)

    results.to_csv(outputs["sensitivity"], index=False, lineterminator="\n")
    null.to_csv(outputs["structured_null"], index=False, lineterminator="\n")
    candidates.to_csv(outputs["candidates"], index=False, lineterminator="\n")
    write_summary(outputs["summary"], discovery, results, null, candidates)

    # Materialized-output checks.
    reread_results = pd.read_csv(outputs["sensitivity"], low_memory=False)
    reread_null = pd.read_csv(outputs["structured_null"], low_memory=False)
    reread_candidates = pd.read_csv(outputs["candidates"], low_memory=False)
    if len(reread_results) != 84 * 3 * 8:
        raise RuntimeError(f"Unexpected detailed result row count: {len(reread_results)}")
    if len(reread_null) != 84 or len(reread_candidates) != 84:
        raise RuntimeError("Unexpected null or candidate-summary row count")
    if not reread_candidates["checks_required"].eq(7).all():
        raise RuntimeError("Candidate summary does not use seven required checks")

    completed = datetime.now(timezone.utc)
    manifest = {
        "purpose": "selected-family sensitivity analysis for the primary exact-ID GEE findings",
        "selection_rule": "all 84 pairs meeting global BH q < 0.05 in the authenticated 87,123-test exact-ID GEE screen",
        "created_by": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "script_version": SCRIPT_VERSION,
        "started_utc": started.isoformat(),
        "completed_utc": completed.isoformat(),
        "runtime_seconds": time.monotonic() - start_clock,
        "parameters": vars(args),
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "statsmodels": statsmodels.__version__,
        },
        "inputs": inputs,
        "integrity_checks": {
            "real_data_only": True,
            "synthetic_scientific_values": False,
            "randomness_limited_to_observed_data_resampling": True,
            "exact_genome_ids": 126,
            "unique_coordinate_sites": int(metadata["site_id"].nunique()),
            "selected_pairs": len(discovery),
            "selected_unique_pfams": int(discovery["pfam"].nunique()),
            "discovery_bonferroni_pairs": int(discovery["sig_bonferroni05"].sum()),
            "maximum_absolute_discovery_r_reproduction_error": reproduction_delta,
        },
        "results": {
            "pairs_passing_all_seven": int(candidates["robust_candidate_all_required_checks"].sum()),
            "checks_passed_distribution": {
                str(int(key)): int(value)
                for key, value in candidates["checks_passed"].value_counts().sort_index().items()
            },
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
    manifest["outputs"]["manifest"] = {
        "path": str(outputs["manifest"].relative_to(ROOT)),
        "realpath": str(outputs["manifest"].resolve()),
        "sha256_before_self_reference": sha256(outputs["manifest"]),
        "bytes_before_self_reference": outputs["manifest"].stat().st_size,
    }
    outputs["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({
        "run_id": run_id,
        "selected_pairs": len(candidates),
        "pairs_passing_all_seven": int(candidates["robust_candidate_all_required_checks"].sum()),
        "checks_passed_distribution": manifest["results"]["checks_passed_distribution"],
        "runtime_seconds": manifest["runtime_seconds"],
        "outputs": {key: str(path.resolve()) for key, path in outputs.items()},
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

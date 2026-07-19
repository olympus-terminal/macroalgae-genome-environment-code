#!/usr/bin/env python3
"""Refine the selected-family GEE structured null with 99,999 permutations.

This standalone refinement leaves the completed 20260719_205317 run intact.
It reuses that run's authenticated non-null sensitivity rows, generates one
shared set of observed-data site-label permutations for all selected Pfams
within each environmental variable, and rebuilds the 84-pair candidate table.

The conditional null permutes intact observed environmental ranks among unique
coordinate sites only within site phylum-composition strata. The two-sided
randomization p-value is twice the smaller empirical tail because the
conditional null need not be centered at zero.

Randomness is used only to permute observed labels; no scientific values are
simulated or synthesized.
"""

from __future__ import annotations

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
PRIMARY_SCRIPT = OUTDIR / "run_gee_primary_sensitivity_20260719_203058.py"
ROBUST_SCRIPT = ANALYSIS_DIR / "run_robustness_20260711_085930.py"
SOURCE_RESULTS = OUTDIR / "GEE_primary_selected84_sensitivity_20260719_205317.csv"
SOURCE_NULL = OUTDIR / "GEE_primary_selected84_structured_null_20260719_205317.csv"
SOURCE_CANDIDATES = OUTDIR / "GEE_primary_selected84_candidate_summary_20260719_205317.csv"

EXPECTED_SHA256 = {
    PRIMARY_SCRIPT: "7bfef1126ee62ba76e2b5f31fa0ed037a84bc189fe17dacd99ffce15cde8f491",
    ROBUST_SCRIPT: "ac5c3dc9676dfc69bddca373679445b8ecbde095e88489f17e34b34e0ba226ca",
    SOURCE_RESULTS: "4cf5588254756d3cc050fab62f0719a941a431b72e6ee47fe8eee461b5478d15",
    SOURCE_NULL: "5cd6dbe1f4d5c0379aef1abe511c9ac7a1c80a3b8d0e562c1689b6baf9908853",
    SOURCE_CANDIDATES: "a855fe9d5e94297caec9556639907a90a2721a8845b97d76a51a91bdb92a336b",
}

N_PERMUTATIONS = 99_999
BATCH_SIZE = 5_000
SEED = 20260719


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def import_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def authenticate_refinement_inputs() -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    for path, expected in EXPECTED_SHA256.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        observed = sha256(path)
        if observed != expected:
            raise RuntimeError(f"Refinement input changed: {path}: {observed}")
        records[path.name] = {
            "path": str(path.relative_to(ROOT)),
            "realpath": str(path.resolve()),
            "sha256": observed,
            "bytes": path.stat().st_size,
        }
    return records


def refined_structured_null(
    robust,
    metadata: pd.DataFrame,
    discovery: pd.DataFrame,
    representations: dict[str, pd.DataFrame],
    rng: np.random.Generator,
) -> pd.DataFrame:
    representation = representations["per_total_pfam_hit"]
    rows: list[dict[str, object]] = []

    for environment, selected in discovery.groupby("env_var", sort=False):
        exposure = pd.to_numeric(metadata[environment], errors="coerce")
        pfams = selected["pfam"].tolist()
        site_frames = []
        for pfam in pfams:
            outcome = representation.loc[metadata.index, pfam]
            valid = outcome.notna() & exposure.notna()
            site_frames.append(
                robust.site_level_frame(
                    outcome.loc[valid], exposure.loc[valid], metadata.loc[valid]
                )
            )

        reference = site_frames[0]
        for pfam, site in zip(pfams[1:], site_frames[1:]):
            if not site.index.equals(reference.index):
                raise RuntimeError(f"Site index differs within {environment}: {pfam}")
            if not np.array_equal(
                site["exposure"].to_numpy(float),
                reference["exposure"].to_numpy(float),
                equal_nan=True,
            ):
                raise RuntimeError(f"Site exposure differs within {environment}: {pfam}")

        y_rank = rankdata(reference["exposure"].to_numpy(float), method="average")
        y_centered = y_rank - y_rank.mean()
        x_rank = np.column_stack([
            rankdata(site["outcome"].to_numpy(float), method="average")
            for site in site_frames
        ])
        x_centered = x_rank - x_rank.mean(axis=0, keepdims=True)
        x_norm = np.linalg.norm(x_centered, axis=0)
        if np.any(x_norm <= 0) or np.linalg.norm(y_centered) <= 0:
            raise RuntimeError(f"Constant site-rank vector in {environment}")
        observed = (y_centered @ x_centered) / (
            np.linalg.norm(y_centered) * x_norm
        )

        group_labels = reference["site_phyla"].to_numpy(object)
        group_indices = [
            np.flatnonzero(group_labels == label)
            for label in pd.unique(group_labels)
        ]
        group_sizes = np.asarray([len(indices) for indices in group_indices])

        lower_counts = np.zeros(len(pfams), dtype=np.int64)
        upper_counts = np.zeros(len(pfams), dtype=np.int64)
        null_values = np.empty((N_PERMUTATIONS, len(pfams)), dtype=float)
        completed = 0
        while completed < N_PERMUTATIONS:
            batch_n = min(BATCH_SIZE, N_PERMUTATIONS - completed)
            permuted = np.empty((batch_n, len(y_rank)), dtype=float)
            for iteration in range(batch_n):
                values = y_rank.copy()
                for indices in group_indices:
                    if len(indices) > 1:
                        values[indices] = y_rank[rng.permutation(indices)]
                permuted[iteration] = values
            permuted -= permuted.mean(axis=1, keepdims=True)
            denominator = np.linalg.norm(permuted, axis=1)[:, None] * x_norm[None, :]
            batch_null = np.divide(
                permuted @ x_centered,
                denominator,
                out=np.full((batch_n, len(pfams)), np.nan),
                where=denominator > 0,
            )
            if not np.isfinite(batch_null).all():
                raise RuntimeError(f"Nonfinite permutation statistic in {environment}")
            null_values[completed : completed + batch_n] = batch_null
            lower_counts += np.sum(batch_null <= observed[None, :], axis=0)
            upper_counts += np.sum(batch_null >= observed[None, :], axis=0)
            completed += batch_n

        lower_p = (1 + lower_counts) / (1 + N_PERMUTATIONS)
        upper_p = (1 + upper_counts) / (1 + N_PERMUTATIONS)
        two_sided_p = np.minimum(1.0, 2.0 * np.minimum(lower_p, upper_p))

        for column, (pfam, site) in enumerate(zip(pfams, site_frames)):
            null_column = null_values[:, column]
            rows.append({
                "pfam": pfam,
                "environment": environment,
                "representation": "per_total_pfam_hit",
                "analysis_level": "unique_coordinate_site_mean",
                "n": len(site),
                "n_sites": len(site),
                "n_genomes_contributing": int(site["n_genomes"].sum()),
                "observed_spearman_rho": float(observed[column]),
                "empirical_two_sided_p": float(two_sided_p[column]),
                "empirical_lower_tail_p": float(lower_p[column]),
                "empirical_upper_tail_p": float(upper_p[column]),
                "lower_tail_exceedances": int(lower_counts[column]),
                "upper_tail_exceedances": int(upper_counts[column]),
                "two_sided_method": "two_times_smaller_randomization_tail",
                "null_mean": float(null_column.mean()),
                "null_sd": float(null_column.std(ddof=1)),
                "null_q025": float(np.percentile(null_column, 2.5)),
                "null_q975": float(np.percentile(null_column, 97.5)),
                "permutations": N_PERMUTATIONS,
                "seed": SEED,
                "shared_permutation_family": environment,
                "n_groups": len(group_indices),
                "n_multisample_groups": int(np.sum(group_sizes > 1)),
                "n_samples_in_multisample_groups": int(group_sizes[group_sizes > 1].sum()),
                "max_group_size": int(group_sizes.max()),
            })
        del null_values

    result = pd.DataFrame(rows)
    result = robust.bh_adjust_by_group(
        result,
        ["representation"],
        p_column="empirical_two_sided_p",
        output_column="selected_GEE_empirical_bh_q",
    )
    return result


def write_summary(
    path: Path,
    candidates: pd.DataFrame,
    old_candidates: pd.DataFrame,
    refined_null: pd.DataFrame,
) -> None:
    distribution = candidates["checks_passed"].value_counts().sort_index(ascending=False)
    old_lookup = old_candidates.set_index(["pfam", "environment"])
    new_lookup = candidates.set_index(["pfam", "environment"])
    changed = []
    for key in new_lookup.index:
        old_value = bool(old_lookup.loc[key, "robust_candidate_all_required_checks"])
        new_value = bool(new_lookup.loc[key, "robust_candidate_all_required_checks"])
        if old_value != new_value:
            changed.append((key[0], key[1], old_value, new_value))

    bonf = candidates[candidates["discovery_bonferroni_significant"]]
    lines = [
        "# Refined primary GEE structured-null sensitivity",
        "",
        "- Selected family: the same 84 global-BH discovery pairs from 87,123 exact-ID tests.",
        f"- Conditional permutations: {N_PERMUTATIONS:,} observed site-label permutations within site phylum-composition strata.",
        "- One shared permutation stream was used for all selected Pfams within each environmental variable.",
        "- Two-sided p-value: twice the smaller conditional randomization tail; BH correction across all 84 selected pairs.",
        f"- Structured-null pairs with selected-family q < 0.05: {int((refined_null['selected_GEE_empirical_bh_q'] < 0.05).sum())}.",
        f"- Pairs retaining direction and selected-family support in all seven required specifications: {int(candidates['robust_candidate_all_required_checks'].sum())}.",
        f"- Bonferroni discovery pairs retaining all seven: {int(bonf['robust_candidate_all_required_checks'].sum())}/{len(bonf)}.",
        "",
        "## Seven-check distribution",
        "",
    ]
    for number, count in distribution.items():
        lines.append(f"- {int(number)}/7: {int(count)} pairs.")
    lines.extend(["", "## Classifications changed from the 9,999-permutation run", ""])
    if changed:
        for pfam, environment, old_value, new_value in changed:
            lines.append(f"- {pfam}--{environment}: all-seven {old_value} -> {new_value}.")
    else:
        lines.append("- None.")
    lines.extend([
        "",
        "## Interpretation boundary",
        "",
        "The permutation conditions on exact coordinate sites and site phylum-composition strata. It does not preserve geographic distance or spatial autocorrelation among distinct sites, and the seven specifications are correlated re-analyses of the same selected data rather than independent validations.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    started = datetime.now(timezone.utc)
    start_clock = time.monotonic()
    refinement_inputs = authenticate_refinement_inputs()
    primary = import_module("primary_gee_sensitivity", PRIMARY_SCRIPT)
    robust = import_module("existing_robustness", ROBUST_SCRIPT)
    primary_inputs = primary.authenticate_inputs()
    metadata, discovery, raw, representations, pev, reproduction_delta = primary.load_data(robust)
    del raw, pev

    results = pd.read_csv(SOURCE_RESULTS, low_memory=False)
    old_null = pd.read_csv(SOURCE_NULL, low_memory=False)
    old_candidates = pd.read_csv(SOURCE_CANDIDATES, low_memory=False)
    if len(results) != 84 * 3 * 8 or len(old_null) != 84 or len(old_candidates) != 84:
        raise RuntimeError("Unexpected source-output row count")
    if results[["pfam", "environment"]].drop_duplicates().shape[0] != 84:
        raise RuntimeError("Source sensitivity table does not contain 84 unique pairs")

    rng = np.random.default_rng(SEED)
    refined_null = refined_structured_null(
        robust, metadata, discovery, representations, rng
    )
    candidates = primary.candidate_summary(results, refined_null, discovery)

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    outputs = {
        "structured_null": OUTDIR / f"GEE_primary_selected84_structured_null_refined99999_{run_id}.csv",
        "candidates": OUTDIR / f"GEE_primary_selected84_candidate_summary_refined99999_{run_id}.csv",
        "summary": OUTDIR / f"GEE_primary_selected84_summary_refined99999_{run_id}.md",
        "manifest": OUTDIR / f"GEE_primary_selected84_manifest_refined99999_{run_id}.json",
    }
    for path in outputs.values():
        if path.exists():
            raise FileExistsError(path)

    refined_null.to_csv(outputs["structured_null"], index=False, lineterminator="\n")
    candidates.to_csv(outputs["candidates"], index=False, lineterminator="\n")
    write_summary(outputs["summary"], candidates, old_candidates, refined_null)

    reread_null = pd.read_csv(outputs["structured_null"], low_memory=False)
    reread_candidates = pd.read_csv(outputs["candidates"], low_memory=False)
    if len(reread_null) != 84 or len(reread_candidates) != 84:
        raise RuntimeError("Refined output row count changed on materialization")
    if not reread_null["permutations"].eq(N_PERMUTATIONS).all():
        raise RuntimeError("Refined null does not record 99,999 permutations")
    expected_two_sided = np.minimum(
        1.0,
        2.0 * np.minimum(
            reread_null["empirical_lower_tail_p"],
            reread_null["empirical_upper_tail_p"],
        ),
    )
    if float(np.max(np.abs(expected_two_sided - reread_null["empirical_two_sided_p"]))) > 1e-15:
        raise RuntimeError("Materialized two-sided p-values do not match conditional tails")
    if not reread_candidates["checks_required"].eq(7).all():
        raise RuntimeError("Refined candidate table does not use seven checks")

    completed = datetime.now(timezone.utc)
    old_all_seven = int(old_candidates["robust_candidate_all_required_checks"].sum())
    new_all_seven = int(candidates["robust_candidate_all_required_checks"].sum())
    manifest = {
        "purpose": "99,999-permutation refinement of the primary exact-ID GEE selected-family structured null",
        "created_by": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "script_version": SCRIPT_VERSION,
        "started_utc": started.isoformat(),
        "completed_utc": completed.isoformat(),
        "runtime_seconds": time.monotonic() - start_clock,
        "parameters": {
            "permutations": N_PERMUTATIONS,
            "batch_size": BATCH_SIZE,
            "seed": SEED,
            "two_sided_method": "two_times_smaller_randomization_tail",
            "permutation_unit": "intact unique-coordinate-site environmental rank",
            "permutation_strata": "site phylum-composition",
            "shared_draws_within": "environmental variable",
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "statsmodels": statsmodels.__version__,
        },
        "inputs": {**refinement_inputs, **primary_inputs},
        "integrity_checks": {
            "real_data_only": True,
            "synthetic_scientific_values": False,
            "randomness_limited_to_permutation_of_observed_site_labels": True,
            "exact_genome_ids": 126,
            "unique_coordinate_sites": int(metadata["site_id"].nunique()),
            "selected_pairs": len(discovery),
            "selected_unique_pfams": int(discovery["pfam"].nunique()),
            "discovery_bonferroni_pairs": int(discovery["sig_bonferroni05"].sum()),
            "maximum_absolute_discovery_r_reproduction_error": reproduction_delta,
            "materialized_outputs_reopened": True,
        },
        "results": {
            "structured_null_q_lt_0_05": int((refined_null["selected_GEE_empirical_bh_q"] < 0.05).sum()),
            "old_9999_pairs_passing_all_seven": old_all_seven,
            "refined_99999_pairs_passing_all_seven": new_all_seven,
            "bonferroni_discovery_pairs_passing_all_seven": int(
                candidates.loc[
                    candidates["discovery_bonferroni_significant"],
                    "robust_candidate_all_required_checks",
                ].sum()
            ),
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
        "permutations": N_PERMUTATIONS,
        "structured_null_q_lt_0_05": manifest["results"]["structured_null_q_lt_0_05"],
        "pairs_passing_all_seven": new_all_seven,
        "bonferroni_pairs_passing_all_seven": manifest["results"]["bonferroni_discovery_pairs_passing_all_seven"],
        "runtime_seconds": manifest["runtime_seconds"],
        "outputs": {label: str(path.resolve()) for label, path in outputs.items()},
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

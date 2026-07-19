#!/usr/bin/env python3
"""Correct selected AEF structured-null p-values using conditional tails.

This program reuses the real-data site aggregation and within-stratum
permutation machinery from ``run_robustness_20260711_085930.py``.  It does not
simulate observations or reconstruct scientific values.  The only random
operation is permutation of observed, intact coordinate-site units within the
observed site phylum-composition strata.

The archived implementation compared ``abs(T_perm)`` with ``abs(T_obs)``.
That comparison is not a valid two-sided test when the conditional permutation
distribution is not centered on zero.  Here the two-sided p-value is

    p_lo  = (1 + number(T_perm <= T_obs)) / (B + 1)
    p_hi  = (1 + number(T_perm >= T_obs)) / (B + 1)
    p_two = min(1, 2 * min(p_lo, p_hi)).

Benjamini-Hochberg correction is applied across the complete selected family
of 68 pairs separately for each abundance representation.  Existing non-null
sensitivity results are read without modification to rebuild the seven-check
candidate table.  Every output is timestamped and existing files are never
overwritten.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import rankdata
import statsmodels
from statsmodels.stats.multitest import multipletests

import run_robustness_20260711_085930 as archived


SCRIPT_VERSION = "2026-07-19.1"
DEFAULT_PERMUTATIONS = 99_999
DEFAULT_SEED = 20_260_711
REPRESENTATIONS = (
    "raw_count",
    "per_total_pfam_hit",
    "per_final_peptide_record",
)
EXPLICIT_PAIRS = (
    ("PF13411", "A18"),
    ("PF01638", "A52"),
    ("PF01638", "A53"),
    ("PF10988", "A36"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate all real-data joins and archived inputs without writing outputs.",
    )
    return parser.parse_args()


def now_stamp() -> str:
    return datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(chunk_size), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, object]:
    path = path.resolve()
    stat = path.stat()
    return {
        "path": str(path),
        "bytes": stat.st_size,
        "mtime_utc": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
        "sha256": sha256_file(path),
    }


def require_file(path: Path) -> Path:
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def selected_pairs() -> list[tuple[str, str]]:
    axes = [f"A{index:02d}" for index in range(64)]
    pairs = list(EXPLICIT_PAIRS) + [("PF00092", axis) for axis in axes]
    if len(pairs) != 68 or len(set(pairs)) != 68:
        raise AssertionError("Selected AEF family must contain 68 unique pairs")
    return pairs


def load_real_inputs(root: Path) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    dict[str, pd.DataFrame],
    pd.DataFrame,
    pd.DataFrame,
    dict[str, Path],
]:
    analysis_dir = root / "ISCIENCE_REVISION_20260711" / "analysis_stats"
    paths = {
        "archived_script": require_file(analysis_dir / "run_robustness_20260711_085930.py"),
        "raw_counts": require_file(
            analysis_dir / "reconstructed_raw_pfam_counts_20260711_131706.csv.gz"
        ),
        "peptide_denominators": require_file(
            analysis_dir / "final_peptide_denominator_manifest_20260711_131706.csv"
        ),
        "archived_non_null": require_file(
            analysis_dir / "archived_AEF_priority_robustness_20260711_131706.csv"
        ),
        "archived_structured_null": require_file(
            analysis_dir / "archived_AEF_within_phylum_null_20260711_131706.csv"
        ),
        "archived_candidates": require_file(
            analysis_dir / "robust_candidate_table_20260711_131706.csv"
        ),
        "archived_manifest": require_file(
            analysis_dir / "run_manifest_20260711_131706.json"
        ),
    }

    source_inputs = archived.load_inputs(root)
    cohort, aef, _legacy = archived.load_canonical_cohort(source_inputs)
    paths["aef_embeddings"] = require_file(source_inputs["aef"])
    paths["master_metadata"] = require_file(source_inputs["master"])

    raw_pfams = pd.read_csv(paths["raw_counts"], compression="gzip")
    if "Genome" not in raw_pfams.columns or raw_pfams["Genome"].duplicated().any():
        raise ValueError("Archived raw Pfam matrix lacks unique Genome identifiers")
    raw_pfams = raw_pfams.set_index("Genome")
    if set(raw_pfams.index) != set(cohort["Genome"]):
        raise ValueError("Archived raw Pfam matrix does not match the canonical cohort")
    raw_pfams = raw_pfams.loc[cohort["Genome"]]

    required_pfams = sorted({pfam for pfam, _axis in selected_pairs()})
    missing_pfams = sorted(set(required_pfams) - set(raw_pfams.columns))
    if missing_pfams:
        raise ValueError(f"Selected Pfams missing from raw matrix: {missing_pfams}")

    peptide = pd.read_csv(paths["peptide_denominators"])
    if peptide["Genome"].duplicated().any():
        raise ValueError("Peptide denominator manifest contains duplicate Genome IDs")
    cohort = cohort.merge(
        peptide[["Genome", "final_peptide_records"]],
        on="Genome",
        how="left",
        validate="one_to_one",
    )
    cohort["total_pfam_hits"] = cohort["Genome"].map(raw_pfams.sum(axis=1))
    cohort["site_id"] = (
        cohort["DD latitude"].round(10).map(lambda value: f"{value:.10f}")
        + "|"
        + cohort["DD longitude"].round(10).map(lambda value: f"{value:.10f}")
    )
    if cohort["site_id"].nunique() != 90:
        raise ValueError("Canonical cohort no longer resolves to 90 exact coordinate sites")

    representations = archived.make_representations(raw_pfams, cohort)
    representations = {
        name: frame.loc[cohort["Genome"], required_pfams]
        for name, frame in representations.items()
    }
    if tuple(representations) != REPRESENTATIONS:
        raise ValueError(f"Unexpected representation order: {tuple(representations)}")

    non_null = pd.read_csv(paths["archived_non_null"])
    legacy_null = pd.read_csv(paths["archived_structured_null"])
    legacy_candidates = pd.read_csv(paths["archived_candidates"])
    expected_pair_set = set(selected_pairs())
    for frame, label in ((non_null, "non-null"), (legacy_null, "structured-null")):
        observed = set(zip(frame["pfam"], frame["latent_axis"]))
        if observed != expected_pair_set:
            raise ValueError(f"Archived {label} selected pair family changed")
    if len(legacy_null) != 68 * len(REPRESENTATIONS):
        raise ValueError("Archived structured-null table is not 68 pairs x 3 representations")

    metadata = cohort.set_index("Genome", drop=False)
    embedding = aef.set_index("Genome").loc[cohort["Genome"]]
    return metadata, embedding, representations, non_null, legacy_null, legacy_candidates, paths


def two_sided_conditional_tail(
    null: np.ndarray,
    observed: np.ndarray,
) -> dict[str, np.ndarray]:
    """Return doubled-smaller-tail p-values and Monte Carlo diagnostics."""
    null = np.asarray(null, dtype=float)
    observed = np.asarray(observed, dtype=float)
    if null.ndim != 2 or observed.shape != (null.shape[1],):
        raise ValueError("Null/observed shape mismatch")
    finite = np.isfinite(null)
    n = finite.sum(axis=0)
    if (n == 0).any() or not np.isfinite(observed).all():
        raise ValueError("Nonfinite observed statistic or empty permutation distribution")
    count_lo = ((null <= observed[None, :]) & finite).sum(axis=0)
    count_hi = ((null >= observed[None, :]) & finite).sum(axis=0)
    p_lo = (1.0 + count_lo) / (1.0 + n)
    p_hi = (1.0 + count_hi) / (1.0 + n)
    smaller = np.minimum(p_lo, p_hi)
    p_two = np.minimum(1.0, 2.0 * smaller)
    # Approximate Monte Carlo standard error for the doubled selected tail.
    mc_se = 2.0 * np.sqrt(smaller * (1.0 - smaller) / (1.0 + n))
    return {
        "finite": n,
        "count_le_observed": count_lo,
        "count_ge_observed": count_hi,
        "p_lower_conditional": p_lo,
        "p_upper_conditional": p_hi,
        "empirical_two_sided_p_conditional_tail": p_two,
        "mc_se_two_sided_approx": mc_se,
    }


def prepare_representation_sites(
    pair_family: Iterable[tuple[str, str]],
    rep_name: str,
    representation: pd.DataFrame,
    metadata: pd.DataFrame,
    embedding: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[tuple[str, str], dict[str, object]]]:
    records: dict[tuple[str, str], dict[str, object]] = {}
    reference_site: pd.DataFrame | None = None
    for pfam, axis in pair_family:
        outcome = representation.loc[metadata.index, pfam]
        exposure = pd.to_numeric(embedding.loc[metadata.index, axis], errors="coerce")
        valid = outcome.notna() & exposure.notna()
        site = archived.site_level_frame(
            outcome.loc[valid], exposure.loc[valid], metadata.loc[valid]
        )
        if reference_site is None:
            reference_site = site[["site_phyla", "n_genomes"]].copy()
        else:
            if not site.index.equals(reference_site.index):
                raise ValueError(
                    f"Site membership differs within representation {rep_name}: {pfam}-{axis}"
                )
            if not site["site_phyla"].equals(reference_site["site_phyla"]):
                raise ValueError("Site phylum-composition strata changed across pairs")
        x_rank = rankdata(site["outcome"].to_numpy(float), method="average")
        y_rank = rankdata(site["exposure"].to_numpy(float), method="average")
        x0 = x_rank - x_rank.mean()
        y0 = y_rank - y_rank.mean()
        denominator = np.linalg.norm(x0) * np.linalg.norm(y0)
        if denominator <= 0:
            raise ValueError(f"Constant site rank vector for {pfam}-{axis}/{rep_name}")
        observed = float(np.dot(x0, y0) / denominator)
        records[(pfam, axis)] = {
            "site": site,
            "x0": x0,
            "y0": y0,
            "x_norm": float(np.linalg.norm(x0)),
            "y_norm": float(np.linalg.norm(y0)),
            "observed": observed,
        }
    if reference_site is None:
        raise ValueError("No selected pairs were prepared")
    return reference_site, records


def corrected_permutation_for_representation(
    pair_family: list[tuple[str, str]],
    rep_name: str,
    representation: pd.DataFrame,
    metadata: pd.DataFrame,
    embedding: pd.DataFrame,
    n_perm: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, dict[str, object]]:
    reference_site, records = prepare_representation_sites(
        pair_family, rep_name, representation, metadata, embedding
    )
    n_sites = len(reference_site)
    # Reuse the archived routine literally to obtain within-stratum random
    # permutation indices.  Applying these indices to each observed outcome is
    # distributionally equivalent to permuting the exposure in the archived
    # code because the same within-stratum permutation group is used.
    permuted_index_float, group_audit = archived.grouped_permutation_matrix(
        np.arange(n_sites, dtype=float),
        reference_site["site_phyla"].to_numpy(),
        n_perm,
        rng,
    )
    permuted_index = permuted_index_float.astype(np.int16, copy=False)
    if not np.array_equal(np.sort(permuted_index, axis=1)[0], np.arange(n_sites)):
        raise AssertionError("Permutation index does not preserve intact site units")
    del permuted_index_float

    rows: list[dict[str, object]] = []
    by_pfam: dict[str, list[str]] = {}
    for pfam, axis in pair_family:
        by_pfam.setdefault(pfam, []).append(axis)

    for pfam, axes in by_pfam.items():
        first = records[(pfam, axes[0])]
        permuted_x = np.asarray(first["x0"])[permuted_index]
        # Permutation preserves the centered-vector norm, but calculate each
        # row explicitly to retain the archived numerical definition.
        permuted_x -= permuted_x.mean(axis=1, keepdims=True)
        x_norm = np.linalg.norm(permuted_x, axis=1)
        y_matrix = np.column_stack([records[(pfam, axis)]["y0"] for axis in axes])
        y_norm = np.asarray([records[(pfam, axis)]["y_norm"] for axis in axes])
        denominator = x_norm[:, None] * y_norm[None, :]
        null = np.divide(
            permuted_x @ y_matrix,
            denominator,
            out=np.full((n_perm, len(axes)), np.nan, dtype=float),
            where=denominator > 0,
        )
        observed = np.asarray([records[(pfam, axis)]["observed"] for axis in axes])
        tails = two_sided_conditional_tail(null, observed)
        null_q025, null_q975 = np.percentile(null, [2.5, 97.5], axis=0)
        for column, axis in enumerate(axes):
            site = records[(pfam, axis)]["site"]
            rows.append(
                {
                    "pfam": pfam,
                    "latent_axis": axis,
                    "representation": rep_name,
                    "analysis_level": "unique_coordinate_site_mean",
                    "permutation_structure": "intact_sites_within_site_phylum_composition",
                    "two_sided_method": "doubled_smaller_conditional_tail_including_ties",
                    "n": len(site),
                    "n_sites": len(site),
                    "n_genomes_contributing": int(site["n_genomes"].sum()),
                    "observed_spearman_rho": observed[column],
                    "count_perm_le_observed": int(tails["count_le_observed"][column]),
                    "count_perm_ge_observed": int(tails["count_ge_observed"][column]),
                    "p_lower_conditional": tails["p_lower_conditional"][column],
                    "p_upper_conditional": tails["p_upper_conditional"][column],
                    "empirical_two_sided_p_conditional_tail": tails[
                        "empirical_two_sided_p_conditional_tail"
                    ][column],
                    "mc_se_two_sided_approx": tails["mc_se_two_sided_approx"][column],
                    "null_mean": float(np.mean(null[:, column])),
                    "null_sd": float(np.std(null[:, column], ddof=1)),
                    "null_q025": float(null_q025[column]),
                    "null_q975": float(null_q975[column]),
                    "permutations": int(tails["finite"][column]),
                    "minimum_attainable_two_sided_p": 2.0 / (n_perm + 1.0),
                    **group_audit,
                }
            )
        del permuted_x, null

    result = pd.DataFrame(rows)
    if len(result) != 68:
        raise AssertionError(f"Expected 68 rows for {rep_name}; observed {len(result)}")
    return result, {"representation": rep_name, "n_sites": n_sites, **group_audit}


def add_bh_and_legacy_comparison(
    corrected: pd.DataFrame,
    legacy_null: pd.DataFrame,
) -> pd.DataFrame:
    corrected = corrected.copy()
    corrected["selected_AEF_empirical_bh_q_conditional_tail"] = np.nan
    for representation, indices in corrected.groupby("representation").groups.items():
        p = corrected.loc[indices, "empirical_two_sided_p_conditional_tail"]
        corrected.loc[indices, "selected_AEF_empirical_bh_q_conditional_tail"] = (
            multipletests(p.to_numpy(float), method="fdr_bh")[1]
        )
    legacy_columns = [
        "pfam",
        "latent_axis",
        "representation",
        "observed_spearman_rho",
        "empirical_two_sided_p",
        "selected_AEF_empirical_bh_q",
    ]
    old = legacy_null[legacy_columns].rename(
        columns={
            "observed_spearman_rho": "archived_observed_spearman_rho",
            "empirical_two_sided_p": "archived_absolute_tail_p",
            "selected_AEF_empirical_bh_q": "archived_absolute_tail_bh_q",
        }
    )
    merged = corrected.merge(
        old,
        on=["pfam", "latent_axis", "representation"],
        how="left",
        validate="one_to_one",
    )
    difference = (
        merged["observed_spearman_rho"] - merged["archived_observed_spearman_rho"]
    ).abs()
    if difference.max() > 1e-12:
        raise ValueError(
            "Recomputed observed site effects differ from archived values; "
            f"maximum absolute difference={difference.max()}"
        )
    merged["conditional_tail_q_lt_0.05"] = (
        merged["selected_AEF_empirical_bh_q_conditional_tail"] < 0.05
    )
    merged["archived_absolute_tail_q_lt_0.05"] = (
        merged["archived_absolute_tail_bh_q"] < 0.05
    )
    merged["significance_status_changed"] = (
        merged["conditional_tail_q_lt_0.05"]
        != merged["archived_absolute_tail_q_lt_0.05"]
    )
    return merged.sort_values(["representation", "pfam", "latent_axis"])


def rebuild_candidates(
    non_null: pd.DataFrame,
    corrected_null: pd.DataFrame,
    legacy_candidates: pd.DataFrame,
) -> pd.DataFrame:
    perm_for_builder = corrected_null[
        ["pfam", "latent_axis", "representation"]
    ].copy()
    perm_for_builder["selected_AEF_empirical_bh_q"] = corrected_null[
        "selected_AEF_empirical_bh_q_conditional_tail"
    ].to_numpy(float)
    rebuilt = archived.build_robust_candidate_table(non_null, perm_for_builder)
    annotation_columns = [
        "pfam",
        "raw_hmm_query_name",
        "verified_interpro_short_name",
        "verified_interpro_name",
        "annotation_source_database",
        "annotation_api_url",
        "annotation_http_status",
        "annotation_response_sha256",
    ]
    annotations = legacy_candidates[annotation_columns].drop_duplicates("pfam")
    rebuilt = rebuilt.merge(annotations, on="pfam", how="left", validate="many_to_one")
    old = legacy_candidates[
        [
            "pfam",
            "latent_axis",
            "checks_passed",
            "robust_candidate_all_required_checks",
            "empirical_selected_set_q_lt_0.05_total_hit_within_phylum",
        ]
    ].rename(
        columns={
            "checks_passed": "archived_checks_passed",
            "robust_candidate_all_required_checks": "archived_robust_candidate_all_required_checks",
            "empirical_selected_set_q_lt_0.05_total_hit_within_phylum": (
                "archived_empirical_q_lt_0.05_total_hit"
            ),
        }
    )
    rebuilt = rebuilt.merge(old, on=["pfam", "latent_axis"], validate="one_to_one")
    rebuilt["all_seven_status_changed"] = (
        rebuilt["robust_candidate_all_required_checks"]
        != rebuilt["archived_robust_candidate_all_required_checks"]
    )
    rebuilt["structured_null_p_method"] = "doubled_smaller_conditional_tail"
    return rebuilt.sort_values(
        ["robust_candidate_all_required_checks", "checks_passed", "pfam", "latent_axis"],
        ascending=[False, False, True, True],
    )


def write_summary(
    path: Path,
    corrected: pd.DataFrame,
    candidates: pd.DataFrame,
    n_perm: int,
) -> None:
    lines = [
        "# Corrected selected-AEF structured-null audit",
        "",
        f"Permutation replicates: {n_perm:,} per abundance representation.",
        (
            "Two-sided p-values use twice the smaller conditional tail: "
            "p_lo=(1+#Tperm<=Tobs)/(B+1), p_hi=(1+#Tperm>=Tobs)/(B+1), "
            "p_two=min(1,2*min(p_lo,p_hi))."
        ),
        "BH correction spans all 68 selected pairs separately within each representation.",
        "",
        "## Aggregate results",
        "",
    ]
    for representation, group in corrected.groupby("representation", sort=False):
        lines.append(
            f"- {representation}: {(group['selected_AEF_empirical_bh_q_conditional_tail'] < 0.05).sum()} "
            f"of {len(group)} pairs have corrected BH q < 0.05."
        )
    lines.extend(
        [
            f"- All-seven robust candidates: {int(candidates['robust_candidate_all_required_checks'].sum())} of {len(candidates)}.",
            f"- Candidate all-seven status changes: {int(candidates['all_seven_status_changed'].sum())}.",
            "",
            "## Priority pairs",
            "",
        ]
    )
    priority = candidates[candidates["pfam"].ne("PF00092")]
    total_null = corrected[corrected["representation"].eq("per_total_pfam_hit")]
    for row in priority.itertuples(index=False):
        null = total_null[
            total_null["pfam"].eq(row.pfam)
            & total_null["latent_axis"].eq(row.latent_axis)
        ].iloc[0]
        lines.append(
            f"- {row.pfam}-{row.latent_axis}: rho={null.observed_spearman_rho:.6g}; "
            f"corrected p={null.empirical_two_sided_p_conditional_tail:.6g}; "
            f"corrected q={null.selected_AEF_empirical_bh_q_conditional_tail:.6g}; "
            f"checks={row.checks_passed}/{row.checks_required}; "
            f"all-seven={bool(row.robust_candidate_all_required_checks)}."
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    if args.permutations < 99:
        raise ValueError("At least 99 real-data permutations are required")
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    output_dir = script_path.parent
    (
        metadata,
        embedding,
        representations,
        non_null,
        legacy_null,
        legacy_candidates,
        input_paths,
    ) = load_real_inputs(root)
    family = selected_pairs()
    if args.validate_only:
        print(
            json.dumps(
                {
                    "status": "validated",
                    "genomes": len(metadata),
                    "sites": int(metadata["site_id"].nunique()),
                    "selected_pairs": len(family),
                    "representations": list(representations),
                    "traceable_peptide_denominators": int(
                        metadata["final_peptide_records"].notna().sum()
                    ),
                },
                indent=2,
            )
        )
        return 0

    run_id = now_stamp()
    outputs = {
        "corrected_null": output_dir
        / f"corrected_AEF_structured_null_conditional_tail_{run_id}.csv",
        "corrected_candidates": output_dir
        / f"corrected_AEF_robust_candidate_table_{run_id}.csv",
        "summary": output_dir
        / f"corrected_AEF_structured_null_summary_{run_id}.md",
        "manifest": output_dir
        / f"corrected_AEF_structured_null_manifest_{run_id}.json",
    }
    existing = [str(path) for path in outputs.values() if path.exists()]
    if existing:
        raise FileExistsError(f"Refusing to overwrite outputs: {existing}")

    seed_sequence = np.random.SeedSequence(args.seed)
    child_seeds = seed_sequence.spawn(len(REPRESENTATIONS))
    frames: list[pd.DataFrame] = []
    group_audits: list[dict[str, object]] = []
    for rep_name, child_seed in zip(REPRESENTATIONS, child_seeds):
        print(
            f"Permuting {len(family)} selected pairs for {rep_name}: "
            f"B={args.permutations:,}",
            flush=True,
        )
        frame, audit = corrected_permutation_for_representation(
            family,
            rep_name,
            representations[rep_name],
            metadata,
            embedding,
            args.permutations,
            np.random.default_rng(child_seed),
        )
        frames.append(frame)
        group_audits.append(audit)

    corrected = add_bh_and_legacy_comparison(pd.concat(frames, ignore_index=True), legacy_null)
    candidates = rebuild_candidates(non_null, corrected, legacy_candidates)
    corrected.to_csv(outputs["corrected_null"], index=False)
    candidates.to_csv(outputs["corrected_candidates"], index=False)
    write_summary(outputs["summary"], corrected, candidates, args.permutations)

    manifest: dict[str, object] = {
        "status": "complete",
        "run_id": run_id,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "script_version": SCRIPT_VERSION,
        "script": file_record(script_path),
        "command": sys.argv,
        "parameters": {
            "permutations": args.permutations,
            "seed": args.seed,
            "selected_pairs": len(family),
            "representations": list(REPRESENTATIONS),
            "bh_family": "68 selected pairs separately within each representation",
            "two_sided_method": "doubled smaller conditional tail including ties",
            "formula": {
                "p_lower": "(1 + count(T_perm <= T_obs)) / (B + 1)",
                "p_upper": "(1 + count(T_perm >= T_obs)) / (B + 1)",
                "p_two": "min(1, 2 * min(p_lower, p_upper))",
            },
            "minimum_attainable_two_sided_p": 2.0 / (args.permutations + 1.0),
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "statsmodels": statsmodels.__version__,
        },
        "inputs": {name: file_record(path) for name, path in input_paths.items()},
        "data_counts": {
            "canonical_genomes": len(metadata),
            "unique_coordinate_sites_raw_total": int(metadata["site_id"].nunique()),
            "traceable_peptide_denominators": int(
                metadata["final_peptide_records"].notna().sum()
            ),
            "selected_pairs": len(family),
            "structured_null_rows": len(corrected),
            "candidate_rows": len(candidates),
        },
        "group_audits": group_audits,
        "results": {
            "corrected_bh_q_lt_0.05_by_representation": {
                representation: int(
                    (group["selected_AEF_empirical_bh_q_conditional_tail"] < 0.05).sum()
                )
                for representation, group in corrected.groupby("representation")
            },
            "legacy_to_corrected_significance_changes": int(
                corrected["significance_status_changed"].sum()
            ),
            "all_seven_robust_candidates": int(
                candidates["robust_candidate_all_required_checks"].sum()
            ),
            "all_seven_status_changes": int(candidates["all_seven_status_changed"].sum()),
            "checks_passed_distribution": {
                str(int(key)): int(value)
                for key, value in candidates["checks_passed"].value_counts().sort_index().items()
            },
        },
        "outputs": {},
    }
    # Record hashes only after all non-manifest outputs are complete.
    manifest["outputs"] = {
        name: file_record(path)
        for name, path in outputs.items()
        if name != "manifest"
    }
    outputs["manifest"].write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print("Corrected structured-null analysis complete.", flush=True)
    for name, path in outputs.items():
        print(f"  {name}: {path.resolve()}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

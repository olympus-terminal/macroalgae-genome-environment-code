#!/usr/bin/env python3
"""Inventory and integrity-check the public code release without test data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_VERSION = "2026-07-19.2"
TEXT_SUFFIXES = {
    ".py",
    ".md",
    ".txt",
    ".cff",
    ".json",
    ".yaml",
    ".yml",
    ".toml",
    ".tsv",
    ".csv",
    ".gitignore",
}
REQUIRED_PATHS = {
    "README.md",
    "LICENSE",
    "CITATION.cff",
    "DATA_INPUTS.md",
    "MANUAL_STEPS.md",
    "SCIENTIFIC_INTEGRITY.md",
    "SOFTWARE_CITATIONS.md",
    "SOFTWARE_ENVIRONMENTS.md",
    "requirements.txt",
    "ISCIENCE_REVISION_20260711/aef/recompute_full_aef_pfam_analysis_20260718_222823.py",
    "ISCIENCE_REVISION_20260711/aef/extract_exact_id_aef_embeddings_20260718_224936.py",
    "ISCIENCE_REVISION_20260711/aef/archived_reported_extractor/README.md",
    "ISCIENCE_REVISION_20260711/aef/archived_reported_extractor/extract_aef_embeddings_20251019.py",
    "ISCIENCE_REVISION_20260711/figure3_126/rebuild_figure3_recorded_metadata_20260718_223224.py",
    "ISCIENCE_REVISION_20260711/figure4_126/rebuild_figure4_126_20260715_232158.py",
    "ISCIENCE_REVISION_20260711/figure5_gee_primary/build_figure5_primary_gee_sensitivity_20260719_210207.py",
    "ISCIENCE_REVISION_20260711/gee_validation/run_exact_id_gee_correlation_validation_20260712_072151.py",
    "ISCIENCE_REVISION_20260711/gee_validation/run_aef_gee_site_alignment_20260719_204821.py",
    "ISCIENCE_REVISION_20260711/gee_validation/run_gee_primary_sensitivity_20260719_203058.py",
    "ISCIENCE_REVISION_20260711/gee_validation/run_gee_structured_null_refinement_20260719_205905.py",
    "ISCIENCE_REVISION_20260711/analysis_stats/run_robustness_20260711_085930.py",
    "ISCIENCE_REVISION_20260711/analysis_stats/correct_aef_structured_null_conditional_tail_20260719_205621.py",
    "ISCIENCE_REVISION_20260711/annotations/fetch_gee_robust_interpro_20260719_212508.py",
    "ISCIENCE_REVISION_20260711/pf00092_midas/test_pf00092_midas_20260715_232919.py",
    "ISCIENCE_REVISION_20260711/supplement/build_corrected_table_s3_20260718_224354.py",
    "ISCIENCE_REVISION_20260711/supplement/build_reader_facing_tables_s2_s3_20260718_235342.py",
    "ISCIENCE_REVISION_20260711/supplement/build_tables_s2_v3_s3_v6_20260719_210611.py",
    "ISCIENCE_REVISION_20260711/supplement/build_tables_s2_v4_s3_v7_20260719_213408.py",
    "ISCIENCE_REVISION_20260711/supplement/build_tables_s2_v5_s3_v8_20260719_214045.py",
    "ISCIENCE_REVISION_20260711/supplement/build_table_s2_v6_sha256_expansion_20260719_215020.py",
    "ISCIENCE_REVISION_20260711/supplement/build_supplemental_legends_v7_20260719_211018.py",
    "ISCIENCE_REVISION_20260711/supplement/build_supplemental_legends_v8_20260719_214546.py",
    "provenance/MANUAL_CODE_AUDIT_V1.0.2.md",
    "provenance/SOURCE_PROVENANCE_V1.0.2.md",
}
FORBIDDEN_BASENAMES = {
    "generate_verified_correlation_tables_20251225_120000.py",
    "rhodophyta_bicluster_analysis_20251208.py",
    "spatial_revision_analysis_20260711_080029.py",
    "earth_engine_spatial_reextraction_20260711_104304.py",
    "analyze_alphaearth_pfam_correlations_20251017.py",
    "create_pfam_alphaearth_biclustered_heatmap_20251112.py",
    "per_phylum_alphaearth_analysis_20251207.py",
    "extract_alphaearth_embeddings_20251019.py",
}
PRIVATE_PATH_PATTERNS = [
    re.compile(r"/Users/[A-Za-z0-9._-]+/"),
    re.compile(r"/home/[A-Za-z0-9._-]+/"),
    re.compile(r"[A-Za-z]:\\\\Users\\\\[^\\\\]+\\\\"),
]
CREDENTIAL_PATTERNS = [
    re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----"),
    re.compile(r"\bgh[pousr]_[A-Za-z0-9]{20,}\b"),
    re.compile(r"\bgithub_pat_[A-Za-z0-9_]{20,}\b"),
    re.compile(r"\bAIza[0-9A-Za-z_-]{30,}\b"),
    re.compile(r"\bAKIA[0-9A-Z]{16}\b"),
]
MAX_RELEASE_FILE_BYTES = 5 * 1024 * 1024
DISALLOWED_RELEASE_SUFFIXES = {
    ".csv",
    ".csv.gz",
    ".docx",
    ".fa",
    ".faa",
    ".fasta",
    ".fna",
    ".hmm",
    ".npy",
    ".npz",
    ".parquet",
    ".pdf",
    ".pkl",
    ".png",
    ".svg",
    ".tif",
    ".tiff",
    ".xlsx",
}
PROHIBITED_SYNTHETIC_API_PATTERNS = [
    re.compile(r"\btorch\.rand(?:n|int)?\s*\("),
    re.compile(r"\bnp\.random\.(?:rand|randn|random|random_sample)\s*\("),
]
GENERATED_PROVENANCE_RE = re.compile(
    r"^provenance/public_release_(?:audit|inventory)_\d{8}_\d{6}\.(?:json|tsv)$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--run-id", default=datetime.now().strftime("%Y%m%d_%H%M%S"))
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def classify(path: Path) -> str:
    relative = path.as_posix()
    if path.suffix == ".py":
        if "/integrity/" in relative or "/provenance/" in relative:
            return "integrity_or_provenance_code"
        if (
            "/supplement/" in relative
            or "/figures/" in relative
            or "/figure5_gee_primary/" in relative
        ):
            return "figure_or_table_support_code"
        return "analysis_code"
    if path.name in {"README.md", "DATA_INPUTS.md", "MANUAL_STEPS.md"}:
        return "workflow_documentation"
    if path.name in {"SOFTWARE_CITATIONS.md", "SOFTWARE_ENVIRONMENTS.md"}:
        return "software_documentation"
    if path.name in {"LICENSE", "CITATION.cff"}:
        return "release_metadata"
    return "other_documentation"


def iter_release_files(root: Path, output_names: set[str]):
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(root)
        if ".git" in relative.parts or "__pycache__" in relative.parts:
            continue
        if (
            path.name.endswith((".pyc", ".pyo"))
            or path.name in output_names
            or GENERATED_PROVENANCE_RE.fullmatch(relative.as_posix())
        ):
            continue
        yield path, relative


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    if not root.is_dir():
        raise NotADirectoryError(root)
    provenance = root / "provenance"
    provenance.mkdir(exist_ok=True)
    tsv_path = provenance / f"public_release_inventory_{args.run_id}.tsv"
    json_path = provenance / f"public_release_audit_{args.run_id}.json"
    if tsv_path.exists() or json_path.exists():
        raise FileExistsError("Refusing to overwrite an existing release audit")

    output_names = {tsv_path.name, json_path.name}
    rows: list[dict[str, object]] = []
    failures: list[str] = []
    warnings: list[str] = []
    observed: set[str] = set()
    python_count = 0

    for path, relative in iter_release_files(root, output_names):
        relative_text = relative.as_posix()
        observed.add(relative_text)
        if path.name in FORBIDDEN_BASENAMES:
            failures.append(f"forbidden superseded script: {relative_text}")
        rows.append(
            {
                "relative_path": relative_text,
                "role": classify(relative),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
        if path.stat().st_size > MAX_RELEASE_FILE_BYTES:
            failures.append(
                f"release file exceeds {MAX_RELEASE_FILE_BYTES} bytes: {relative_text}"
            )
        compound_suffix = "".join(path.suffixes).lower()
        if (
            path.suffix.lower() in DISALLOWED_RELEASE_SUFFIXES
            or compound_suffix in DISALLOWED_RELEASE_SUFFIXES
        ):
            failures.append(f"scientific data/result file included: {relative_text}")
        is_text = path.suffix.lower() in TEXT_SUFFIXES or path.name in {
            "LICENSE",
            ".gitignore",
        }
        if not is_text:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError as exc:
            failures.append(f"non-UTF-8 text file {relative_text}: {exc}")
            continue
        if path.resolve() != Path(__file__).resolve():
            for pattern in PRIVATE_PATH_PATTERNS:
                if pattern.search(text):
                    failures.append(f"private absolute path in {relative_text}")
                    break
        for pattern in CREDENTIAL_PATTERNS:
            if pattern.search(text):
                failures.append(f"credential-like string in {relative_text}")
                break
        if path.suffix == ".py":
            python_count += 1
            try:
                compile(text, relative_text, "exec")
            except SyntaxError as exc:
                failures.append(f"syntax error in {relative_text}: {exc}")
            if path.resolve() != Path(__file__).resolve() and (
                "np.random" in text or "random." in text
            ):
                warnings.append(
                    f"random API present for manual review (resampling permitted): {relative_text}"
                )
            if path.resolve() != Path(__file__).resolve() and "np.linspace" in text:
                warnings.append(
                    f"np.linspace present for manual review (plot bins/ticks permitted): {relative_text}"
                )
            for pattern in PROHIBITED_SYNTHETIC_API_PATTERNS:
                if pattern.search(text):
                    failures.append(
                        f"prohibited unscoped synthetic random API in {relative_text}"
                    )
                    break

    missing = sorted(REQUIRED_PATHS - observed)
    failures.extend(f"missing required release file: {path}" for path in missing)
    rows.sort(key=lambda row: str(row["relative_path"]))
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["relative_path", "role", "bytes", "sha256"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    audit = {
        "script_version": SCRIPT_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "root_directory_name": root.name,
        "files_in_inventory": len(rows),
        "python_files_compiled_in_memory": python_count,
        "inventory": {
            "path": tsv_path.relative_to(root).as_posix(),
            "sha256": sha256(tsv_path),
        },
        "checks": {
            "required_files_present": not missing,
            "forbidden_superseded_scripts_absent": not any(
                item.startswith("forbidden superseded script") for item in failures
            ),
            "private_absolute_paths_absent": not any(
                item.startswith("private absolute path") for item in failures
            ),
            "credential_like_strings_absent": not any(
                item.startswith("credential-like string") for item in failures
            ),
            "python_syntax_valid": not any(
                item.startswith("syntax error") for item in failures
            ),
            "no_large_or_scientific_result_files": not any(
                item.startswith("release file exceeds")
                or item.startswith("scientific data/result file included")
                for item in failures
            ),
            "prohibited_synthetic_random_apis_absent": not any(
                item.startswith("prohibited unscoped synthetic random API")
                for item in failures
            ),
        },
        "warnings": sorted(set(warnings)),
        "failures": failures,
        "passed": not failures,
    }
    json_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(tsv_path)
    print(json_path)
    print(json.dumps({"passed": audit["passed"], "failures": failures}, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

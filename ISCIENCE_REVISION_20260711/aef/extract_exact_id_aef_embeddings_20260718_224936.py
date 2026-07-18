#!/usr/bin/env python3
"""Extract one year of AlphaEarth embeddings for the exact 126-genome cohort.

The reconciled analysis manifest is the sole cohort authority. The script
selects rows marked ``safe_for_aef_pfam_analysis``, preserves each canonical
``Genome`` identifier through Earth Engine, and requires a complete A00--A63
vector for every selected genome before writing output.

The annual collection is filtered to the explicit calendar year supplied with
``--year``. This prevents the ambiguous unfiltered all-years mosaic used by the
superseded extractor. ``--validate-only`` validates the real manifest and exits
without importing Earth Engine, contacting a remote service, or writing files.

No synthetic, simulated, imputed, or hard-coded result values are produced.

Created: 2026-07-18
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCRIPT_VERSION = "2026-07-18.1"
DATASET_ID = "GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL"
EXPECTED_MANIFEST_SHA256 = (
    "5880f930192d4cbb11a7563825c853e458fb1690b01c9cf7cc83323e7541bd67"
)
EXPECTED_MANIFEST_ROWS = 131
EXPECTED_COHORT_SIZE = 126
EXPECTED_PHYLUM_COUNTS = {
    "Rhodophyta": 70,
    "Ochrophyta": 43,
    "Chlorophyta": 13,
}
AXES = [f"A{index:02d}" for index in range(64)]
REQUEST_SCALE_M = 10
TILE_SCALE = 4
MIN_YEAR = 2017
MAX_YEAR = 2024

REQUIRED_MANIFEST_COLUMNS = {
    "master_row",
    "Genome",
    "Species",
    "Phylum",
    "DD latitude",
    "DD longitude",
    "aef_present",
    "aef_id_number",
    "aef_coordinate_status",
    "safe_for_raw_pfam_analysis",
    "safe_for_aef_pfam_analysis",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        required=True,
        help="Reconciled 131-row analysis manifest pinned by SHA-256.",
    )
    parser.add_argument(
        "--year",
        type=int,
        required=True,
        help="One explicit Satellite Embedding V1 calendar year (2017--2024).",
    )
    parser.add_argument(
        "--output-parent",
        type=Path,
        default=Path("ISCIENCE_REVISION_20260711/aef"),
        help="Parent for a new timestamped run directory.",
    )
    parser.add_argument(
        "--project",
        default=None,
        help=(
            "Optional Earth Engine-enabled Google Cloud project. The value is "
            "used for initialization but is not stored in outputs."
        ),
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate the pinned manifest and exact 126-row cohort, then exit.",
    )
    return parser.parse_args()


def sha256(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_bool(value: object, field: str, genome: str) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"true", "1"}:
        return True
    if normalized in {"false", "0"}:
        return False
    raise ValueError(f"Invalid Boolean {field}={value!r} for {genome}")


def parse_finite_float(value: object, field: str, genome: str) -> float:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid numeric {field}={value!r} for {genome}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"Non-finite {field}={value!r} for {genome}")
    return parsed


def parse_integer(value: object, field: str, genome: str) -> int:
    parsed = parse_finite_float(value, field, genome)
    if not parsed.is_integer():
        raise ValueError(f"Non-integer {field}={value!r} for {genome}")
    return int(parsed)


def cohort_digest(rows: list[dict[str, Any]]) -> str:
    fields = [
        "master_row",
        "Genome",
        "Species",
        "Phylum",
        "aef_id_number",
        "DD latitude",
        "DD longitude",
    ]
    canonical = [{field: row[field] for field in fields} for row in rows]
    encoded = json.dumps(
        canonical, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def load_validated_cohort(manifest_path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    actual_hash = sha256(manifest_path)
    if actual_hash != EXPECTED_MANIFEST_SHA256:
        raise ValueError(
            "Manifest SHA-256 mismatch: "
            f"observed {actual_hash}, expected {EXPECTED_MANIFEST_SHA256}"
        )

    with manifest_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        headers = set(reader.fieldnames or [])
        missing = sorted(REQUIRED_MANIFEST_COLUMNS - headers)
        if missing:
            raise ValueError(f"Manifest is missing required columns: {missing}")
        manifest_rows = list(reader)

    if len(manifest_rows) != EXPECTED_MANIFEST_ROWS:
        raise ValueError(
            f"Expected {EXPECTED_MANIFEST_ROWS} manifest rows, found {len(manifest_rows)}"
        )
    all_genomes = [str(row["Genome"]).strip() for row in manifest_rows]
    if any(not genome for genome in all_genomes):
        raise ValueError("Manifest contains a blank Genome identifier")
    duplicate_genomes = sorted(
        genome for genome, count in Counter(all_genomes).items() if count > 1
    )
    if duplicate_genomes:
        raise ValueError(f"Manifest contains duplicate Genome IDs: {duplicate_genomes}")

    selected: list[dict[str, Any]] = []
    for source in manifest_rows:
        genome = str(source["Genome"]).strip()
        if not parse_bool(
            source["safe_for_aef_pfam_analysis"],
            "safe_for_aef_pfam_analysis",
            genome,
        ):
            continue
        if not parse_bool(source["aef_present"], "aef_present", genome):
            raise ValueError(f"Selected genome is not marked aef_present: {genome}")
        if not parse_bool(
            source["safe_for_raw_pfam_analysis"],
            "safe_for_raw_pfam_analysis",
            genome,
        ):
            raise ValueError(f"Selected genome is not safe for raw Pfam analysis: {genome}")
        if str(source["aef_coordinate_status"]).strip() != "EXACT_NUMERIC_MATCH_TO_MASTER":
            raise ValueError(f"Selected genome lacks exact coordinate authentication: {genome}")

        species = str(source["Species"]).strip()
        phylum = str(source["Phylum"]).strip()
        if not species or not phylum:
            raise ValueError(f"Selected genome has blank Species or Phylum: {genome}")
        latitude = parse_finite_float(source["DD latitude"], "DD latitude", genome)
        longitude = parse_finite_float(source["DD longitude"], "DD longitude", genome)
        if not -90 <= latitude <= 90:
            raise ValueError(f"Latitude outside [-90, 90] for {genome}: {latitude}")
        if not -180 <= longitude <= 180:
            raise ValueError(f"Longitude outside [-180, 180] for {genome}: {longitude}")

        selected.append(
            {
                "master_row": parse_integer(source["master_row"], "master_row", genome),
                "Genome": genome,
                "Species": species,
                "Phylum": phylum,
                "aef_id_number": parse_integer(
                    source["aef_id_number"], "aef_id_number", genome
                ),
                "DD latitude": latitude,
                "DD longitude": longitude,
            }
        )

    selected.sort(key=lambda row: row["master_row"])
    if len(selected) != EXPECTED_COHORT_SIZE:
        raise ValueError(
            f"Expected {EXPECTED_COHORT_SIZE} selected genomes, found {len(selected)}"
        )
    for field in ("master_row", "Genome", "aef_id_number"):
        values = [row[field] for row in selected]
        if len(values) != len(set(values)):
            raise ValueError(f"Selected cohort contains duplicate {field} values")
    observed_phyla = dict(Counter(row["Phylum"] for row in selected))
    if observed_phyla != EXPECTED_PHYLUM_COUNTS:
        raise ValueError(
            f"Unexpected selected phylum composition: {observed_phyla}; "
            f"expected {EXPECTED_PHYLUM_COUNTS}"
        )

    audit = {
        "status": "validated",
        "script_version": SCRIPT_VERSION,
        "manifest_file": manifest_path.name,
        "manifest_sha256": actual_hash,
        "manifest_rows": len(manifest_rows),
        "selected_genomes": len(selected),
        "selected_coordinate_sites": len(
            {(row["DD latitude"], row["DD longitude"]) for row in selected}
        ),
        "selected_phylum_counts": observed_phyla,
        "cohort_sha256": cohort_digest(selected),
        "cohort_filter": "safe_for_aef_pfam_analysis == true",
    }
    return selected, audit


def initialize_earth_engine(project: str | None):
    try:
        import ee
    except ImportError as exc:
        raise RuntimeError(
            "The Earth Engine Python API is required for extraction. "
            "Install earthengine-api or use --validate-only."
        ) from exc
    try:
        if project:
            ee.Initialize(project=project)
        else:
            ee.Initialize()
    except Exception as exc:
        raise RuntimeError(
            "Earth Engine authentication and an enabled Cloud project are required. "
            "No extraction output was written."
        ) from exc
    return ee


def extract_embeddings(
    cohort: list[dict[str, Any]], year: int, project: str | None
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    ee = initialize_earth_engine(project)
    start_date = f"{year}-01-01"
    end_date_exclusive = f"{year + 1}-01-01"
    collection = (
        ee.ImageCollection(DATASET_ID)
        .filterDate(start_date, end_date_exclusive)
        .select(AXES)
    )
    image_count = int(collection.size().getInfo())
    if image_count < 1:
        raise RuntimeError(f"No {DATASET_ID} images were returned for {year}")
    annual_image = collection.mosaic().select(AXES)
    returned_axes = annual_image.bandNames().getInfo()
    if returned_axes != AXES:
        raise RuntimeError(f"Unexpected AEF band order: {returned_axes}")

    features = []
    for input_row, row in enumerate(cohort, start=1):
        properties = {
            "input_row": input_row,
            "master_row": row["master_row"],
            "genome_id": row["Genome"],
            "species": row["Species"],
            "phylum": row["Phylum"],
            "aef_id_number": row["aef_id_number"],
            "input_latitude": row["DD latitude"],
            "input_longitude": row["DD longitude"],
            "embedding_year": year,
        }
        point = ee.Geometry.Point([row["DD longitude"], row["DD latitude"]])
        features.append(ee.Feature(point, properties))

    sampled = annual_image.reduceRegions(
        collection=ee.FeatureCollection(features),
        reducer=ee.Reducer.mean(),
        scale=REQUEST_SCALE_M,
        tileScale=TILE_SCALE,
    )
    response = sampled.getInfo()
    returned_features = response.get("features", [])
    if len(returned_features) != EXPECTED_COHORT_SIZE:
        raise RuntimeError(
            f"Earth Engine returned {len(returned_features)} features; "
            f"expected {EXPECTED_COHORT_SIZE}"
        )

    properties_by_id: dict[str, dict[str, Any]] = {}
    for feature in returned_features:
        properties = feature.get("properties", {})
        genome = str(properties.get("genome_id", "")).strip()
        if not genome or genome in properties_by_id:
            raise RuntimeError(f"Missing or duplicate returned genome ID: {genome!r}")
        properties_by_id[genome] = properties

    expected_ids = {row["Genome"] for row in cohort}
    returned_ids = set(properties_by_id)
    if returned_ids != expected_ids:
        raise RuntimeError(
            "Returned exact-ID set mismatch: "
            f"missing={sorted(expected_ids-returned_ids)}, "
            f"extra={sorted(returned_ids-expected_ids)}"
        )

    output_rows: list[dict[str, Any]] = []
    for source in cohort:
        genome = source["Genome"]
        properties = properties_by_id[genome]
        for property_name, expected in (
            ("master_row", source["master_row"]),
            ("aef_id_number", source["aef_id_number"]),
            ("embedding_year", year),
        ):
            observed = parse_integer(properties.get(property_name), property_name, genome)
            if observed != expected:
                raise RuntimeError(
                    f"Returned {property_name} changed for {genome}: {observed} != {expected}"
                )
        for property_name, expected in (
            ("input_latitude", source["DD latitude"]),
            ("input_longitude", source["DD longitude"]),
        ):
            observed = parse_finite_float(
                properties.get(property_name), property_name, genome
            )
            if not math.isclose(observed, expected, rel_tol=0, abs_tol=1e-12):
                raise RuntimeError(
                    f"Returned {property_name} changed for {genome}: {observed} != {expected}"
                )

        output = {
            "Genome": genome,
            "Species": source["Species"],
            "Phylum": source["Phylum"],
            "ID number": source["aef_id_number"],
            "manifest_master_row": source["master_row"],
            "DD latitude": source["DD latitude"],
            "DD longitude": source["DD longitude"],
            "embedding_year": year,
            "dataset_id": DATASET_ID,
        }
        for axis in AXES:
            output[axis] = parse_finite_float(properties.get(axis), axis, genome)
        output_rows.append(output)

    extraction_audit = {
        "dataset_id": DATASET_ID,
        "embedding_year": year,
        "filter_start": start_date,
        "filter_end_exclusive": end_date_exclusive,
        "collection_image_count": image_count,
        "requested_scale_m": REQUEST_SCALE_M,
        "tile_scale": TILE_SCALE,
        "reducer": "mean",
        "returned_genomes": len(output_rows),
        "returned_axes": len(AXES),
        "earth_engine_project_supplied": bool(project),
        "earthengine_api_version": getattr(ee, "__version__", "unknown"),
    }
    return output_rows, extraction_audit


def write_outputs(
    rows: list[dict[str, Any]],
    cohort_audit: dict[str, Any],
    extraction_audit: dict[str, Any],
    output_parent: Path,
) -> tuple[Path, Path]:
    generated = datetime.now(timezone.utc)
    timestamp = generated.strftime("%Y%m%d_%H%M%S")
    run_dir = output_parent / f"exact_id_aef_extraction_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=False)
    csv_path = run_dir / "exact_id_aef_embeddings.csv"
    manifest_path = run_dir / "run_manifest.json"
    fields = [
        "Genome",
        "Species",
        "Phylum",
        "ID number",
        "manifest_master_row",
        "DD latitude",
        "DD longitude",
        "embedding_year",
        "dataset_id",
        *AXES,
    ]
    with csv_path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    script_path = Path(__file__).resolve()
    manifest = {
        "status": "completed",
        "generated_utc": generated.isoformat(),
        "purpose": "year-explicit exact-ID AEF extraction for the 126-genome cohort",
        "script": script_path.name,
        "script_version": SCRIPT_VERSION,
        "script_sha256": sha256(script_path),
        "cohort_validation": cohort_audit,
        "extraction": extraction_audit,
        "output": {
            "file": csv_path.name,
            "sha256": sha256(csv_path),
            "rows": len(rows),
            "columns": len(fields),
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "earthengine_api": extraction_audit["earthengine_api_version"],
        },
        "scientific_integrity": {
            "synthetic_or_simulated_values": False,
            "imputed_embedding_values": False,
            "all_returned_genome_ids_exactly_validated": True,
            "all_aef_axes_complete_and_finite": True,
        },
    }
    with manifest_path.open("x", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return csv_path, manifest_path


def main() -> int:
    args = parse_args()
    if not MIN_YEAR <= args.year <= MAX_YEAR:
        raise ValueError(f"--year must be between {MIN_YEAR} and {MAX_YEAR}")
    cohort, cohort_audit = load_validated_cohort(args.manifest)
    if args.validate_only:
        validation = {
            **cohort_audit,
            "status": "validated_only",
            "embedding_year": args.year,
            "dataset_id": DATASET_ID,
            "files_written": 0,
        }
        print(json.dumps(validation, indent=2, ensure_ascii=False))
        return 0

    rows, extraction_audit = extract_embeddings(cohort, args.year, args.project)
    csv_path, manifest_path = write_outputs(
        rows, cohort_audit, extraction_audit, args.output_parent
    )
    print(f"Extracted genomes: {len(rows)}")
    print(f"Embedding axes: {len(AXES)}")
    print(f"Results: {csv_path.resolve()}")
    print(f"Manifest: {manifest_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

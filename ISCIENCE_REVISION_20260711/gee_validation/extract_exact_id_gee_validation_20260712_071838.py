#!/usr/bin/env python3
"""Re-extract the original 13 GEE variables for the exact 126-genome cohort.

This is a bounded validation/correction of the submitted workflow, not a new
discovery analysis. It reproduces the retained JavaScript image collections,
date filters, reducers, combined stack, and 1-km reduceRegions request, while
carrying the canonical genome ID through every row. No values are simulated,
imputed, or reconstructed from summary statistics.

Generated: 2026-07-12
Requires: authenticated Google Earth Engine Python API.
"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import ee


ROOT = Path(__file__).resolve().parents[2]
INPUT = (
    ROOT
    / "ISCIENCE_REVISION_20260711"
    / "spatial"
    / "coordinate_confidence_audit_20260711_105047.csv"
)
OUTDIR = ROOT / "ISCIENCE_REVISION_20260711" / "gee_validation"
OUTPUT = OUTDIR / "exact_id_gee_environmental_extraction_20260712_071838.csv"
MANIFEST = OUTDIR / "exact_id_gee_environmental_extraction_manifest_20260712_071838.json"

EXPECTED_INPUT_SHA256 = "36ab856da339da891b279bf7fffc9cbc8ccb91de37cb273d2beb9c8ba8f73da1"
START_DATE = "2020-01-01"
END_DATE_EXCLUSIVE = "2023-12-31"  # exact retained-script filterDate end
REQUESTED_SCALE_M = 1000
TILE_SCALE = 4

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
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_cohort() -> list[dict[str, str]]:
    if not INPUT.is_file():
        raise FileNotFoundError(INPUT)
    actual = sha256(INPUT)
    if actual != EXPECTED_INPUT_SHA256:
        raise RuntimeError(f"Cohort hash mismatch: {actual} != {EXPECTED_INPUT_SHA256}")
    with INPUT.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    ids = [row["genome_id"] for row in rows]
    if len(rows) != 126 or len(set(ids)) != 126:
        raise RuntimeError("Input is not the exact 126-row unique-genome AEF cohort")
    for row in rows:
        float(row["latitude"])
        float(row["longitude"])
        if row["clean_metadata_match"].lower() != "true":
            raise RuntimeError(f"Metadata mismatch flag for {row['genome_id']}")
    return rows


def build_stack() -> tuple[ee.Image, dict[str, int]]:
    sst_collection = (
        ee.ImageCollection("NASA/OCEANDATA/MODIS-Aqua/L3SMI")
        .select(["sst"])
        .filterDate(START_DATE, END_DATE_EXCLUSIVE)
    )
    sst_mean = sst_collection.mean().rename("sst_mean_c")
    sst_max = sst_collection.max().rename("sst_max_c")
    sst_min = sst_collection.min().rename("sst_min_c")
    sst_range = sst_max.subtract(sst_min).rename("sst_annual_range_c")
    sst_summer = (
        sst_collection.filter(ee.Filter.calendarRange(6, 8, "month"))
        .mean()
        .rename("sst_summer_c")
    )
    sst_winter = (
        sst_collection.filter(ee.Filter.calendarRange(12, 2, "month"))
        .mean()
        .rename("sst_winter_c")
    )

    chlorophyll_collection = (
        ee.ImageCollection("NASA/OCEANDATA/MODIS-Aqua/L3SMI")
        .select(["chlor_a"])
        .filterDate(START_DATE, END_DATE_EXCLUSIVE)
    )
    chlorophyll_mean = chlorophyll_collection.mean().rename("chlorophyll_mean_mg_m3")
    chlorophyll_max = chlorophyll_collection.max().rename("chlorophyll_max_mg_m3")
    chlorophyll_std = chlorophyll_collection.reduce(ee.Reducer.stdDev()).rename(
        "chlorophyll_std_mg_m3"
    )

    poc_collection = (
        ee.ImageCollection("NASA/OCEANDATA/MODIS-Aqua/L3SMI")
        .select(["poc"])
        .filterDate(START_DATE, END_DATE_EXCLUSIVE)
    )
    poc_mean = poc_collection.mean().rename("poc_mean_mg_m3")

    etopo = ee.Image("NOAA/NGDC/ETOPO1").select("bedrock")
    bathymetry = etopo.multiply(-1).rename("depth_meters")
    land = etopo.gt(0)
    distance_to_coast = (
        land.fastDistanceTransform()
        .sqrt()
        .multiply(ee.Image.pixelArea().sqrt())
        .divide(1000)
        .rename("distance_coast_km")
    )

    rrs_collection = (
        ee.ImageCollection("NASA/OCEANDATA/MODIS-Aqua/L3SMI")
        .select(["Rrs_443", "Rrs_555"])
        .filterDate(START_DATE, END_DATE_EXCLUSIVE)
    )
    rrs_mean = rrs_collection.mean()
    water_clarity = (
        rrs_mean.select("Rrs_443")
        .divide(rrs_mean.select("Rrs_555"))
        .rename("water_clarity_ratio")
    )

    stack = (
        sst_mean.addBands(sst_max)
        .addBands(sst_min)
        .addBands(sst_range)
        .addBands(sst_summer)
        .addBands(sst_winter)
        .addBands(chlorophyll_mean)
        .addBands(chlorophyll_max)
        .addBands(chlorophyll_std)
        .addBands(poc_mean)
        .addBands(bathymetry)
        .addBands(distance_to_coast)
        .addBands(water_clarity)
    )
    image_counts = {
        "sst": int(sst_collection.size().getInfo()),
        "chlorophyll": int(chlorophyll_collection.size().getInfo()),
        "poc": int(poc_collection.size().getInfo()),
        "rrs": int(rrs_collection.size().getInfo()),
    }
    return stack, image_counts


def build_features(rows: list[dict[str, str]]) -> ee.FeatureCollection:
    features = []
    for index, row in enumerate(rows, start=1):
        latitude = float(row["latitude"])
        longitude = float(row["longitude"])
        properties = {
            "input_row": index,
            "genome_id": row["genome_id"],
            "species": row["species"],
            "phylum": row["phylum"],
            "input_latitude": latitude,
            "input_longitude": longitude,
        }
        features.append(ee.Feature(ee.Geometry.Point([longitude, latitude]), properties))
    return ee.FeatureCollection(features)


def materialize(rows: list[dict[str, str]], collection_info: dict) -> list[dict[str, object]]:
    features = collection_info.get("features", [])
    by_id: dict[str, dict] = {}
    for feature in features:
        properties = feature.get("properties", {})
        genome_id = properties.get("genome_id")
        if not genome_id or genome_id in by_id:
            raise RuntimeError(f"Missing or duplicate returned genome ID: {genome_id}")
        by_id[genome_id] = properties
    expected_ids = {row["genome_id"] for row in rows}
    if set(by_id) != expected_ids:
        raise RuntimeError(
            f"Returned ID set mismatch: missing={sorted(expected_ids-set(by_id))}, "
            f"extra={sorted(set(by_id)-expected_ids)}"
        )

    output_rows: list[dict[str, object]] = []
    for index, source in enumerate(rows, start=1):
        props = by_id[source["genome_id"]]
        returned_lat = float(props["input_latitude"])
        returned_lon = float(props["input_longitude"])
        if abs(returned_lat - float(source["latitude"])) > 1e-10:
            raise RuntimeError(f"Latitude changed for {source['genome_id']}")
        if abs(returned_lon - float(source["longitude"])) > 1e-10:
            raise RuntimeError(f"Longitude changed for {source['genome_id']}")
        out = {
            "input_row": index,
            "genome_id": source["genome_id"],
            "species": source["species"],
            "phylum": source["phylum"],
            "latitude": source["latitude"],
            "longitude": source["longitude"],
        }
        for column in ENVIRONMENTAL_COLUMNS:
            out[column] = props.get(column, "")
        output_rows.append(out)
    return output_rows


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    cohort = read_cohort()
    started = datetime.now(timezone.utc)
    ee.Initialize()
    stack, image_counts = build_stack()
    band_names = stack.bandNames().getInfo()
    if band_names != ENVIRONMENTAL_COLUMNS:
        raise RuntimeError(f"Unexpected band order: {band_names}")

    features = build_features(cohort)
    extracted = stack.reduceRegions(
        collection=features,
        reducer=ee.Reducer.first(),
        scale=REQUESTED_SCALE_M,
        tileScale=TILE_SCALE,
    )
    collection_info = extracted.getInfo()
    output_rows = materialize(cohort, collection_info)

    fields = [
        "input_row",
        "genome_id",
        "species",
        "phylum",
        "latitude",
        "longitude",
        *ENVIRONMENTAL_COLUMNS,
    ]
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)

    # Re-read the actual CSV and compute only audit summaries from materialized values.
    with OUTPUT.open(newline="", encoding="utf-8") as handle:
        check_rows = list(csv.DictReader(handle))
    if len(check_rows) != 126 or len({row["genome_id"] for row in check_rows}) != 126:
        raise RuntimeError("Materialized extraction is not 126 unique genome IDs")
    nonmissing = {
        column: sum(row[column] not in {"", "NA", "nan", "None"} for row in check_rows)
        for column in ENVIRONMENTAL_COLUMNS
    }
    missing_patterns = Counter(
        tuple(column for column in ENVIRONMENTAL_COLUMNS if row[column] == "")
        for row in check_rows
    )
    completed = datetime.now(timezone.utc)
    manifest = {
        "purpose": "bounded exact-genome-ID validation/correction of the submitted GEE extraction",
        "not_a_new_discovery_analysis": True,
        "generated_by": str(Path(__file__).relative_to(ROOT)),
        "generator_sha256": sha256(Path(__file__)),
        "started_utc": started.isoformat(),
        "completed_utc": completed.isoformat(),
        "runtime_seconds": (completed - started).total_seconds(),
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "earthengine_api": ee.__version__,
        },
        "input": {
            "path": str(INPUT.relative_to(ROOT)),
            "sha256": sha256(INPUT),
            "rows": len(cohort),
            "unique_genome_ids": len({row["genome_id"] for row in cohort}),
            "unique_coordinate_sites": len(
                {(row["latitude"], row["longitude"]) for row in cohort}
            ),
        },
        "earth_engine_specification": {
            "ocean_collection": "NASA/OCEANDATA/MODIS-Aqua/L3SMI",
            "bathymetry_image": "NOAA/NGDC/ETOPO1",
            "start_date_inclusive": START_DATE,
            "end_date_exclusive": END_DATE_EXCLUSIVE,
            "requested_reduce_regions_scale_m": REQUESTED_SCALE_M,
            "requested_scale_does_not_change_native_product_resolution": True,
            "tile_scale": TILE_SCALE,
            "reducer": "ee.Reducer.first()",
            "image_counts": image_counts,
            "bands": ENVIRONMENTAL_COLUMNS,
            "distance_to_coast_formula": "sqrt(fastDistanceTransform(ETOPO1 bedrock > 0)) * sqrt(pixelArea) / 1000",
        },
        "output": {
            "path": str(OUTPUT.relative_to(ROOT)),
            "sha256": sha256(OUTPUT),
            "rows": len(check_rows),
            "unique_genome_ids": len({row["genome_id"] for row in check_rows}),
            "nonmissing_by_variable": nonmissing,
            "missing_pattern_counts": {
                "|".join(pattern) if pattern else "COMPLETE_13": count
                for pattern, count in sorted(missing_patterns.items(), key=lambda item: str(item[0]))
            },
        },
        "integrity_checks": {
            "input_hash_matches_expected": True,
            "input_126_unique_ids": True,
            "returned_id_set_exact": True,
            "coordinates_unchanged_at_1e-10": True,
            "no_species_name_join": True,
            "no_imputation": True,
            "no_synthetic_values": True,
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), "manifest": str(MANIFEST), "nonmissing": nonmissing}, indent=2))


if __name__ == "__main__":
    main()

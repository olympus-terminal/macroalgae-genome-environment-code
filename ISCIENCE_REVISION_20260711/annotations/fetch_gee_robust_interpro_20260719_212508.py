#!/usr/bin/env python3
"""Fetch InterPro/Pfam records for the 49 all-seven primary GEE pairs.

Accessions are read from the authenticated GEE candidate table rather than
entered as analysis results.  Missing API annotations remain explicitly
unavailable.  Raw API responses and their hashes are retained for provenance.

Created: 2026-07-19.
"""

from __future__ import annotations

import csv
import hashlib
import json
import time
from datetime import datetime
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


RUN_TAG = "20260719_212508"
ROOT = Path(__file__).resolve().parents[2]
CANDIDATES = (
    ROOT
    / "ISCIENCE_REVISION_20260711/gee_validation/"
    "GEE_primary_selected84_candidate_summary_refined99999_20260719_210049.csv"
)
CANDIDATES_SHA256 = "85bd71ff3128219949e74b168a8714f9b65f2020b5854a0ee4ff0f05059d2e76"
OUTDIR = Path(__file__).resolve().parent
CSV_OUT = OUTDIR / f"GEE_all_seven_InterPro_annotations_{RUN_TAG}.csv"
JSON_OUT = OUTDIR / f"GEE_all_seven_InterPro_annotations_{RUN_TAG}.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def read_accessions() -> tuple[list[str], int]:
    if sha256_file(CANDIDATES) != CANDIDATES_SHA256:
        raise RuntimeError(f"Authenticated candidate input changed: {CANDIDATES}")
    robust_rows = []
    with CANDIDATES.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            value = row["robust_candidate_all_required_checks"].strip().lower()
            if value in {"true", "1"}:
                robust_rows.append(row)
            elif value not in {"false", "0"}:
                raise RuntimeError(f"Unexpected robust-candidate flag: {value!r}")
    if len(robust_rows) != 49:
        raise RuntimeError(f"Expected 49 all-seven rows; observed {len(robust_rows)}")
    accessions = sorted({row["pfam"] for row in robust_rows})
    if len(accessions) != 30:
        raise RuntimeError(f"Expected 30 distinct all-seven Pfams; observed {len(accessions)}")
    return accessions, len(robust_rows)


def fetch(accession: str) -> tuple[bytes, int]:
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}/"
    request = Request(url, headers={"User-Agent": "Macroalgae-iScience-revision/1.0"})
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            with urlopen(request, timeout=60) as response:
                return response.read(), response.status
        except (HTTPError, URLError, TimeoutError) as error:
            last_error = error
            if attempt < 2:
                time.sleep(2**attempt)
    raise RuntimeError(f"InterPro request failed for {accession}: {last_error}")


def main() -> None:
    for path in (CSV_OUT, JSON_OUT):
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite {path}")

    accessions, robust_pair_count = read_accessions()
    records = []
    raw_responses = {}
    for accession in accessions:
        url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}/"
        payload, status = fetch(accession)
        if status != 200:
            raise RuntimeError(f"InterPro returned HTTP {status} for {accession}")
        document = json.loads(payload)
        metadata = document.get("metadata")
        if not isinstance(metadata, dict) or metadata.get("accession") != accession:
            raise RuntimeError(f"Unexpected InterPro payload for {accession}")
        name = metadata.get("name") or {}
        records.append(
            {
                "accession": accession,
                "short_name": name.get("short") or "annotation not available",
                "name": name.get("name") or "annotation not available",
                "entry_type": metadata.get("type") or "annotation not available",
                "integrated_interpro": metadata.get("integrated") or "annotation not available",
                "source_database": metadata.get("source_database") or "annotation not available",
                "api_url": url,
                "http_status": status,
                "response_sha256": sha256_bytes(payload),
            }
        )
        raw_responses[accession] = document

    with CSV_OUT.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    JSON_OUT.write_text(
        json.dumps(
            {
                "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
                "program": str(Path(__file__).resolve()),
                "program_sha256": sha256_file(Path(__file__).resolve()),
                "candidate_input": str(CANDIDATES.resolve()),
                "candidate_input_sha256": CANDIDATES_SHA256,
                "robust_pair_count": robust_pair_count,
                "distinct_pfam_count": len(accessions),
                "random_or_synthetic_data_used": False,
                "records": records,
                "raw_api_responses": raw_responses,
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    print(CSV_OUT.resolve())
    print(JSON_OUT.resolve())


if __name__ == "__main__":
    main()

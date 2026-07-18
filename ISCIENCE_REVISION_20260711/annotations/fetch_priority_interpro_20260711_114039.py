#!/usr/bin/env python3
"""Fetch verified InterPro/Pfam names for manuscript priority accessions.

The program writes only responses returned by the public InterPro API. It does
not infer missing annotations and refuses to overwrite an existing output.
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path
from urllib.request import Request, urlopen


ACCESSIONS = (
    "PF00092",
    "PF00223",
    "PF00346",
    "PF01513",
    "PF01638",
    "PF05605",
    "PF10988",
    "PF12094",
    "PF13411",
    "PF14317",
    "PF14450",
    "PF16087",
)
RUN_TAG = "20260711_114039"
OUTDIR = Path(__file__).resolve().parent
CSV_OUT = OUTDIR / f"priority_interpro_annotations_{RUN_TAG}.csv"
JSON_OUT = OUTDIR / f"priority_interpro_annotations_{RUN_TAG}.json"


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    for path in (CSV_OUT, JSON_OUT):
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite {path}")

    records = []
    raw_responses = {}
    for accession in ACCESSIONS:
        url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}/"
        request = Request(url, headers={"User-Agent": "Macroalgae-iScience-revision/1.0"})
        with urlopen(request, timeout=60) as response:
            payload = response.read()
            status = response.status
        if status != 200:
            raise RuntimeError(f"InterPro returned HTTP {status} for {accession}")
        document = json.loads(payload)
        metadata = document.get("metadata")
        if not isinstance(metadata, dict) or metadata.get("accession") != accession:
            raise RuntimeError(f"Unexpected InterPro payload for {accession}")
        name = metadata.get("name") or {}
        record = {
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
        records.append(record)
        raw_responses[accession] = document

    generated_at = datetime.now().astimezone().isoformat(timespec="seconds")
    with CSV_OUT.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    JSON_OUT.write_text(
        json.dumps(
            {
                "generated_at": generated_at,
                "program": str(Path(__file__).resolve()),
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

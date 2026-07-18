#!/usr/bin/env python3
"""Audit retained HMMER output headers used by the canonical Pfam reconstruction.

This script reads only the 131 real tblout paths in the canonical raw-source
manifest. It writes one row per genome plus a JSON summary and refuses to
overwrite either output.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUN_ID = "20260711_123302"
INPUT = ROOT / "ISCIENCE_REVISION_20260711/integrity/reconciled_analysis_manifest_20260711_110650.csv"
OUTPUT_CSV = ROOT / f"ISCIENCE_REVISION_20260711/integrity/hmm_header_audit_{RUN_ID}.csv"
OUTPUT_JSON = ROOT / f"ISCIENCE_REVISION_20260711/integrity/hmm_header_audit_{RUN_ID}.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_header(path: Path) -> dict[str, str | bool]:
    fields: dict[str, str | bool] = {
        "program": "",
        "version": "",
        "pipeline_mode": "",
        "query_file": "",
        "option_settings": "",
        "sequence_e_threshold": "",
    }
    patterns = {
        "program": re.compile(r"#\s+Program:\s*(.+?)\s*$"),
        "version": re.compile(r"#\s+Version:\s*([^\s]+)"),
        "pipeline_mode": re.compile(r"#\s+Pipeline\s+mode:\s*([^\s]+)"),
        "query_file": re.compile(r"#\s+Query\s+file:\s*(.+?)\s*$"),
        "option_settings": re.compile(r"#\s+Option\s+settings:\s*(.+?)\s*$"),
    }
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            normalized = " ".join(line.replace("\t", " ").split())
            for key, pattern in patterns.items():
                match = pattern.search(normalized)
                if match:
                    fields[key] = match.group(1)
    option = str(fields["option_settings"])
    threshold = re.search(r"(?:^|\s)-E\s+([^\s]+)", option)
    fields["sequence_e_threshold"] = threshold.group(1) if threshold else ""
    fields["query_basename"] = Path(str(fields["query_file"])).name
    fields["header_complete"] = all(
        fields[key] for key in ("program", "version", "pipeline_mode", "query_file", "option_settings")
    )
    return fields


def main() -> None:
    for output in (OUTPUT_CSV, OUTPUT_JSON):
        if output.exists():
            raise FileExistsError(f"Refusing to overwrite {output}")

    rows: list[dict[str, str | bool]] = []
    with INPUT.open(newline="", encoding="utf-8") as handle:
        for source in csv.DictReader(handle):
            tblout = ROOT / source["pfam_source_path"]
            if not tblout.is_file():
                raise FileNotFoundError(tblout)
            parsed = parse_header(tblout)
            rows.append(
                {
                    "Genome": source["Genome"],
                    "raw_tblout_realpath": str(tblout),
                    "raw_tblout_sha256": source["pfam_source_sha256"],
                    **parsed,
                }
            )

    if len(rows) != 131:
        raise ValueError(f"Expected 131 source rows, observed {len(rows)}")
    if not all(bool(row["header_complete"]) for row in rows):
        missing = [str(row["Genome"]) for row in rows if not bool(row["header_complete"])]
        raise ValueError(f"Incomplete HMMER headers: {missing}")

    fieldnames = list(rows[0])
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "run_id": RUN_ID,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "program": str(Path(__file__).resolve()),
        "program_sha256": sha256(Path(__file__).resolve()),
        "input_manifest": str(INPUT),
        "input_manifest_sha256": sha256(INPUT),
        "n_genomes": len(rows),
        "program_counts": dict(sorted(Counter(str(row["program"]) for row in rows).items())),
        "version_counts": dict(sorted(Counter(str(row["version"]) for row in rows).items())),
        "pipeline_mode_counts": dict(sorted(Counter(str(row["pipeline_mode"]) for row in rows).items())),
        "query_basename_counts": dict(sorted(Counter(str(row["query_basename"]) for row in rows).items())),
        "sequence_e_threshold_counts": dict(
            sorted(Counter(str(row["sequence_e_threshold"]) for row in rows).items())
        ),
        "output_csv": str(OUTPUT_CSV),
        "output_csv_sha256": sha256(OUTPUT_CSV),
        "synthetic_or_hardcoded_results": False,
    }
    with OUTPUT_JSON.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

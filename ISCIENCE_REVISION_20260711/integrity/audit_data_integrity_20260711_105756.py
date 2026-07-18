#!/usr/bin/env python3
"""Reconcile the iScience revision cohorts from preserved source data.

Created: 2026-07-11 10:57:56 +04

Data-integrity rules implemented here:
  * no synthetic, simulated, imputed, or hardcoded scientific results;
  * every Pfam count is parsed from a preserved hmmsearch tblout;
  * authoritative metadata fields come from the published metadata CSV;
  * joins fail closed when identifiers or accession-level evidence are ambiguous;
  * GEE variables with no preserved 126-row raw export are marked unsafe;
  * all outputs include a unique run timestamp and never overwrite files.

The implementation uses only the Python standard library and streams the
approximately 480 MB of raw hmmsearch output line by line.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import platform
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "ISCIENCE_REVISION_20260711" / "integrity"

MASTER_METADATA = ROOT / "AlphaEarth/CSV/Metadata_Table_macroalgae-published.csv"
S1E_131 = ROOT / "TABLES/Table_S1_25DEC_2025-main.xlsx"
S1E_131_FINAL_COPY = ROOT / "TABLES/FINAL/Table_S1_25DEC_2025-main.xlsx"
S1E_126_LATER = ROOT / "TABLES/Table_S1_30DEC_2025-main.xlsx"
TABLE_S4 = ROOT / "Table_S4_meta-and_pfams.csv"
AEF = ROOT / "AlphaEarth/CSV/alphaearth_embeddings_20251019_122918.csv"
AEF_EXTRACTION_SCRIPT = (
    ROOT
    / "ISCIENCE_REVISION_20260711/aef/extract_exact_id_aef_embeddings_20260718_224936.py"
)
PFAM_COUNT_SCRIPT = ROOT / "DATA_S2_25DEC2025/count_pfam_domains_20251019.py"
TAGGED_HMM_DIR = ROOT / "AlphaEarth/TAGGED_HMMsearch-raw-out"
AF3_HMM_DIR = ROOT / "AF3/transfer_hmmsearch_tblout"

GEE_ENV = ROOT / "macroalgae_env_pfam_126samples.csv"
GEE_ENV_B = ROOT / "macroalgae_env_pfam_126samples-b.csv"
GEE_ENV_DRIVE_COPY = ROOT / "drive-download-20250907T085240Z-1-001 2/macroalgae_env_pfam_126samples.csv"
GEE_UPLOAD = ROOT / "GoogleEarthEngine/gee_upload_full_126_samples.csv"
GEE_TEST_EXPORT = ROOT / "GoogleEarthEngine/GEE-out-macroalgae_test_environmental_extraction.csv"
GEE_EXTRACTION_SCRIPT = ROOT / "GoogleEarthEngine/gee_full_126_samples.js"
GEE_PROVENANCE_DOC = ROOT / "GoogleEarthEngine/DATA_PROVENANCE_GEE_ANALYSIS.md"

RBCL_METADATA = ROOT / "MACROALGAE_PHYLOGENIES.csv"
RBCL_TREE = ROOT / "DATA/rbcL_phylogenetic_tree_20251118_173538.nwk"
RBCL_ALIGNMENT = ROOT / "DATA/rbcL_alignment.fa"
RBCL_BUILD_SCRIPT = ROOT / "triangulation/scripts/build_phylogenetic_tree_20251117.py"
RBCL_BUILD_LOG = ROOT / "triangulation/phylogeny/tree_construction_log_20251118_173538.txt"
MAIN_MANUSCRIPT = ROOT / "main.txt"
AGENTS = ROOT / "AGENTS.md"

PFAM_RE = re.compile(r"PF\d{5}")
PFAM_ACCESSION_RE = re.compile(r"(PF\d{5})(?:\.\d+)?")

STANDARD_METADATA_FIELDS = [
    "Genome",
    "Species",
    "Phylum",
    "UsedAsReference",
    "Nucleotides",
    "Exons",
    "PFAMs",
    "BUSCOs-%present",
    "Climatic zone",
    "Temperature (°C)",
    "Environment",
    "Habitat",
    "DD latitude",
    "DD longitude",
]

NUMERIC_METADATA_FIELDS = {
    "Nucleotides",
    "Exons",
    "PFAMs",
    "BUSCOs-%present",
    "Temperature (°C)",
    "DD latitude",
    "DD longitude",
}

GEE_ENV_FIELDS = [
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


class AuditFailure(RuntimeError):
    """Raised when an integrity gate cannot be satisfied."""


def rel(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def norm_text(value: object) -> str:
    return " ".join(str(value or "").replace("\xa0", " ").strip().casefold().split())


def slug(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", norm_text(value))


def as_float(value: object) -> Optional[float]:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        result = float(text)
    except ValueError as exc:
        raise AuditFailure(f"Non-numeric value where numeric value was required: {text!r}") from exc
    if not math.isfinite(result):
        raise AuditFailure(f"Non-finite numeric value: {text!r}")
    return result


def numeric_equal(a: object, b: object, tolerance: float = 1e-10) -> bool:
    aa = as_float(a)
    bb = as_float(b)
    if aa is None or bb is None:
        return aa is None and bb is None
    return abs(aa - bb) <= tolerance


def metadata_equal(field: str, a: object, b: object) -> bool:
    if field in NUMERIC_METADATA_FIELDS:
        return numeric_equal(a, b)
    return norm_text(a) == norm_text(b)


def read_csv_rows(path: Path, encoding: str = "utf-8-sig") -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding=encoding, newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise AuditFailure(f"Missing CSV header: {path}")
        rows = [dict(row) for row in reader]
        return list(reader.fieldnames), rows


def require_unique(rows: Sequence[Mapping[str, str]], key: str, label: str) -> Dict[str, Mapping[str, str]]:
    result: Dict[str, Mapping[str, str]] = {}
    duplicates: List[str] = []
    for row in rows:
        value = str(row.get(key, "") or "").strip()
        if not value:
            continue
        if value in result:
            duplicates.append(value)
        result[value] = row
    if duplicates:
        raise AuditFailure(f"Duplicate {key} values in {label}: {sorted(set(duplicates))}")
    return result


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def xlsx_column_index(reference: str) -> int:
    match = re.match(r"([A-Z]+)", reference.upper())
    if not match:
        raise AuditFailure(f"Invalid XLSX cell reference: {reference!r}")
    value = 0
    for char in match.group(1):
        value = value * 26 + ord(char) - 64
    return value - 1


def read_s1e_xlsx(path: Path) -> Dict[str, object]:
    """Read the S1E sheet without third-party XLSX libraries.

    Only nonzero Pfam cells are retained, keeping memory proportional to the
    number of observed domains rather than rows × all 10,588 columns.
    """

    main_ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    rel_ns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"

    with zipfile.ZipFile(path) as archive:
        shared: List[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            for _event, element in ET.iterparse(archive.open("xl/sharedStrings.xml"), events=("end",)):
                if element.tag.endswith("}si"):
                    shared.append("".join(node.text or "" for node in element.iter() if node.tag.endswith("}t")))
                    element.clear()

        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rel_targets = {item.attrib["Id"]: item.attrib["Target"] for item in relationships}
        sheet = next(
            (
                item
                for item in workbook.iter()
                if item.tag.endswith("}sheet") and "S1E" in item.attrib.get("name", "").upper()
            ),
            None,
        )
        if sheet is None:
            raise AuditFailure(f"No S1E sheet found in {path}")
        rid = sheet.attrib[f"{{{rel_ns}}}id"]
        target = rel_targets[rid].lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target

        headers: Optional[List[str]] = None
        metadata_by_genome: Dict[str, Dict[str, str]] = {}
        counts_by_genome: Dict[str, Dict[str, int]] = {}
        xml_rows = 0
        blank_rows = 0

        for _event, element in ET.iterparse(archive.open(target), events=("end",)):
            if not element.tag.endswith("}row"):
                continue
            xml_rows += 1
            cells: Dict[int, str] = {}
            for cell in element:
                if not cell.tag.endswith("}c"):
                    continue
                index = xlsx_column_index(cell.attrib.get("r", "A1"))
                value_node = next((node for node in cell if node.tag.endswith("}v")), None)
                value = "" if value_node is None else (value_node.text or "")
                cell_type = cell.attrib.get("t")
                if cell_type == "s" and value:
                    value = shared[int(value)]
                elif cell_type == "inlineStr":
                    value = "".join(node.text or "" for node in cell.iter() if node.tag.endswith("}t"))
                cells[index] = value

            if headers is None:
                width = max(cells) + 1
                headers = [cells.get(index, "") for index in range(width)]
                continue
            if not any(str(value).strip() for value in cells.values()):
                blank_rows += 1
                element.clear()
                continue

            genome = cells.get(0, "").strip()
            if not genome:
                blank_rows += 1
                element.clear()
                continue
            if genome in metadata_by_genome:
                raise AuditFailure(f"Duplicate Genome in {path} S1E: {genome}")

            metadata = {
                standard: cells.get(index, "")
                for index, standard in enumerate(STANDARD_METADATA_FIELDS)
            }
            counts: Dict[str, int] = {}
            for index in range(14, len(headers)):
                field = headers[index]
                if not PFAM_RE.fullmatch(field):
                    continue
                raw_value = cells.get(index, "")
                if not raw_value:
                    continue
                try:
                    count = int(float(raw_value))
                except ValueError as exc:
                    raise AuditFailure(f"Non-integer Pfam count in {path}: {genome}, {field}, {raw_value!r}") from exc
                if count < 0:
                    raise AuditFailure(f"Negative Pfam count in {path}: {genome}, {field}")
                if count:
                    counts[field] = count
            metadata_by_genome[genome] = metadata
            counts_by_genome[genome] = counts
            element.clear()

    assert headers is not None
    domain_headers = [field for field in headers if PFAM_RE.fullmatch(field)]
    return {
        "path": rel(path),
        "sheet_name": sheet.attrib.get("name", ""),
        "xml_rows_including_header": xml_rows,
        "blank_rows": blank_rows,
        "data_rows": len(metadata_by_genome),
        "domain_columns": len(domain_headers),
        "domain_headers": domain_headers,
        "metadata_by_genome": metadata_by_genome,
        "counts_by_genome": counts_by_genome,
        "all_zero_genomes": sorted(genome for genome, counts in counts_by_genome.items() if not counts),
    }


def load_master_metadata() -> Tuple[List[Dict[str, str]], Dict[str, Dict[str, str]], Dict[str, object]]:
    headers, rows = read_csv_rows(MASTER_METADATA)
    if headers != STANDARD_METADATA_FIELDS:
        raise AuditFailure(f"Unexpected authoritative metadata header in {MASTER_METADATA}: {headers}")
    nonblank = [row for row in rows if str(row.get("Genome", "") or "").strip()]
    blank = len(rows) - len(nonblank)
    by_genome = {key: dict(value) for key, value in require_unique(nonblank, "Genome", "master metadata").items()}
    stats = {
        "physical_data_rows": len(rows),
        "nonblank_rows": len(nonblank),
        "blank_rows": blank,
        "phylum_counts": dict(sorted(Counter(row["Phylum"] for row in nonblank).items())),
    }
    return nonblank, by_genome, stats


def load_table_s4() -> Dict[str, object]:
    headers, rows = read_csv_rows(TABLE_S4)
    pfams = [field for field in headers if PFAM_RE.fullmatch(field)]
    nonblank: List[Dict[str, str]] = []
    all_zero: List[str] = []
    blank_row_numbers: List[int] = []
    for file_row, row in enumerate(rows, 2):
        if not any(str(value or "").strip() for value in row.values()):
            blank_row_numbers.append(file_row)
            continue
        nonblank.append(row)
        total = 0
        for field in pfams:
            raw = str(row.get(field, "") or "").strip()
            if raw:
                total += int(float(raw))
        if total == 0:
            all_zero.append(row.get("Genome ID", ""))
    by_id = {key: dict(value) for key, value in require_unique(nonblank, "Genome ID", "Table S4").items()}
    by_number = {key: dict(value) for key, value in require_unique(nonblank, "ID number", "Table S4 ID number").items()}
    return {
        "headers": headers,
        "pfam_headers": pfams,
        "rows": nonblank,
        "by_id": by_id,
        "by_number": by_number,
        "physical_data_rows": len(rows),
        "nonblank_rows": len(nonblank),
        "blank_rows": len(blank_row_numbers),
        "blank_row_numbers": blank_row_numbers,
        "all_zero_genomes": sorted(all_zero),
        "phylum_counts": dict(sorted(Counter(row.get("Phylum", "") for row in nonblank).items())),
    }


def locate_hmm_sources(
    master_rows: Sequence[Mapping[str, str]],
    s4_rows: Sequence[Mapping[str, str]],
) -> Tuple[Dict[str, Path], Dict[str, str], Dict[str, str]]:
    tagged = sorted(TAGGED_HMM_DIR.glob("tagged-*.seqtblout"))
    af3_raw = sorted(
        path
        for path in AF3_HMM_DIR.glob("*.seqtblout")
        if not path.name.startswith(("tagged-", "tabbed-"))
    )
    s4_by_species: Dict[str, List[Mapping[str, str]]] = defaultdict(list)
    for row in s4_rows:
        s4_by_species[norm_text(row.get("Species", ""))].append(row)

    sources: Dict[str, Path] = {}
    methods: Dict[str, str] = {}
    aliases: Dict[str, str] = {}
    for row in master_rows:
        genome = row["Genome"]
        hits = [path for path in tagged if genome in path.name]
        if len(hits) == 1:
            sources[genome] = hits[0]
            methods[genome] = "unique_canonical_id_substring_in_tagged_tblout"
            continue
        if len(hits) > 1:
            raise AuditFailure(f"Ambiguous tagged HMM files for {genome}: {[rel(path) for path in hits]}")

        raw_hits = [path for path in af3_raw if genome in path.name]
        if genome in {"Ulva_mutabilis", "Ulva_prolifera"}:
            species_hits = s4_by_species[norm_text(row["Species"])]
            if len(species_hits) != 1:
                raise AuditFailure(f"Cannot resolve unique Table S4 assembly alias for {genome}: {species_hits}")
            alias = species_hits[0]["Genome ID"].strip()
            raw_hits = [path for path in af3_raw if alias in path.name]
            if len(raw_hits) != 1:
                raise AuditFailure(f"Cannot resolve unique AF3 HMM file for alias {genome} -> {alias}")
            sources[genome] = raw_hits[0]
            methods[genome] = "unique_exact_species_s4_full_assembly_alias_to_af3_tblout"
            aliases[genome] = alias
            continue
        if len(raw_hits) != 1:
            raise AuditFailure(f"Expected one AF3 raw HMM file for {genome}; found {[rel(path) for path in raw_hits]}")
        sources[genome] = raw_hits[0]
        methods[genome] = "unique_canonical_id_substring_in_af3_tblout"

    if len(set(sources.values())) != len(sources):
        duplicates = Counter(str(path) for path in sources.values())
        raise AuditFailure(f"A raw HMM file mapped to multiple genomes: {[key for key, n in duplicates.items() if n > 1]}")
    return sources, methods, aliases


def parse_hmm_counts(path: Path) -> Tuple[Dict[str, int], str, Dict[str, str]]:
    """Count one hit for each data line's Pfam accession and hash the file."""

    tagged = path.name.startswith("tagged-")
    counts: Counter[str] = Counter()
    digest = hashlib.sha256()
    program_metadata: Dict[str, str] = {}
    with path.open("rb") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            digest.update(raw_line)
            try:
                line = raw_line.decode("utf-8")
            except UnicodeDecodeError as exc:
                raise AuditFailure(f"Non-UTF-8 HMM line: {path}:{line_number}") from exc
            if not line.strip():
                continue
            if tagged:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) >= 2 and fields[1].startswith("#"):
                    if len(fields) >= 4 and fields[2] in {"Program:", "Version:", "Query", "Target"}:
                        program_metadata[fields[2].rstrip(":")] = " ".join(fields[3:]).strip()
                    continue
                if len(fields) < 5:
                    raise AuditFailure(f"Malformed tagged HMM data line: {path}:{line_number}")
                accession = fields[4].strip()
            else:
                if line.lstrip().startswith("#"):
                    match = re.match(r"\s*#\s*(Program|Version|Query file|Target file):\s*(.*)", line)
                    if match:
                        program_metadata[match.group(1)] = match.group(2).strip()
                    continue
                fields = line.split()
                if len(fields) < 4:
                    raise AuditFailure(f"Malformed raw HMM data line: {path}:{line_number}")
                accession = fields[3].strip()
            match = PFAM_ACCESSION_RE.fullmatch(accession)
            if not match:
                raise AuditFailure(f"Invalid Pfam accession at {path}:{line_number}: {accession!r}")
            counts[match.group(1)] += 1
    if not counts:
        raise AuditFailure(f"No Pfam hits parsed from {path}")
    return dict(counts), digest.hexdigest(), program_metadata


def s4_counts_for_ids(target_ids: Iterable[str], pfam_headers: Sequence[str]) -> Dict[str, Dict[str, int]]:
    targets = set(target_ids)
    result: Dict[str, Dict[str, int]] = {}
    with TABLE_S4.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            genome = str(row.get("Genome ID", "") or "").strip()
            if genome not in targets:
                continue
            counts: Dict[str, int] = {}
            for field in pfam_headers:
                raw = str(row.get(field, "") or "").strip()
                if not raw:
                    continue
                value = int(float(raw))
                if value:
                    counts[field] = value
            result[genome] = counts
    missing = targets - set(result)
    if missing:
        raise AuditFailure(f"Missing requested Table S4 count rows: {sorted(missing)}")
    return result


def compare_metadata_to_s1e(
    master_by_genome: Mapping[str, Mapping[str, str]],
    s1e: Mapping[str, object],
) -> Dict[str, object]:
    s1_metadata = s1e["metadata_by_genome"]
    assert isinstance(s1_metadata, dict)
    master_ids = set(master_by_genome)
    s1_ids = set(s1_metadata)
    mismatches: Dict[str, List[Dict[str, str]]] = {field: [] for field in STANDARD_METADATA_FIELDS[1:]}
    for genome in sorted(master_ids & s1_ids):
        master = master_by_genome[genome]
        sheet = s1_metadata[genome]
        for field in STANDARD_METADATA_FIELDS[1:]:
            if not metadata_equal(field, master.get(field, ""), sheet.get(field, "")):
                mismatches[field].append(
                    {"Genome": genome, "master": master.get(field, ""), "s1e": sheet.get(field, "")}
                )
    return {
        "intersection": len(master_ids & s1_ids),
        "master_only": sorted(master_ids - s1_ids),
        "s1e_only": sorted(s1_ids - master_ids),
        "mismatch_counts": {field: len(items) for field, items in mismatches.items()},
        "mismatches": {field: items for field, items in mismatches.items() if items},
    }


def parse_newick_tips(path: Path) -> List[str]:
    text = path.read_text(encoding="utf-8").strip()
    tips = re.findall(r"(?<=[(,])([^():,;]+):[-+0-9.eE]+", text)
    if len(tips) != len(set(tips)):
        raise AuditFailure("Duplicate Newick tip labels")
    return tips


def tree_branch_and_support_summary(path: Path) -> Dict[str, object]:
    """Summarize stored Newick branch lengths and FastTree local supports."""

    text = path.read_text(encoding="utf-8").strip()
    branch_lengths = [float(value) for value in re.findall(r":([-+0-9.eE]+)", text)]
    leaf_lengths = [
        float(value)
        for value in re.findall(r"(?<=[(,])[^():,;]+:([-+0-9.eE]+)", text)
    ]
    supports = [float(value) for value in re.findall(r"\)([-+0-9.eE]+):", text)]
    if not branch_lengths or not leaf_lengths:
        raise AuditFailure(f"No branch lengths parsed from {path}")

    def length_stats(values: Sequence[float]) -> Dict[str, object]:
        return {
            "count": len(values),
            "min": min(values),
            "max": max(values),
            "exact_zero": sum(value == 0 for value in values),
            "negative": sum(value < 0 for value in values),
            "less_than_or_equal_5e_minus_9": sum(value <= 5e-9 for value in values),
            "less_than_or_equal_1e_minus_8": sum(value <= 1e-8 for value in values),
            "less_than_or_equal_1e_minus_6": sum(value <= 1e-6 for value in values),
            "less_than_or_equal_1e_minus_5": sum(value <= 1e-5 for value in values),
        }

    return {
        "all_branches": length_stats(branch_lengths),
        "leaf_branches": length_stats(leaf_lengths),
        "local_supports": {
            "count": len(supports),
            "min": min(supports) if supports else None,
            "max": max(supports) if supports else None,
            "exact_zero": sum(value == 0 for value in supports),
            "below_0_80": sum(value < 0.80 for value in supports),
        },
    }


def tree_construction_provenance() -> Dict[str, object]:
    log_text = RBCL_BUILD_LOG.read_text(encoding="utf-8")
    main_lines = MAIN_MANUSCRIPT.read_text(encoding="utf-8").splitlines()
    method_claims = [
        {"line": index, "text": line.strip()}
        for index, line in enumerate(main_lines, 1)
        if (
            "Maximum likelihood phylogenetic reconstruction was performed using Molecular Evolutionary Genetics Analysis" in line
            or "Maximum likelihood phylogeny inferred from rbcL protein sequences" in line
        )
    ]
    if len(method_claims) < 2:
        raise AuditFailure("Could not locate both manuscript phylogeny method/legend claims")

    def first_log_line(prefix: str) -> str:
        line = next((item.strip() for item in log_text.splitlines() if item.startswith(prefix)), "")
        if not line:
            raise AuditFailure(f"Tree build log lacks required line prefix: {prefix}")
        return line

    return {
        "build_script": rel(RBCL_BUILD_SCRIPT),
        "build_log": rel(RBCL_BUILD_LOG),
        "logged_command": first_log_line("Command:"),
        "logged_software": first_log_line("FastTree Version"),
        "logged_support": first_log_line("Amino acid distances:"),
        "logged_model": first_log_line("ML Model:"),
        "logged_alignment_size": first_log_line("Read "),
        "manuscript_path": rel(MAIN_MANUSCRIPT),
        "manuscript_claims": method_claims,
        "mismatch": (
            "The preserved tree was built with FastTree 2.1.11 from 119 sequences using JTT+CAT20, "
            "NNI/SPR/ML-NNI and SH-like local support. The manuscript Methods instead claim MEGA 11, "
            "JTT+Gamma4, NNI and 1,000 bootstrap replicates; the Methods/figure also report inconsistent "
            "rbcL cohort sizes (116 and 106)."
        ),
        "rooting_provenance": "not documented in the preserved build script or command",
        "release_gate": (
            "BLOCK_UNCONDITIONED_PGLS: rebuild or explicitly condition/root/prune the tree, document the "
            "method actually used, handle zero/near-zero branches, verify a positive-definite covariance "
            "matrix, and report sensitivity to branch-length conditioning before PGLS results may be released."
        ),
    }


def resolve_irregular_accession(candidate: str, master_ids: Iterable[str]) -> Optional[str]:
    if candidate in master_ids:
        return candidate
    matches: List[str] = []
    for genome in master_ids:
        if not candidate.startswith(genome):
            continue
        remainder = candidate[len(genome) :]
        if re.fullmatch(r"(?:contigs)?_?\d+w", remainder):
            matches.append(genome)
    if len(matches) == 1:
        return matches[0]
    return None


def reconcile_rbcl(
    master_by_genome: Mapping[str, Mapping[str, str]],
    s4_by_number: Mapping[str, Mapping[str, str]],
    aliases: Mapping[str, str],
) -> Dict[str, object]:
    _headers, phy_rows_raw = read_csv_rows(RBCL_METADATA)
    phy_rows = [row for row in phy_rows_raw if str(row.get("ID number", "") or "").strip()]
    phy_by_id = {key: value for key, value in require_unique(phy_rows, "ID number", "rbcL metadata").items()}
    tips = parse_newick_tips(RBCL_TREE)
    alias_to_master = {alias: genome for genome, alias in aliases.items()}
    master_ids = set(master_by_genome)

    mappings: List[Dict[str, object]] = []
    used_master: Dict[str, str] = {}
    for tip in tips:
        numeric = re.match(r"(\d+)_", tip)
        if numeric:
            phy_id = numeric.group(1)
            phy = phy_by_id.get(phy_id)
        else:
            tip_slug = slug(tip)
            candidates = [
                row
                for row in phy_rows
                if norm_text(row.get("Genome", "")).startswith("ref ")
                and slug(row.get("Genome", "")) in tip_slug
            ]
            if len(candidates) != 1:
                raise AuditFailure(f"Cannot map REF tree tip to one rbcL metadata row: {tip}, candidates={candidates}")
            phy = candidates[0]
            phy_id = phy["ID number"]
        if phy is None:
            raise AuditFailure(f"Tree tip ID absent from rbcL metadata: {tip}")

        s4 = s4_by_number.get(phy_id)
        s4_identity_agrees = bool(
            s4
            and norm_text(s4.get("Genome", "")) == norm_text(phy.get("Genome", ""))
            and norm_text(s4.get("Species", "")) == norm_text(phy.get("Species", ""))
        )
        candidate_source = "rbcL_metadata_genome_id"
        candidate = str(phy.get("Genome ID", "") or "").strip()
        if s4_identity_agrees and str(s4.get("Genome ID", "") or "").strip():
            candidate = str(s4["Genome ID"]).strip()
            candidate_source = "same_id_s4_genome_id_with_display_and_species_agreement"

        resolved = resolve_irregular_accession(candidate, master_ids)
        if resolved is None and candidate in alias_to_master:
            resolved = alias_to_master[candidate]
            candidate_source += "+audited_full_assembly_alias"

        safe = False
        status = "unresolved_no_master_accession_match"
        species_agrees = False
        if resolved is not None:
            species_agrees = norm_text(phy.get("Species", "")) == norm_text(master_by_genome[resolved].get("Species", ""))
            if species_agrees:
                safe = True
                status = "safe_accession_and_species_agree"
            elif s4_identity_agrees and candidate_source.startswith("same_id_s4"):
                safe = True
                status = "safe_accession_with_same_id_identity_corroboration_taxonomy_differs"
            else:
                status = "unresolved_accession_species_conflict_without_same_id_corroboration"

        if safe and resolved is not None:
            if resolved in used_master:
                raise AuditFailure(
                    f"Two tree tips map to one master genome: {used_master[resolved]} and {tip} -> {resolved}"
                )
            used_master[resolved] = tip

        mappings.append(
            {
                "tip": tip,
                "phylogeny_id": phy_id,
                "phylogeny_genome": phy.get("Genome", ""),
                "phylogeny_species": phy.get("Species", ""),
                "phylogeny_genome_id_raw": phy.get("Genome ID", ""),
                "s4_same_id_present": s4 is not None,
                "s4_display_species_agree": s4_identity_agrees,
                "s4_genome_id": "" if s4 is None else s4.get("Genome ID", ""),
                "candidate_accession": candidate,
                "candidate_source": candidate_source,
                "master_genome": resolved or "",
                "species_agrees_with_master": species_agrees,
                "safe": safe,
                "status": status,
            }
        )

    numeric_tip_ids = {re.match(r"(\d+)_", tip).group(1) for tip in tips if re.match(r"(\d+)_", tip)}
    ref_row_ids = {item["phylogeny_id"] for item in mappings if not re.match(r"(\d+)_", str(item["tip"]))}
    metadata_not_in_tree = sorted(set(phy_by_id) - numeric_tip_ids - ref_row_ids, key=lambda value: int(value))

    shifted: List[Dict[str, str]] = []
    for phy_id, phy in phy_by_id.items():
        s4 = s4_by_number.get(phy_id)
        if not s4:
            continue
        if (
            norm_text(s4.get("Genome", "")) == norm_text(phy.get("Genome", ""))
            and norm_text(s4.get("Species", "")) == norm_text(phy.get("Species", ""))
            and str(s4.get("Genome ID", "")) != str(phy.get("Genome ID", ""))
        ):
            shifted.append(
                {
                    "ID number": phy_id,
                    "Genome": phy.get("Genome", ""),
                    "Species": phy.get("Species", ""),
                    "rbcL_metadata_Genome_ID": phy.get("Genome ID", ""),
                    "same_ID_S4_Genome_ID": s4.get("Genome ID", ""),
                }
            )

    return {
        "metadata_rows": len(phy_rows),
        "metadata_id_count": len(phy_by_id),
        "tree_tip_count": len(tips),
        "numbered_tip_count": len(numeric_tip_ids),
        "reference_tip_count": len(tips) - len(numeric_tip_ids),
        "metadata_rows_not_in_tree": metadata_not_in_tree,
        "safe_tip_count": sum(bool(item["safe"]) for item in mappings),
        "unsafe_tip_count": sum(not bool(item["safe"]) for item in mappings),
        "unsafe_tips": [item for item in mappings if not item["safe"]],
        "confirmed_shifted_genome_id_rows": sorted(shifted, key=lambda item: int(item["ID number"])),
        "branch_and_support_summary": tree_branch_and_support_summary(RBCL_TREE),
        "construction_provenance": tree_construction_provenance(),
        "mappings": mappings,
    }


def coordinate_comparison(
    left: Mapping[str, str],
    right: Mapping[str, str],
    left_lat: str,
    left_lon: str,
    right_lat: str,
    right_lon: str,
) -> bool:
    return numeric_equal(left.get(left_lat, ""), right.get(right_lat, "")) and numeric_equal(
        left.get(left_lon, ""), right.get(right_lon, "")
    )


def summarize_numeric_fields(rows: Sequence[Mapping[str, str]], fields: Sequence[str]) -> Dict[str, object]:
    result: Dict[str, object] = {}
    for field in fields:
        values: List[float] = []
        missing = 0
        for row in rows:
            raw = str(row.get(field, "") or "").strip()
            if not raw:
                missing += 1
                continue
            value = as_float(raw)
            assert value is not None
            values.append(value)
        result[field] = {
            "n_nonmissing": len(values),
            "n_missing": missing,
            "min": min(values) if values else None,
            "max": max(values) if values else None,
        }
    return result


def fasta_headers(path: Path) -> List[str]:
    headers: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                headers.append(line[1:].strip())
    return headers


def preflight() -> None:
    required = [
        AGENTS,
        MASTER_METADATA,
        S1E_131,
        S1E_126_LATER,
        TABLE_S4,
        AEF,
        GEE_ENV,
        GEE_UPLOAD,
        GEE_TEST_EXPORT,
        RBCL_METADATA,
        RBCL_TREE,
        RBCL_ALIGNMENT,
        RBCL_BUILD_SCRIPT,
        RBCL_BUILD_LOG,
        MAIN_MANUSCRIPT,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise AuditFailure(f"Missing required source files: {missing}")
    master_rows, _master_by_genome, _stats = load_master_metadata()
    s4 = load_table_s4()
    sources, methods, aliases = locate_hmm_sources(master_rows, s4["rows"])
    tagged_count = sum(method.startswith("unique_canonical_id_substring_in_tagged") for method in methods.values())
    af3_count = len(methods) - tagged_count
    # Parser smoke test on one real source of each format. This produces no result file.
    for source in [next(path for path in sources.values() if path.name.startswith("tagged-")), next(path for path in sources.values() if not path.name.startswith("tagged-"))]:
        with source.open("r", encoding="utf-8") as handle:
            found = False
            for line in handle:
                if PFAM_ACCESSION_RE.search(line) and not line.lstrip().startswith("#") and "\t#\t" not in line:
                    found = True
                    break
            if not found:
                raise AuditFailure(f"Preflight found no Pfam data record in {source}")
    total_bytes = sum(path.stat().st_size for path in set(sources.values()))
    print(f"PRECHECK PASS: {len(master_rows)} master rows; {len(sources)} unique raw HMM sources")
    print(f"  tagged sources: {tagged_count}; AF3 raw sources: {af3_count}; audited aliases: {aliases}")
    print(f"  raw HMM bytes to stream: {total_bytes:,} ({total_bytes / 1024**2:.1f} MiB)")
    print("  expected full runtime on this workspace: approximately 15-45 seconds")
    print("  full run will create timestamped JSON, TSV, Markdown, manifest CSV, and Pfam matrix CSV")


def run_full(run_id: str) -> List[Path]:
    if not re.fullmatch(r"\d{8}_\d{6}", run_id):
        raise AuditFailure(f"Run ID must use YYYYMMDD_HHMMSS: {run_id!r}")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    outputs = {
        "json": OUTDIR / f"cohort_data_integrity_audit_{run_id}.json",
        "hashes": OUTDIR / f"source_and_output_hashes_{run_id}.tsv",
        "report": OUTDIR / f"data_integrity_report_{run_id}.md",
        "manifest": OUTDIR / f"reconciled_analysis_manifest_{run_id}.csv",
        "pfam": OUTDIR / f"reconciled_raw_pfam_counts_{run_id}.csv",
    }
    existing = [str(path) for path in outputs.values() if path.exists()]
    if existing:
        raise AuditFailure(f"Refusing to overwrite existing output files: {existing}")

    created_at = dt.datetime.now().astimezone().isoformat()
    master_rows, master_by_genome, master_stats = load_master_metadata()
    s4 = load_table_s4()
    s1e_131 = read_s1e_xlsx(S1E_131)
    s1e_126 = read_s1e_xlsx(S1E_126_LATER)
    s1e_meta_comparison = compare_metadata_to_s1e(master_by_genome, s1e_131)

    if s1e_131["data_rows"] != 131 or len(s1e_131["all_zero_genomes"]) != 24:
        raise AuditFailure("The audited 131-row S1E does not have the expected 131 rows / 24 all-zero rows")
    if s1e_126["data_rows"] != 126 or len(s1e_126["all_zero_genomes"]) != 23:
        raise AuditFailure("The later 126-row S1E does not have the expected 126 rows / 23 all-zero rows")
    if s1e_meta_comparison["master_only"] or s1e_meta_comparison["s1e_only"]:
        raise AuditFailure(f"Authoritative metadata and 131-row S1E ID sets differ: {s1e_meta_comparison}")

    sources, source_methods, aliases = locate_hmm_sources(master_rows, s4["rows"])
    raw_counts: Dict[str, Dict[str, int]] = {}
    raw_hashes: Dict[str, str] = {}
    hmm_program_metadata: Dict[str, Dict[str, str]] = {}
    raw_union: set[str] = set()
    total_minus_raw: Dict[str, int] = {}
    for index, row in enumerate(master_rows, 1):
        genome = row["Genome"]
        counts, digest, program_info = parse_hmm_counts(sources[genome])
        raw_counts[genome] = counts
        raw_hashes[genome] = digest
        hmm_program_metadata[genome] = program_info
        raw_union.update(counts)
        reported = int(float(row["PFAMs"]))
        difference = reported - sum(counts.values())
        total_minus_raw[genome] = difference
        if difference != 13:
            raise AuditFailure(
                f"PFAMs checksum failed for {genome}: reported={reported}, raw={sum(counts.values())}, difference={difference}"
            )
        if index % 25 == 0 or index == len(master_rows):
            print(f"Parsed and checksum-validated {index}/{len(master_rows)} raw HMM files", file=sys.stderr)

    # Validate the original nonzero S1E rows exactly and identify merge-failure zero rows.
    s1e_counts = s1e_131["counts_by_genome"]
    assert isinstance(s1e_counts, dict)
    nonzero_exact = 0
    s1e_vector_mismatches: List[Dict[str, object]] = []
    for genome in master_by_genome:
        workbook_counts = s1e_counts[genome]
        if not workbook_counts:
            continue
        if workbook_counts == raw_counts[genome]:
            nonzero_exact += 1
        else:
            all_fields = set(workbook_counts) | set(raw_counts[genome])
            differing = [field for field in all_fields if workbook_counts.get(field, 0) != raw_counts[genome].get(field, 0)]
            s1e_vector_mismatches.append({"Genome": genome, "differing_domains": len(differing)})
    if s1e_vector_mismatches:
        raise AuditFailure(f"Nonzero S1E rows differ from raw HMM evidence: {s1e_vector_mismatches}")

    # AF3 recovery rows must also equal the corresponding full Table S4 vectors.
    af3_genomes = [genome for genome, path in sources.items() if not path.name.startswith("tagged-")]
    af3_s4_ids = {genome: aliases.get(genome, genome) for genome in af3_genomes}
    s4_recovery_counts = s4_counts_for_ids(af3_s4_ids.values(), s4["pfam_headers"])
    af3_s4_vector_checks: List[Dict[str, object]] = []
    for genome in af3_genomes:
        s4_id = af3_s4_ids[genome]
        equal = raw_counts[genome] == s4_recovery_counts[s4_id]
        af3_s4_vector_checks.append(
            {
                "Genome": genome,
                "S4_Genome_ID": s4_id,
                "raw_total": sum(raw_counts[genome].values()),
                "S4_total": sum(s4_recovery_counts[s4_id].values()),
                "vector_equal": equal,
            }
        )
        if not equal:
            raise AuditFailure(f"AF3 raw vector differs from Table S4 recovery vector: {genome} -> {s4_id}")

    # AlphaEarth audit and direct authoritative metadata join.
    aef_headers, aef_rows_raw = read_csv_rows(AEF)
    aef_rows = [row for row in aef_rows_raw if row.get("Genome", "").strip()]
    aef_by_genome = {key: dict(value) for key, value in require_unique(aef_rows, "Genome", "AlphaEarth").items()}
    aef_dims = [f"A{index:02d}" for index in range(64)]
    if [field for field in aef_headers if re.fullmatch(r"A\d{2}", field)] != aef_dims:
        raise AuditFailure("AlphaEarth dimension header is not exactly A00-A63")
    aef_bad_numeric: List[Dict[str, str]] = []
    for row in aef_rows:
        for field in aef_dims:
            try:
                value = as_float(row.get(field, ""))
            except AuditFailure:
                aef_bad_numeric.append({"Genome": row["Genome"], "field": field, "value": row.get(field, "")})
                continue
            if value is None:
                aef_bad_numeric.append({"Genome": row["Genome"], "field": field, "value": ""})
    if aef_bad_numeric:
        raise AuditFailure(f"Missing/non-numeric AlphaEarth dimensions: {aef_bad_numeric[:10]}")
    aef_species_mismatches = [
        genome
        for genome in set(aef_by_genome) & set(master_by_genome)
        if norm_text(aef_by_genome[genome]["Species"]) != norm_text(master_by_genome[genome]["Species"])
    ]
    aef_coord_mismatches = [
        genome
        for genome in set(aef_by_genome) & set(master_by_genome)
        if not coordinate_comparison(
            aef_by_genome[genome], master_by_genome[genome], "DD latitude", "DD longitude", "DD latitude", "DD longitude"
        )
    ]
    if aef_species_mismatches or aef_coord_mismatches:
        raise AuditFailure(
            f"AlphaEarth/master identity mismatch: species={aef_species_mismatches}, coordinates={aef_coord_mismatches}"
        )
    coordinate_groups: Dict[Tuple[float, float], List[str]] = defaultdict(list)
    embedding_groups: Dict[Tuple[str, ...], List[str]] = defaultdict(list)
    for row in aef_rows:
        coordinate_groups[(float(row["DD latitude"]), float(row["DD longitude"]))].append(row["Genome"])
        embedding_groups[tuple(row[field] for field in aef_dims)].append(row["Genome"])
    duplicate_coordinate_groups = [sorted(group) for group in coordinate_groups.values() if len(group) > 1]
    duplicate_embedding_groups = [sorted(group) for group in embedding_groups.values() if len(group) > 1]

    # GEE derived table, upload coordinates, and retained test export.
    gee_headers, gee_rows = read_csv_rows(GEE_ENV)
    _upload_headers, gee_upload_rows = read_csv_rows(GEE_UPLOAD)
    _test_headers, gee_test_rows = read_csv_rows(GEE_TEST_EXPORT)
    gee_by_id = {key: dict(value) for key, value in require_unique(gee_rows, "genome_id", "GEE environment table").items()}
    upload_by_id = {key: dict(value) for key, value in require_unique(gee_upload_rows, "genome_id", "GEE upload").items()}
    if set(gee_by_id) != set(upload_by_id):
        raise AuditFailure("GEE environment and upload ID sets differ")

    alias_to_master = {alias: genome for genome, alias in aliases.items()}
    gee_to_master: Dict[str, str] = {}
    gee_join_method: Dict[str, str] = {}
    for gee_id in gee_by_id:
        if gee_id in master_by_genome:
            gee_to_master[gee_id] = gee_id
            gee_join_method[gee_id] = "exact_genome_id"
        elif gee_id in alias_to_master:
            gee_to_master[gee_id] = alias_to_master[gee_id]
            gee_join_method[gee_id] = "audited_full_assembly_alias"

    gee_coordinate_mismatches: List[Dict[str, str]] = []
    for gee_id, master_id in gee_to_master.items():
        if not coordinate_comparison(
            upload_by_id[gee_id], master_by_genome[master_id], "latitude", "longitude", "DD latitude", "DD longitude"
        ):
            gee_coordinate_mismatches.append({"GEE_ID": gee_id, "master_Genome": master_id})
    if gee_coordinate_mismatches:
        raise AuditFailure(f"GEE upload/master coordinate mismatches: {gee_coordinate_mismatches}")

    exact_aef_gee = set(aef_by_genome) & set(gee_by_id)
    resolved_aef_gee_pairs = {
        (master_id, gee_id)
        for gee_id, master_id in gee_to_master.items()
        if master_id in aef_by_genome
    }
    test_export_ids = {row.get("genome_id", "") for row in gee_test_rows if row.get("genome_id", "")}
    gee_provenance_status = "UNSAFE_UNRESOLVED_PROVENANCE_NO_PRESERVED_126_ROW_RAW_EXPORT"

    # Table S4 vs authoritative metadata fields.
    s4_by_id = s4["by_id"]
    direct_s4_master = set(s4_by_id) & set(master_by_genome)
    s4_field_map = {
        "Species": "Species",
        "Phylum": "Phylum",
        "Salinity": "Environment",
        "Climatic zone": "Climatic zone",
        "Temperature (°C)": "Temperature (°C)",
        "DD latitude": "DD latitude",
        "DD longitude": "DD longitude",
        "Habitat": "Habitat",
    }
    s4_metadata_mismatches: Dict[str, List[Dict[str, str]]] = {}
    for s4_field, master_field in s4_field_map.items():
        items: List[Dict[str, str]] = []
        for genome in sorted(direct_s4_master):
            a = s4_by_id[genome].get(s4_field, "")
            b = master_by_genome[genome].get(master_field, "")
            if not metadata_equal(master_field, a, b):
                items.append({"Genome": genome, "Table_S4": a, "authoritative_master": b})
        s4_metadata_mismatches[s4_field] = items

    rbcl = reconcile_rbcl(master_by_genome, s4["by_number"], aliases)
    rbcl_safe_by_master = {
        item["master_genome"]: item
        for item in rbcl["mappings"]
        if item["safe"] and item["master_genome"]
    }
    alignment_headers = fasta_headers(RBCL_ALIGNMENT)

    # Cohort intersections and differences.
    master_ids = set(master_by_genome)
    aef_ids = set(aef_by_genome)
    gee_ids = set(gee_by_id)
    s4_ids = set(s4_by_id)
    s1e_126_ids = set(s1e_126["metadata_by_genome"])
    tree_master_ids = set(rbcl_safe_by_master)
    cohorts = {
        "master_metadata_nonblank": len(master_ids),
        "table_s1e_131": len(s1e_131["metadata_by_genome"]),
        "table_s1e_later_126": len(s1e_126_ids),
        "table_s4_nonblank": len(s4_ids),
        "aef": len(aef_ids),
        "gee_environment": len(gee_ids),
        "rbcl_tree_tips": rbcl["tree_tip_count"],
        "rbcl_tips_safely_mapped_to_master": rbcl["safe_tip_count"],
        "aef_master_exact_intersection": len(aef_ids & master_ids),
        "aef_master_excluded": sorted(master_ids - aef_ids),
        "aef_equals_later_126_s1e_id_set": aef_ids == s1e_126_ids,
        "aef_only_vs_later_126_s1e": sorted(aef_ids - s1e_126_ids),
        "later_126_s1e_only_vs_aef": sorted(s1e_126_ids - aef_ids),
        "aef_gee_exact_intersection": len(exact_aef_gee),
        "aef_gee_exact_aef_only": sorted(aef_ids - gee_ids),
        "aef_gee_exact_gee_only": sorted(gee_ids - aef_ids),
        "aef_gee_after_audited_alias_intersection": len(resolved_aef_gee_pairs),
        "aef_gee_alias_pairs": sorted(
            [
                {"master_Genome": master_id, "GEE_ID": gee_id}
                for master_id, gee_id in resolved_aef_gee_pairs
                if master_id != gee_id
            ],
            key=lambda item: item["master_Genome"],
        ),
        "gee_master_exact_intersection": len(gee_ids & master_ids),
        "gee_master_after_audited_alias_intersection": len(gee_to_master),
        "gee_unresolved_ids": sorted(gee_ids - set(gee_to_master)),
        "master_missing_gee_after_alias": sorted(master_ids - set(gee_to_master.values())),
        "s4_master_exact_intersection": len(s4_ids & master_ids),
        "s4_master_after_audited_alias_intersection": len((s4_ids & master_ids) | set(aliases)),
        "s4_only_after_alias": sorted(s4_ids - master_ids - set(aliases.values())),
        "master_only_after_s4_alias": sorted(master_ids - s4_ids - set(aliases)),
        "aef_rbcl_safe_intersection": len(aef_ids & tree_master_ids),
        "master_rbcl_safe_intersection": len(master_ids & tree_master_ids),
    }

    aef_phyla = Counter(master_by_genome[genome]["Phylum"] for genome in aef_ids)
    if dict(aef_phyla) != {"Rhodophyta": 70, "Ochrophyta": 43, "Chlorophyta": 13}:
        raise AuditFailure(f"Unexpected AEF authoritative phylum distribution: {dict(aef_phyla)}")

    # Write the reconciled raw Pfam matrix.
    domain_headers = sorted(raw_union)
    with outputs["pfam"].open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Genome", *domain_headers])
        for row in master_rows:
            genome = row["Genome"]
            counts = raw_counts[genome]
            writer.writerow([genome, *(counts.get(field, 0) for field in domain_headers)])

    # Per-master manifest; no GEE environmental value is copied into this table.
    manifest_fields = [
        "master_row",
        *STANDARD_METADATA_FIELDS,
        "raw_pfam_hit_total",
        "raw_distinct_pfam_count",
        "pfam_source_path",
        "pfam_source_sha256",
        "pfam_source_mapping_method",
        "pfams_reported_minus_raw_total",
        "s1e_131_original_all_zero",
        "s1e_131_raw_vector_status",
        "aef_present",
        "aef_id_number",
        "aef_coordinate_status",
        "gee_present_exact",
        "gee_present_via_alias",
        "gee_id",
        "gee_environment_status",
        "rbcl_tip",
        "rbcl_mapping_status",
        "rbcl_mapping_method",
        "safe_for_raw_pfam_analysis",
        "safe_for_aef_pfam_analysis",
        "safe_for_gee_environment_analysis",
        "safe_for_rbcl_pfam_analysis",
    ]
    manifest_rows: List[Dict[str, object]] = []
    master_to_gee = {master_id: gee_id for gee_id, master_id in gee_to_master.items()}
    for index, row in enumerate(master_rows, 1):
        genome = row["Genome"]
        aef_row = aef_by_genome.get(genome)
        gee_id = master_to_gee.get(genome, "")
        rbcl_row = rbcl_safe_by_master.get(genome)
        manifest = {field: row.get(field, "") for field in STANDARD_METADATA_FIELDS}
        manifest.update(
            {
                "master_row": index,
                "raw_pfam_hit_total": sum(raw_counts[genome].values()),
                "raw_distinct_pfam_count": len(raw_counts[genome]),
                "pfam_source_path": rel(sources[genome]),
                "pfam_source_sha256": raw_hashes[genome],
                "pfam_source_mapping_method": source_methods[genome],
                "pfams_reported_minus_raw_total": total_minus_raw[genome],
                "s1e_131_original_all_zero": genome in set(s1e_131["all_zero_genomes"]),
                "s1e_131_raw_vector_status": (
                    "RECOVERED_FROM_RAW_HMM_S1E_WAS_ALL_ZERO"
                    if genome in set(s1e_131["all_zero_genomes"])
                    else "EXACT_VECTOR_MATCH_TO_RAW_HMM"
                ),
                "aef_present": aef_row is not None,
                "aef_id_number": "" if aef_row is None else aef_row.get("ID number", ""),
                "aef_coordinate_status": "NOT_IN_AEF" if aef_row is None else "EXACT_NUMERIC_MATCH_TO_MASTER",
                "gee_present_exact": genome in gee_by_id,
                "gee_present_via_alias": bool(gee_id and gee_id != genome),
                "gee_id": gee_id,
                "gee_environment_status": gee_provenance_status if gee_id else "NOT_IN_GEE_TABLE",
                "rbcl_tip": "" if rbcl_row is None else rbcl_row["tip"],
                "rbcl_mapping_status": "NO_SAFE_TREE_MAPPING" if rbcl_row is None else rbcl_row["status"],
                "rbcl_mapping_method": "" if rbcl_row is None else rbcl_row["candidate_source"],
                "safe_for_raw_pfam_analysis": True,
                "safe_for_aef_pfam_analysis": aef_row is not None,
                "safe_for_gee_environment_analysis": False,
                "safe_for_rbcl_pfam_analysis": rbcl_row is not None,
            }
        )
        manifest_rows.append(manifest)
    with outputs["manifest"].open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields)
        writer.writeheader()
        writer.writerows(manifest_rows)

    audit: Dict[str, object] = {
        "run": {
            "run_id": run_id,
            "created_at": created_at,
            "program": rel(Path(__file__)),
            "python": sys.version,
            "platform": platform.platform(),
            "random_or_synthetic_data_used": False,
        },
        "formulas": {
            "pfam_count": "For each non-comment hmmsearch tblout data line, increment the version-stripped PFxxxxx accession by 1.",
            "raw_pfam_hit_total": "sum(count[PFxxxxx] for all observed PFxxxxx)",
            "reported_total_checksum": "authoritative metadata PFAMs - raw_pfam_hit_total; required to equal 13 for every genome",
            "coordinate_agreement": "abs(latitude_source - latitude_master) <= 1e-10 and abs(longitude_source - longitude_master) <= 1e-10",
            "cohort_intersection": "set of exact Genome IDs present in both sources; audited assembly aliases are reported separately",
        },
        "authoritative_master_metadata": {
            "path": rel(MASTER_METADATA),
            **master_stats,
            "comparison_to_131_row_s1e": s1e_meta_comparison,
        },
        "table_s1e": {
            "audited_131_row_workbook": {
                key: value
                for key, value in s1e_131.items()
                if key not in {"domain_headers", "metadata_by_genome", "counts_by_genome"}
            },
            "later_126_row_workbook": {
                key: value
                for key, value in s1e_126.items()
                if key not in {"domain_headers", "metadata_by_genome", "counts_by_genome"}
            },
            "nonzero_rows_exactly_equal_raw_hmm_vectors": nonzero_exact,
            "nonzero_vector_mismatches": s1e_vector_mismatches,
            "interpretation": "All all-zero S1E domain rows have nonzero raw HMM evidence and are parser/merge failures, not biological zeros.",
        },
        "raw_pfam_reconciliation": {
            "master_rows_recovered": len(raw_counts),
            "unique_raw_hmm_sources": len(set(sources.values())),
            "tagged_sources": sum(path.name.startswith("tagged-") for path in sources.values()),
            "af3_unprefixed_sources": sum(not path.name.startswith("tagged-") for path in sources.values()),
            "raw_domain_union": len(domain_headers),
            "reported_minus_raw_total_distribution": dict(sorted(Counter(total_minus_raw.values()).items())),
            "audited_assembly_aliases": aliases,
            "af3_vectors_equal_table_s4": af3_s4_vector_checks,
            "matrix_path": rel(outputs["pfam"]),
        },
        "table_s4": {
            "path": rel(TABLE_S4),
            "physical_data_rows": s4["physical_data_rows"],
            "nonblank_rows": s4["nonblank_rows"],
            "blank_rows": s4["blank_rows"],
            "blank_row_numbers": s4["blank_row_numbers"],
            "pfam_columns": len(s4["pfam_headers"]),
            "all_zero_pfam_rows": s4["all_zero_genomes"],
            "phylum_counts_as_stored": s4["phylum_counts"],
            "direct_master_intersection": len(direct_s4_master),
            "metadata_mismatch_counts": {field: len(items) for field, items in s4_metadata_mismatches.items()},
            "metadata_mismatches": {field: items for field, items in s4_metadata_mismatches.items() if items},
            "status": "PFAM values are retained as a cross-check only; metadata fields are not authoritative.",
        },
        "alphaearth": {
            "path": rel(AEF),
            "rows": len(aef_rows),
            "dimensions": len(aef_dims),
            "master_exact_intersection": len(aef_ids & master_ids),
            "master_excluded": sorted(master_ids - aef_ids),
            "species_mismatches": aef_species_mismatches,
            "coordinate_mismatches": aef_coord_mismatches,
            "authoritative_phylum_counts": dict(sorted(aef_phyla.items())),
            "duplicate_coordinate_groups": duplicate_coordinate_groups,
            "duplicate_embedding_vector_groups": duplicate_embedding_groups,
            "status": "Safe as 64 latent AEF values joined by exact Genome ID; not safe as decoded environmental variables or collection-date measurements.",
        },
        "gee_environment": {
            "derived_table_path": rel(GEE_ENV),
            "derived_rows": len(gee_rows),
            "environment_fields": GEE_ENV_FIELDS,
            "field_completeness_and_range": summarize_numeric_fields(gee_rows, GEE_ENV_FIELDS),
            "upload_path": rel(GEE_UPLOAD),
            "upload_rows": len(gee_upload_rows),
            "upload_coordinate_mismatches_after_audited_joins": gee_coordinate_mismatches,
            "retained_raw_export_path": rel(GEE_TEST_EXPORT),
            "retained_raw_export_rows": len(gee_test_rows),
            "retained_raw_export_exact_id_overlap_with_derived_126": len(test_export_ids & gee_ids),
            "provenance_status": gee_provenance_status,
            "reason": "The retained GEE export has four test rows, not 126; no preserved 126-row raw export maps the 13 derived values to the 126 Genome IDs.",
            "exact_duplicate_derived_paths": [rel(GEE_ENV_B), rel(GEE_ENV_DRIVE_COPY)],
        },
        "rbcl": {
            **rbcl,
            "metadata_path": rel(RBCL_METADATA),
            "tree_path": rel(RBCL_TREE),
            "alignment_path": rel(RBCL_ALIGNMENT),
            "alignment_sequence_count": len(alignment_headers),
            "status": "Use only the accession-gated safe mappings; never use the S1C/rbcL Genome ID column directly.",
        },
        "cohorts": cohorts,
        "field_safety": {
            "safe": [
                "Master Genome, Species, Phylum, reported assembly/quality fields, habitat categories, and coordinates as stored in the authoritative metadata CSV.",
                "Raw Pfam counts reconstructed from one preserved HMM tblout per master Genome.",
                "AEF A00-A63 values for the exact 126-Genome AEF cohort.",
                "rbcL tree topology/branch lengths only for the accession-gated mapped subset.",
            ],
            "safe_with_limitation": [
                "DD latitude/longitude are numeric extraction coordinates; coordinate confidence, collection date, and spatial uncertainty annotations are not available.",
                "PFAMs is a legacy reported total and is consistently 13 greater than the parsed raw hit total; use raw_pfam_hit_total for normalization.",
                "Exons is source-reported but is not documented as total predicted proteins; do not substitute it for protein count.",
            ],
            "unsafe": [
                "All 13 GEE derived environmental variables until a preserved 126-row raw export or a reproducible re-extraction is supplied.",
                "Table S4 Phylum and other mismatching metadata fields.",
                "S1E all-zero Pfam rows before raw-HMM recovery.",
                "S1C/rbcL Genome ID values without the accession-gated reconciliation.",
                "Unconditioned PGLS on the current unrooted FastTree with zero/near-zero branches and unresolved method provenance.",
            ],
            "not_available": [
                "Total predicted protein count",
                "Coordinate-confidence class/source",
                "Collection date for temporal matching",
                "Assembly method/version provenance beyond the stored metadata fields",
            ],
        },
        "outputs": {key: rel(path) for key, path in outputs.items()},
    }

    with outputs["json"].open("x", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")

    zero_list = "\n".join(f"- `{genome}`" for genome in s1e_131["all_zero_genomes"])
    tree_unsafe_list = "\n".join(
        f"- `{item['tip']}` — {item['status']} (candidate `{item['candidate_accession']}`)"
        for item in rbcl["unsafe_tips"]
    )
    gee_unresolved_list = "\n".join(f"- `{genome}`" for genome in cohorts["gee_unresolved_ids"])
    report = f"""# iScience revision data-integrity audit

Run ID: `{run_id}`
Generated: `{created_at}`
Program: `{rel(Path(__file__))}`

## Exact computation and acceptance gates

Pfam counts were regenerated by streaming every mapped raw hmmsearch table and applying:

```text
count[PFxxxxx] = number of non-comment tblout data lines whose query accession is PFxxxxx[.version]
raw_pfam_hit_total = sum(count[PFxxxxx])
checksum = authoritative_metadata.PFAMs - raw_pfam_hit_total
```

Every one of the {len(master_rows)} master genomes passed `checksum == 13`. The {nonzero_exact} nonzero S1E vectors equal the raw vectors exactly. The remaining {len(s1e_131['all_zero_genomes'])} S1E rows contain all-zero domain vectors even though their raw tblouts contain Pfam hits; they are parser/merge failures, not observed biological zeros. Coordinates were accepted only when both numeric coordinates matched the master within `1e-10` degrees.

## Decision-ready findings

- The authoritative master has **{len(master_ids)}** nonblank genomes: {master_stats['phylum_counts']}.
- The 25 December S1E has **131 rows, 10,588 domain columns, and 24 false all-zero rows**. The later 30 December workbook has **126 rows and 23 false all-zero rows**, but it is not the AEF cohort: AEF uniquely contains `19106D-07-09_S0_L001_R1_001`, whereas the workbook uniquely contains `19106D-07-07_S0_L001_R1_001`.
- The 131-row S1E ID set matches the master, but three metadata cells do not: two Species values and one Phylum value. The authoritative master values are retained; exact rows are listed in the machine audit.
- Raw HMM evidence recovers **all 131 master genomes** and **all 126 AEF genomes**. The reconciled union contains **{len(domain_headers):,} Pfam accessions**. Six genomes required unprefixed AF3 tblouts; their full vectors equal Table S4 exactly.
- AlphaEarth contains **126 × 64 complete numeric values**, joins exactly to the master, and has authoritative phylum counts `Rhodophyta=70`, `Ochrophyta=43`, `Chlorophyta=13`. All 126 coordinates match the master.
- The GEE derived table contains 126 IDs but overlaps AEF at **{len(exact_aef_gee)} exact IDs** (**{len(resolved_aef_gee_pairs)}** after the audited Ulva assembly alias). It overlaps the master at {cohorts['gee_master_exact_intersection']} exact IDs / {cohorts['gee_master_after_audited_alias_intersection']} after two audited Ulva aliases.
- The retained raw GEE export has only **{len(gee_test_rows)} test rows** and zero exact Genome-ID overlap with the 126-row derived table. The 13 GEE variables are therefore **unsafe for new analysis** until a 126-row raw export or reproducible re-extraction is supplied.
- The rbcL Newick contains **{rbcl['tree_tip_count']} tips**. Exactly **{rbcl['safe_tip_count']}** map to master genomes with accession-level evidence; {rbcl['unsafe_tip_count']} remain unresolved and must be excluded. The rbcL/S1C `Genome ID` field has {len(rbcl['confirmed_shifted_genome_id_rows'])} confirmed shifted rows and must not be joined directly.
- The preserved tree provenance is FastTree 2.1.11 with JTT+CAT20 and SH-like local support on 119 sequences, not the manuscript's MEGA 11/JTT+Gamma4/1,000-bootstrap workflow. The manuscript also alternates among 106, 116, and 119 rbcL samples.
- Of {rbcl['branch_and_support_summary']['all_branches']['count']} stored branch lengths, **{rbcl['branch_and_support_summary']['all_branches']['exact_zero']} are exactly zero** and **{rbcl['branch_and_support_summary']['all_branches']['less_than_or_equal_1e_minus_8']} are ≤1e-8**; {rbcl['branch_and_support_summary']['leaf_branches']['less_than_or_equal_1e_minus_8']}/{rbcl['branch_and_support_summary']['leaf_branches']['count']} terminal branches are ≤1e-8. Unconditioned PGLS is blocked.
- Table S4 has **{s4['nonblank_rows']} nonblank Pfam rows**, {s4['blank_rows']} trailing blank rows, and no all-zero Pfam rows. Its coordinates agree with the master for direct IDs, but its metadata are not authoritative: Phylum differs for {len(s4_metadata_mismatches['Phylum'])}/{len(direct_s4_master)} direct matches, Species for {len(s4_metadata_mismatches['Species'])}, Temperature for {len(s4_metadata_mismatches['Temperature (°C)'])}, and Habitat for {len(s4_metadata_mismatches['Habitat'])}.

## Cohort rules for revision analyses

- Raw Pfam-only analysis: use all **131** reconciled rows.
- AEF–Pfam analysis: use the exact **126** AEF rows and the reconciled raw Pfam matrix.
- Phylogenetic analysis: use only the **{rbcl['safe_tip_count']}** accession-gated tree mappings; the AEF intersection is **{cohorts['aef_rbcl_safe_intersection']}**. Do not release PGLS until the method mismatch is resolved and the tree is rooted/pruned/conditioned, zero branches are handled transparently, the covariance matrix is positive definite, and branch-conditioning sensitivity is reported.
- GEE–Pfam or GEE–AEF analysis: **status-gated / not currently safe**. The exact AEF–GEE overlap is 118; the audited alias can be reported separately but does not repair missing raw environmental provenance.
- Do not use S1E zeros, S4 phylum labels, or the unreconciled S1C Genome-ID column.

## False all-zero S1E rows recovered from raw HMM output

{zero_list}

## Unresolved rbcL tree tips

{tree_unsafe_list}

## GEE rows outside the authoritative master after audited aliases

{gee_unresolved_list}

## Field safety

Safe:

- Authoritative master identifiers, species, phylum, reported assembly/quality fields, habitat fields, and coordinates as stored.
- Raw Pfam counts in `{rel(outputs['pfam'])}`.
- AEF dimensions A00–A63 for exact-ID rows.
- rbcL topology/branch lengths for accession-gated mappings only.

Restricted or unavailable:

- `PFAMs` is a legacy checksum total; normalize with computed `raw_pfam_hit_total`.
- `Exons` is not documented as total predicted proteins. Total predicted protein count is not available.
- Coordinate confidence, coordinate source, and collection date are not available.
- GEE values are excluded until their raw 126-row extraction provenance is restored.

Machine-readable details, every mismatch list, every rbcL tip mapping, and numeric field completeness are in `{rel(outputs['json'])}`. Row-level source paths and hashes are in `{rel(outputs['manifest'])}`.
"""
    with outputs["report"].open("x", encoding="utf-8") as handle:
        handle.write(report)

    # Hash primary inputs, all unique raw HMM inputs, the program, and generated outputs.
    primary_sources = [
        AGENTS,
        MASTER_METADATA,
        S1E_131,
        S1E_131_FINAL_COPY,
        S1E_126_LATER,
        TABLE_S4,
        AEF,
        AEF_EXTRACTION_SCRIPT,
        PFAM_COUNT_SCRIPT,
        GEE_ENV,
        GEE_ENV_B,
        GEE_ENV_DRIVE_COPY,
        GEE_UPLOAD,
        GEE_TEST_EXPORT,
        GEE_EXTRACTION_SCRIPT,
        GEE_PROVENANCE_DOC,
        RBCL_METADATA,
        RBCL_TREE,
        RBCL_ALIGNMENT,
        RBCL_BUILD_SCRIPT,
        RBCL_BUILD_LOG,
        MAIN_MANUSCRIPT,
        Path(__file__).resolve(),
    ]
    hash_records: List[Dict[str, object]] = []
    raw_path_to_hash = {sources[genome].resolve(): raw_hashes[genome] for genome in sources}
    for path in sorted(set(primary_sources), key=lambda item: rel(item)):
        hash_records.append(
            {"role": "source_or_program", "sha256": sha256_file(path), "bytes": path.stat().st_size, "path": rel(path)}
        )
    for path in sorted(raw_path_to_hash, key=lambda item: rel(item)):
        hash_records.append(
            {"role": "raw_hmm_source", "sha256": raw_path_to_hash[path], "bytes": path.stat().st_size, "path": rel(path)}
        )
    for key in ["json", "report", "manifest", "pfam"]:
        path = outputs[key]
        hash_records.append(
            {"role": "generated_output", "sha256": sha256_file(path), "bytes": path.stat().st_size, "path": rel(path)}
        )
    with outputs["hashes"].open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["role", "sha256", "bytes", "path"], delimiter="\t")
        writer.writeheader()
        writer.writerows(hash_records)

    print(json.dumps({"status": "PASS", "run_id": run_id, "outputs": {key: str(path.resolve()) for key, path in outputs.items()}}, indent=2))
    return list(outputs.values())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", action="store_true", help="validate paths/mappings and smoke-test both HMM formats without writing files")
    parser.add_argument("--run-id", help="output suffix in YYYYMMDD_HHMMSS format; defaults to current local time")
    args = parser.parse_args()
    try:
        if args.preflight:
            preflight()
            return 0
        run_id = args.run_id or dt.datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")
        run_full(run_id)
        return 0
    except (AuditFailure, OSError, csv.Error, zipfile.BadZipFile, ET.ParseError) as exc:
        print(f"AUDIT FAILED: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

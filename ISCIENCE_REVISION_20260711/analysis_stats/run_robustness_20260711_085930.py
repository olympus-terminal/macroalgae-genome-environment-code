#!/usr/bin/env python3
"""Reproducible robustness analyses for the iScience revision.

Scientific inputs are existing project data only. Random-number generation is
used solely for paired bootstrap resampling and structured permutation of real
observations. No synthetic observations or hard-coded result values are used.

Primary analyses
----------------
1. Reconstruct the canonical 126-genome Pfam matrix from raw HMM-search tblout
   files, with a per-genome provenance/hash manifest.
2. Status-gate headline GEE reanalysis until a complete, provenance-verified
   environmental export is supplied; the incomplete retained exports are not
   used for inference.
3. Add assembly-size, BUSCO, peptide-search-space, and phylum covariates;
   BUSCO >=50% and >=70% sensitivity subsets; coordinate-site inference;
   structured site permutations within phylum-composition strata; and
   explicit-tree eigenvector models.
4. Re-evaluate PF00092 Arabian Gulf enrichment across the same representations.
5. Run a clearly labelled archived-latent-feature sensitivity for priority AEF
   pairs without assigning environmental meaning to latent axes.

The script writes timestamped outputs and never overwrites a prior run.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import platform
import re
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import rankdata
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from Bio import Phylo
import Bio


SCRIPT_VERSION = "2026-07-18.1-public-release"
PFAM_RE = re.compile(r"\b(PF\d{5})(?:\.\d+)?\b")
VALID_PHYLA = ("Rhodophyta", "Ochrophyta", "Chlorophyta")
DEFAULT_SEED = 20260711
DEFAULT_PERMUTATIONS = 9999
DEFAULT_BOOTSTRAPS = 999
PRIMARY_BUSCO_THRESHOLD = 50.0
STRICT_BUSCO_THRESHOLD = 70.0
GEO_BLOCK_DEGREES = 20.0
EXPECTED_MASTER_MINUS_RAW_TOTAL = 13.0
GEE_ENV_COLUMNS = [
    "sst_mean_c", "sst_min_c", "sst_max_c", "sst_annual_range_c",
    "sst_summer_c", "sst_winter_c", "chlorophyll_mean_mg_m3",
    "chlorophyll_max_mg_m3", "chlorophyll_std_mg_m3", "poc_mean_mg_m3",
    "depth_meters", "distance_coast_km", "water_clarity_ratio",
]
ANALYSIS_PFAMS = ("PF00092", "PF01638", "PF10988", "PF13411")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validate-only", action="store_true",
                        help="Validate and parse every input but skip statistical outputs.")
    parser.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    parser.add_argument("--bootstraps", type=int, default=DEFAULT_BOOTSTRAPS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--tree-map", type=Path, default=None,
                        help="Optional integrity-reviewed CSV with columns tip, Genome, status.")
    parser.add_argument("--verified-gee-csv", type=Path, default=None,
                        help="Optional provenance-restored AEF-126 baseline GEE export; archived GEE is never used.")
    parser.add_argument("--verified-gee-provenance", type=Path, default=None,
                        help="JSON provenance record containing csv_sha256 and status=VERIFIED_REEXTRACTION.")
    return parser.parse_args()


def now_stamp() -> str:
    return datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def relative(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def normalize_text(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def safe_numeric(series: pd.Series, name: str) -> pd.Series:
    result = pd.to_numeric(series, errors="coerce")
    if result.isna().any():
        bad = series[result.isna()].head(5).tolist()
        raise ValueError(f"{name} has nonnumeric/missing values; examples: {bad}")
    return result


def git_metadata(root: Path) -> dict:
    def run(*args: str) -> str:
        try:
            return subprocess.run(
                args, cwd=root, check=True, capture_output=True, text=True
            ).stdout.strip()
        except Exception as exc:  # provenance aid; analysis must not fail if not a git repo
            return f"unavailable: {type(exc).__name__}"

    return {
        "commit": run("git", "rev-parse", "HEAD"),
        "status_porcelain": run("git", "status", "--porcelain"),
    }


def input_record(path: Path, root: Path, with_hash: bool = True) -> dict:
    stat = path.stat()
    return {
        "path": relative(path, root),
        "realpath": str(path.resolve()),
        "bytes": stat.st_size,
        "mtime": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
        "sha256": sha256_file(path) if with_hash else None,
    }


def locate_raw_tblout(
    genome: str,
    species: str,
    primary_dir: Path,
    fallback_dir: Path,
    legacy_table: pd.DataFrame,
) -> tuple[Path, str, str | None]:
    """Return a unique raw tblout and a documented mapping rule.

    Priority is a unique filename containing the canonical genome ID. If no
    primary file exists, an unprefixed fallback raw file containing the exact
    canonical ID is used. The final fallback is a unique exact-species alias
    from the legacy table, used only when both the metadata row and raw filename
    are unique. This rule currently resolves the short Ulva identifier and is
    recorded in the manifest rather than silently renamed.
    """
    primary = [p for p in primary_dir.glob("*.seqtblout") if genome in p.name]
    if len(primary) == 1:
        return primary[0], "exact_canonical_id_primary", None
    if len(primary) > 1:
        raise ValueError(f"Ambiguous primary tblout for {genome}: {primary}")

    fallback = [
        p for p in fallback_dir.glob("*.seqtblout")
        if genome in p.name and not p.name.startswith(("tagged-", "tabbed-"))
    ]
    if len(fallback) == 1:
        return fallback[0], "exact_canonical_id_fallback", None
    if len(fallback) > 1:
        raise ValueError(f"Ambiguous fallback tblout for {genome}: {fallback}")

    exact_species = legacy_table[
        legacy_table["Species"].map(normalize_text) == normalize_text(species)
    ].copy()
    aliases = exact_species["Genome ID"].dropna().astype(str).unique().tolist()
    alias_hits: list[tuple[Path, str]] = []
    for alias in aliases:
        hits = [
            p for p in fallback_dir.glob("*.seqtblout")
            if alias in p.name and not p.name.startswith(("tagged-", "tabbed-"))
        ]
        alias_hits.extend((p, alias) for p in hits)
    unique_hits = {(str(p.resolve()), alias) for p, alias in alias_hits}
    if len(unique_hits) == 1:
        hit, alias = next(iter(unique_hits))
        return Path(hit), "exact_species_accession_alias_fallback", alias
    raise FileNotFoundError(
        f"No unique traceable raw tblout for {genome} ({species}); "
        f"primary={primary}, fallback={fallback}, alias_hits={alias_hits}"
    )


def parse_tblout(path: Path) -> tuple[Counter, str, int]:
    """Stream a HMM-search tblout, count Pfam accessions, and hash exact bytes."""
    counts: Counter = Counter()
    digest = hashlib.sha256()
    parsed_lines = 0
    with path.open("rb") as handle:
        for raw_line in handle:
            digest.update(raw_line)
            line = raw_line.decode("utf-8", errors="replace")
            if line.lstrip().startswith("#") or "\t#\t" in line:
                continue
            match = PFAM_RE.search(line)
            if match:
                counts[match.group(1)] += 1
                parsed_lines += 1
    if parsed_lines == 0:
        raise ValueError(f"No Pfam records parsed from nonempty tblout {path}")
    return counts, digest.hexdigest(), parsed_lines


def locate_final_peptide_fasta(genome: str, fasta_dir: Path) -> Path | None:
    pattern = "*_5x5L_BLEACHd.aa.fa_BLEACH_6.aa.fa"
    hits = [p for p in fasta_dir.glob(pattern) if genome in p.name]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        raise ValueError(f"Ambiguous final peptide FASTA for {genome}: {hits}")
    return None


def count_fasta_records(path: Path) -> tuple[int, str]:
    """Count exact FASTA records and hash the file in one streaming pass."""
    count = 0
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for line in handle:
            digest.update(line)
            if line.startswith(b">"):
                count += 1
    if count <= 0:
        raise ValueError(f"No FASTA records found in {path}")
    return count, digest.hexdigest()


def load_inputs(root: Path) -> dict[str, Path]:
    inputs = {
        "aef": root / "AlphaEarth/CSV/alphaearth_embeddings_20251019_122918.csv",
        "master": root / "AlphaEarth/CSV/Metadata_Table_macroalgae-published.csv",
        "final_s1": root / "TABLES/FINAL/Table_S1_25DEC_2025-main.xlsx",
        "legacy_counts": root / "Table_S4_meta-and_pfams.csv",
        "tree": root / "DATA/rbcL_phylogenetic_tree_20251118_173538.nwk",
        "phy_metadata": root / "MACROALGAE_PHYLOGENIES.csv",
        "raw_primary": root / "AlphaEarth/TAGGED_HMMsearch-raw-out",
        "raw_fallback": root / "AF3/transfer_hmmsearch_tblout",
        "peptide_fastas": root / "AF3/base_seqs",
        "integrity_manifest": root / "ISCIENCE_REVISION_20260711/integrity/reconciled_analysis_manifest_20260711_110650.csv",
        "integrity_report": root / "ISCIENCE_REVISION_20260711/integrity/data_integrity_report_20260711_110650.md",
        "priority_annotations": root / "ISCIENCE_REVISION_20260711/annotations/priority_interpro_annotations_20260711_114039.csv",
    }
    missing = [str(p) for p in inputs.values() if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Required inputs missing: {missing}")
    return inputs


def load_canonical_cohort(inputs: Mapping[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    aef = pd.read_csv(inputs["aef"])
    master = pd.read_csv(inputs["master"], encoding="utf-8-sig")
    master = master[master["Genome"].notna()].copy()
    legacy = pd.read_csv(inputs["legacy_counts"], low_memory=False)

    if aef["Genome"].duplicated().any() or master["Genome"].duplicated().any():
        raise ValueError("Canonical AEF or master Genome IDs are not unique")
    cohort = aef[["Genome", "Species", "DD latitude", "DD longitude"]].merge(
        master,
        on="Genome",
        how="left",
        suffixes=("_aef", ""),
        validate="one_to_one",
    )
    if len(cohort) != len(aef) or cohort["Phylum"].isna().any():
        raise ValueError("Master metadata does not completely cover canonical AEF cohort")
    observed_phyla = set(cohort["Phylum"])
    if observed_phyla != set(VALID_PHYLA):
        raise ValueError(f"Unresolved phylum labels in authoritative master: {observed_phyla}")

    for column in ["Nucleotides", "Exons", "PFAMs", "BUSCOs-%present", "DD latitude", "DD longitude"]:
        cohort[column] = safe_numeric(cohort[column], column)
    if (cohort[["Nucleotides", "Exons", "PFAMs"]] <= 0).any().any():
        raise ValueError("Nonpositive assembly/search-space values in master metadata")

    cohort["Species"] = cohort["Species"].astype(str)
    cohort["canonical_order"] = np.arange(len(cohort))
    return cohort, aef, legacy


def validate_verified_gee_ingest(
    csv_path: Path | None,
    provenance_path: Path | None,
    canonical_genomes: set[str],
) -> tuple[pd.DataFrame | None, dict]:
    """Gate a future Earth Engine re-extraction without consulting archived GEE."""
    if csv_path is None and provenance_path is None:
        return None, {"status": "NOT_SUPPLIED_GEE_REEXTRACTION_REMAINS_GATED"}
    if csv_path is None or provenance_path is None:
        raise ValueError("--verified-gee-csv and --verified-gee-provenance must be supplied together")
    csv_path = csv_path.resolve()
    provenance_path = provenance_path.resolve()
    if not csv_path.exists() or not provenance_path.exists():
        raise FileNotFoundError(f"Verified GEE inputs missing: {csv_path}, {provenance_path}")
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    if str(provenance.get("status", "")).upper() != "VERIFIED_REEXTRACTION":
        raise ValueError("GEE provenance JSON must declare status=VERIFIED_REEXTRACTION")
    digest = sha256_file(csv_path)
    if provenance.get("csv_sha256") != digest:
        raise ValueError("GEE CSV SHA256 does not match provenance JSON")
    data = pd.read_csv(csv_path)
    id_column = "Genome" if "Genome" in data.columns else "genome_id" if "genome_id" in data.columns else None
    if id_column is None:
        raise ValueError("Verified GEE CSV requires Genome or genome_id")
    # Re-extraction may contain sensitivity scenarios; baseline is explicitly
    # selected rather than mixing buffers/jitters with point estimates.
    if "scenario" in data.columns:
        baseline = data[data["scenario"].astype(str).str.lower().eq("point")].copy()
        if "replicate" in baseline.columns:
            baseline = baseline[pd.to_numeric(baseline["replicate"], errors="coerce").fillna(0).eq(0)]
    else:
        baseline = data.copy()
    missing_columns = [c for c in GEE_ENV_COLUMNS if c not in baseline.columns]
    if missing_columns:
        raise ValueError(f"Verified GEE baseline lacks required variables: {missing_columns}")
    if baseline[id_column].duplicated().any():
        raise ValueError("Verified GEE baseline has duplicate genome IDs")
    ids = set(baseline[id_column].astype(str))
    if ids != canonical_genomes:
        raise ValueError(
            f"Verified GEE baseline must exactly equal canonical AEF-126 IDs; "
            f"missing={sorted(canonical_genomes-ids)[:10]}, extra={sorted(ids-canonical_genomes)[:10]}"
        )
    numeric = baseline[GEE_ENV_COLUMNS].apply(pd.to_numeric, errors="coerce")
    audit = {
        "status": "VERIFIED_REEXTRACTION_INGESTED_READY_FOR_SITE_AWARE_GEE_ANALYSIS",
        "csv_realpath": str(csv_path),
        "provenance_realpath": str(provenance_path),
        "csv_sha256": digest,
        "baseline_rows": len(baseline),
        "complete_counts": numeric.notna().sum().to_dict(),
        "follow_up": "run site-aware GEE association module; do not use archived derived GEE tables",
    }
    return baseline, audit


def reconstruct_raw_pfam_matrix(
    cohort: pd.DataFrame,
    inputs: Mapping[str, Path],
    root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    counters: dict[str, Counter] = {}
    manifest_rows: list[dict] = []
    legacy = pd.read_csv(inputs["legacy_counts"], low_memory=False)

    print(f"Parsing raw HMM-search outputs for {len(cohort)} canonical genomes...", flush=True)
    for position, row in cohort.iterrows():
        genome = str(row["Genome"])
        path, rule, alias = locate_raw_tblout(
            genome,
            str(row["Species"]),
            inputs["raw_primary"],
            inputs["raw_fallback"],
            legacy,
        )
        counts, digest, parsed = parse_tblout(path)
        counters[genome] = counts
        manifest_rows.append({
            "Genome": genome,
            "Species": row["Species"],
            "raw_tblout": relative(path, root),
            "mapping_rule": rule,
            "mapping_alias": alias,
            "sha256": digest,
            "bytes": path.stat().st_size,
            "parsed_pfam_hits": parsed,
            "unique_pfams": len(counts),
        })
        if (position + 1) % 20 == 0 or position + 1 == len(cohort):
            print(f"  parsed {position + 1}/{len(cohort)}", flush=True)

    all_pfams = sorted({pfam for counts in counters.values() for pfam in counts})
    matrix = np.zeros((len(cohort), len(all_pfams)), dtype=np.int32)
    col_index = {pfam: j for j, pfam in enumerate(all_pfams)}
    for i, genome in enumerate(cohort["Genome"]):
        for pfam, count in counters[str(genome)].items():
            matrix[i, col_index[pfam]] = count
    pfam_df = pd.DataFrame(matrix, index=cohort["Genome"], columns=all_pfams)
    pfam_df.index.name = "Genome"
    manifest = pd.DataFrame(manifest_rows)

    raw_totals = pfam_df.sum(axis=1).astype(float)
    master_totals = cohort.set_index("Genome")["PFAMs"].astype(float)
    differences = master_totals - raw_totals
    if not differences.eq(EXPECTED_MASTER_MINUS_RAW_TOTAL).all():
        mismatch = pd.DataFrame({"raw": raw_totals, "master": master_totals, "difference": differences})
        raise ValueError(
            f"Raw/master checksum offset must equal {EXPECTED_MASTER_MINUS_RAW_TOTAL:g} "
            f"for every genome:\n{mismatch.head(20)}"
        )
    if (raw_totals <= 0).any():
        raise ValueError("At least one canonical genome has no raw Pfam hits")
    manifest["master_PFAMs_field"] = manifest["Genome"].map(master_totals)
    manifest["master_minus_raw_total"] = manifest["Genome"].map(differences)
    return pfam_df, manifest


def derive_peptide_denominators(
    cohort: pd.DataFrame,
    inputs: Mapping[str, Path],
    root: Path,
) -> pd.DataFrame:
    rows: list[dict] = []
    print("Counting uniquely traceable final-BLEACH peptide FASTA records...", flush=True)
    for i, row in cohort.iterrows():
        genome = str(row["Genome"])
        path = locate_final_peptide_fasta(genome, inputs["peptide_fastas"])
        if path is None:
            rows.append({
                "Genome": genome,
                "peptide_fasta": None,
                "peptide_fasta_realpath": None,
                "peptide_fasta_sha256": None,
                "final_peptide_records": np.nan,
                "status": "missing_no_unique_local_fasta",
            })
            continue
        count, digest = count_fasta_records(path)
        rows.append({
            "Genome": genome,
            "peptide_fasta": relative(path, root),
            "peptide_fasta_realpath": str(path.resolve()),
            "peptide_fasta_sha256": digest,
            "final_peptide_records": count,
            "status": "unique_exact_genome_filename_match",
        })
        if (i + 1) % 20 == 0 or i + 1 == len(cohort):
            print(f"  processed {i + 1}/{len(cohort)} FASTA mappings", flush=True)
    result = pd.DataFrame(rows)
    # The direct count is the analysis denominator. The legacy Exons field is
    # retained only as an audit comparison and is never relabelled as genes.
    exons = cohort.set_index("Genome")["Exons"]
    result["master_Exons_field"] = result["Genome"].map(exons)
    result["peptide_minus_Exons"] = result["final_peptide_records"] - result["master_Exons_field"]
    return result


def extract_raw_hmm_query_names(
    raw_manifest: pd.DataFrame,
    target_pfams: Iterable[str],
) -> dict[str, str]:
    targets = set(target_pfams)
    observed: dict[str, set[str]] = {pfam: set() for pfam in targets}
    for path_text in raw_manifest["raw_tblout_realpath"]:
        path = Path(path_text)
        with path.open("rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.lstrip().startswith("#") or "\t#\t" in line:
                    continue
                fields = line.rstrip().split("\t") if "\t" in line else line.split()
                for i, field in enumerate(fields):
                    match = re.fullmatch(r"(PF\d{5})(?:\.\d+)?", field.strip())
                    if match and match.group(1) in targets and i > 0:
                        observed[match.group(1)].add(fields[i - 1].strip())
        if all(len(names) == 1 for names in observed.values()):
            break
    result = {}
    for pfam, names in observed.items():
        if len(names) != 1:
            raise ValueError(f"Raw HMM query name unresolved/ambiguous for {pfam}: {sorted(names)}")
        result[pfam] = next(iter(names))
    return result


def add_verified_pfam_annotations(
    frame: pd.DataFrame,
    raw_query_names: Mapping[str, str],
    annotation_table: pd.DataFrame,
    fixed_pfam: str | None = None,
) -> pd.DataFrame:
    """Attach accession-level labels without inferring biological function."""
    required = {
        "accession", "short_name", "name", "source_database", "api_url",
        "http_status", "response_sha256",
    }
    missing = required - set(annotation_table.columns)
    if missing:
        raise ValueError(f"Verified annotation table lacks columns: {sorted(missing)}")
    if annotation_table["accession"].duplicated().any():
        raise ValueError("Verified annotation accessions are not unique")

    result = frame.copy()
    if fixed_pfam is not None:
        if "pfam" in result.columns and not result["pfam"].eq(fixed_pfam).all():
            raise ValueError("Fixed Pfam conflicts with result table")
        if "pfam" not in result.columns:
            result.insert(0, "pfam", fixed_pfam)
    if "pfam" not in result.columns:
        raise ValueError("Result table has no Pfam accession column")

    annotations = annotation_table.set_index("accession")
    accessions = set(result["pfam"].dropna().astype(str))
    absent = accessions - set(annotations.index)
    if absent:
        raise ValueError(f"Verified annotations missing result accessions: {sorted(absent)}")
    absent_raw = accessions - set(raw_query_names)
    if absent_raw:
        raise ValueError(f"Raw HMM query names missing result accessions: {sorted(absent_raw)}")
    selected_annotations = annotations.loc[sorted(accessions)]
    if not pd.to_numeric(selected_annotations["http_status"], errors="coerce").eq(200).all():
        raise ValueError("A selected InterPro annotation lacks a successful recorded response")
    if not selected_annotations["response_sha256"].astype(str).str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("A selected InterPro annotation lacks a valid response SHA-256")
    mismatches = {
        accession: (raw_query_names[accession], annotations.at[accession, "short_name"])
        for accession in accessions
        if raw_query_names[accession] != annotations.at[accession, "short_name"]
    }
    if mismatches:
        raise ValueError(f"Raw HMM and verified InterPro short names disagree: {mismatches}")

    insert_at = list(result.columns).index("pfam") + 1
    columns = [
        ("raw_hmm_query_name", None),
        ("verified_interpro_short_name", "short_name"),
        ("verified_interpro_name", "name"),
        ("annotation_source_database", "source_database"),
        ("annotation_api_url", "api_url"),
        ("annotation_http_status", "http_status"),
        ("annotation_response_sha256", "response_sha256"),
    ]
    for offset, (target, source) in enumerate(columns):
        values = (
            result["pfam"].map(raw_query_names)
            if source is None
            else result["pfam"].map(annotations[source])
        )
        result.insert(insert_at + offset, target, values)
    return result


def audit_consolidated_table(
    cohort: pd.DataFrame,
    raw_pfams: pd.DataFrame,
    inputs: Mapping[str, Path],
) -> pd.DataFrame:
    consolidated = pd.read_excel(inputs["final_s1"], sheet_name="TABLE_S1E-pfam_counts_with_met")
    consolidated = consolidated.set_index("Genome")
    missing = set(cohort["Genome"]) - set(consolidated.index)
    if missing:
        raise ValueError(f"Latest consolidated S1E misses canonical rows: {sorted(missing)}")
    c_pfams = [c for c in consolidated.columns if re.fullmatch(r"PF\d{5}", str(c))]
    common = sorted(set(c_pfams) & set(raw_pfams.columns))
    rows = []
    for genome in cohort["Genome"]:
        c = pd.to_numeric(consolidated.loc[genome, c_pfams], errors="coerce").fillna(0).to_numpy(float)
        raw_common = raw_pfams.loc[genome, common].to_numpy(float)
        con_common = pd.to_numeric(consolidated.loc[genome, common], errors="coerce").fillna(0).to_numpy(float)
        rows.append({
            "Genome": genome,
            "consolidated_pfam_sum": float(c.sum()),
            "raw_pfam_sum": float(raw_pfams.loc[genome].sum()),
            "consolidated_all_zero": bool(c.sum() == 0),
            "common_pfam_columns": len(common),
            "common_cell_mismatches": int(np.count_nonzero(raw_common != con_common)),
            "common_vectors_identical": bool(np.array_equal(raw_common, con_common)),
        })
    return pd.DataFrame(rows)


def rank_z(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    ranked = rankdata(values, method="average")
    sd = ranked.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return np.full(len(values), np.nan)
    return (ranked - ranked.mean()) / sd


def spearman_effect(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 5 or np.unique(x[valid]).size < 2 or np.unique(y[valid]).size < 2:
        return np.nan, np.nan
    result = stats.spearmanr(x[valid], y[valid])
    return float(result.statistic), float(result.pvalue)


def bootstrap_spearman_ci(
    x: np.ndarray,
    y: np.ndarray,
    n_boot: int,
    rng: np.random.Generator,
) -> tuple[float, float, int]:
    """Paired bootstrap on real observed pairs only."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    n = len(x)
    if n < 5 or np.unique(x).size < 2 or np.unique(y).size < 2 or n_boot <= 0:
        return np.nan, np.nan, 0
    estimates = np.empty(n_boot, dtype=float)
    kept = 0
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if np.unique(x[idx]).size < 2 or np.unique(y[idx]).size < 2:
            continue
        xr = rankdata(x[idx], method="average")
        yr = rankdata(y[idx], method="average")
        xr -= xr.mean()
        yr -= yr.mean()
        denominator = np.linalg.norm(xr) * np.linalg.norm(yr)
        if denominator > 0:
            estimates[kept] = float(np.dot(xr, yr) / denominator)
            kept += 1
    if kept < max(50, int(0.8 * n_boot)):
        return np.nan, np.nan, kept
    low, high = np.percentile(estimates[:kept], [2.5, 97.5])
    return float(low), float(high), kept


def standardize_column(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce").astype(float)
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.zeros(len(values)), index=values.index)
    return (values - values.mean()) / sd


def covariate_frame(
    metadata: pd.DataFrame,
    include_phylum: bool,
    include_peptide: bool,
    extra: pd.DataFrame | None = None,
) -> pd.DataFrame:
    cov = pd.DataFrame(index=metadata.index)
    cov["log10_assembly_nt"] = np.log10(pd.to_numeric(metadata["Nucleotides"], errors="coerce"))
    cov["BUSCO_percent"] = pd.to_numeric(metadata["BUSCOs-%present"], errors="coerce")
    cov["log10_total_pfam_hits"] = np.log10(pd.to_numeric(metadata["total_pfam_hits"], errors="coerce"))
    if include_peptide:
        cov["log10_final_peptide_records"] = np.log10(
            pd.to_numeric(metadata["final_peptide_records"], errors="coerce")
        )
    if include_phylum:
        dummies = pd.get_dummies(metadata["Phylum"], prefix="phylum", drop_first=True, dtype=float)
        cov = pd.concat([cov, dummies], axis=1)
    if extra is not None:
        cov = pd.concat([cov, extra], axis=1)
    for column in cov.columns:
        if not column.startswith("phylum_"):
            cov[column] = standardize_column(cov[column])
    return cov


def partial_rank_hc3(
    outcome: pd.Series,
    exposure: pd.Series,
    metadata: pd.DataFrame,
    include_phylum: bool = True,
    include_peptide: bool = True,
    extra_covariates: pd.DataFrame | None = None,
) -> dict:
    frame = pd.DataFrame({"outcome": outcome, "exposure": exposure}, index=metadata.index)
    cov = covariate_frame(metadata, include_phylum, include_peptide, extra_covariates)
    frame = pd.concat([frame, cov], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    if len(frame) < max(12, len(cov.columns) + 6):
        return {
            "n": len(frame), "effect": np.nan, "se": np.nan, "p_value": np.nan,
            "ci_low": np.nan, "ci_high": np.nan, "df_resid": np.nan,
            "condition_number": np.nan, "status": "insufficient_complete_cases",
            "complete_case_index": frame.index.tolist(),
        }
    y = rank_z(frame["outcome"].to_numpy())
    x_rank = rank_z(frame["exposure"].to_numpy())
    if np.isnan(y).all() or np.isnan(x_rank).all():
        return {
            "n": len(frame), "effect": np.nan, "se": np.nan, "p_value": np.nan,
            "ci_low": np.nan, "ci_high": np.nan, "df_resid": np.nan,
            "condition_number": np.nan, "status": "constant_rank_variable",
            "complete_case_index": frame.index.tolist(),
        }
    design = frame[cov.columns].astype(float).copy()
    design.insert(0, "exposure_rank_z", x_rank)
    design = sm.add_constant(design, has_constant="add")
    rank = np.linalg.matrix_rank(design.to_numpy())
    if rank < design.shape[1]:
        return {
            "n": len(frame), "effect": np.nan, "se": np.nan, "p_value": np.nan,
            "ci_low": np.nan, "ci_high": np.nan, "df_resid": np.nan,
            "condition_number": float(np.linalg.cond(design.to_numpy())),
            "status": "rank_deficient_design",
            "complete_case_index": frame.index.tolist(),
        }
    fit = sm.OLS(y, design).fit(cov_type="HC3")
    ci = fit.conf_int().loc["exposure_rank_z"]
    return {
        "n": int(fit.nobs),
        "effect": float(fit.params["exposure_rank_z"]),
        "se": float(fit.bse["exposure_rank_z"]),
        "p_value": float(fit.pvalues["exposure_rank_z"]),
        "ci_low": float(ci.iloc[0]),
        "ci_high": float(ci.iloc[1]),
        "df_resid": float(fit.df_resid),
        "condition_number": float(np.linalg.cond(design.to_numpy())),
        "status": "ok",
        "complete_case_index": frame.index.tolist(),
    }


def partial_rank_cluster(
    outcome: pd.Series,
    exposure: pd.Series,
    metadata: pd.DataFrame,
    cluster: pd.Series,
    include_phylum: bool = True,
    include_peptide: bool = True,
    extra_covariates: pd.DataFrame | None = None,
) -> dict:
    """Partial-rank OLS with coordinate-site clustered covariance."""
    frame = pd.DataFrame(
        {"outcome": outcome, "exposure": exposure, "cluster": cluster},
        index=metadata.index,
    )
    cov = covariate_frame(metadata, include_phylum, include_peptide, extra_covariates)
    frame = pd.concat([frame, cov], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    n_clusters = frame["cluster"].nunique()
    if len(frame) < max(12, len(cov.columns) + 6) or n_clusters < 10:
        return {
            "n": len(frame), "n_clusters": n_clusters, "effect": np.nan,
            "se": np.nan, "p_value": np.nan, "ci_low": np.nan,
            "ci_high": np.nan, "df_resid": np.nan, "condition_number": np.nan,
            "status": "insufficient_complete_cases_or_clusters",
            "complete_case_index": frame.index.tolist(),
        }
    y = rank_z(frame["outcome"].to_numpy())
    x_rank = rank_z(frame["exposure"].to_numpy())
    design = frame[cov.columns].astype(float).copy()
    design.insert(0, "exposure_rank_z", x_rank)
    design = sm.add_constant(design, has_constant="add")
    if np.linalg.matrix_rank(design.to_numpy()) < design.shape[1]:
        return {
            "n": len(frame), "n_clusters": n_clusters, "effect": np.nan,
            "se": np.nan, "p_value": np.nan, "ci_low": np.nan,
            "ci_high": np.nan, "df_resid": np.nan,
            "condition_number": float(np.linalg.cond(design.to_numpy())),
            "status": "rank_deficient_design",
            "complete_case_index": frame.index.tolist(),
        }
    fit = sm.OLS(y, design).fit(
        cov_type="cluster",
        cov_kwds={"groups": frame["cluster"], "use_correction": True},
    )
    ci = fit.conf_int().loc["exposure_rank_z"]
    return {
        "n": int(fit.nobs),
        "n_clusters": int(n_clusters),
        "effect": float(fit.params["exposure_rank_z"]),
        "se": float(fit.bse["exposure_rank_z"]),
        "p_value": float(fit.pvalues["exposure_rank_z"]),
        "ci_low": float(ci.iloc[0]),
        "ci_high": float(ci.iloc[1]),
        "df_resid": float(fit.df_resid),
        "condition_number": float(np.linalg.cond(design.to_numpy())),
        "status": "ok",
        "complete_case_index": frame.index.tolist(),
    }


def site_level_frame(
    outcome: pd.Series,
    exposure: pd.Series,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            "outcome": outcome,
            "exposure": exposure,
            "site_id": metadata["site_id"],
            "Phylum": metadata["Phylum"],
        },
        index=metadata.index,
    ).dropna(subset=["outcome", "exposure", "site_id"])
    exposure_nunique = frame.groupby("site_id")["exposure"].nunique(dropna=False)
    if (exposure_nunique > 1).any():
        bad = exposure_nunique[exposure_nunique > 1].index.tolist()[:5]
        raise ValueError(f"AEF exposure differs within exact coordinate sites: {bad}")
    site_phyla = metadata.groupby("site_id")["Phylum"].apply(
        lambda x: ";".join(sorted(set(map(str, x))))
    )
    site = frame.groupby("site_id").agg(
        outcome=("outcome", "mean"),
        exposure=("exposure", "first"),
        n_genomes=("outcome", "size"),
    )
    site["site_phyla"] = site.index.map(site_phyla)
    return site


def bh_adjust_by_group(
    frame: pd.DataFrame,
    group_columns: Sequence[str],
    p_column: str = "p_value",
    output_column: str = "selected_set_bh_q",
) -> pd.DataFrame:
    frame = frame.copy()
    frame[output_column] = np.nan
    for _, idx in frame.groupby(list(group_columns), dropna=False).groups.items():
        p = pd.to_numeric(frame.loc[idx, p_column], errors="coerce")
        valid = p.notna()
        if valid.any():
            q = multipletests(p[valid], method="fdr_bh")[1]
            frame.loc[p[valid].index, output_column] = q
    return frame


def make_representations(
    raw_pfams: pd.DataFrame,
    metadata: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    indexed = metadata.set_index("Genome")
    total = raw_pfams.sum(axis=1).astype(float)
    result = {
        "raw_count": raw_pfams.astype(float),
        "per_total_pfam_hit": raw_pfams.div(total, axis=0),
    }
    peptide = pd.to_numeric(indexed["final_peptide_records"], errors="coerce")
    result["per_final_peptide_record"] = raw_pfams.div(peptide, axis=0)
    return result


def select_gee_pairs(inputs: Mapping[str, Path], available_env: set[str]) -> pd.DataFrame:
    archived = pd.read_excel(inputs["headline"], sheet_name="Table_S2C_GEE_All_Correlations")
    archived["sig_fdr05"] = archived["sig_fdr05"].fillna(False).astype(bool)
    headline = archived[
        archived["sig_fdr05"] & archived["env_var"].isin(available_env)
    ][["pfam", "env_var", "spearman_r", "p_value", "p_fdr"]].copy()
    headline["selection_source"] = "archived_GEE_FDR_headline"

    pf00092 = archived[
        (archived["pfam"] == "PF00092") & archived["env_var"].isin(available_env)
    ].sort_values("p_value").head(3)[["pfam", "env_var", "spearman_r", "p_value", "p_fdr"]].copy()
    pf00092["selection_source"] = "archived_top3_PF00092_priority"
    selected = pd.concat([headline, pf00092], ignore_index=True)
    selected = selected.sort_values(["pfam", "env_var", "selection_source"]).drop_duplicates(
        ["pfam", "env_var"], keep="first"
    )
    selected = selected.rename(columns={
        "spearman_r": "archived_spearman_r",
        "p_value": "archived_p_value",
        "p_fdr": "archived_p_fdr",
    }).reset_index(drop=True)
    if selected.empty:
        raise ValueError("No archived headline/PF00092 GEE pairs available")
    return selected


def exact_gee_intersection(cohort: pd.DataFrame, gee_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    gee = pd.read_csv(gee_path)
    if gee["genome_id"].duplicated().any():
        raise ValueError("GEE genome_id values are not unique")
    metadata_columns = [
        "Genome", "Species", "Phylum", "Nucleotides", "BUSCOs-%present",
        "DD latitude", "DD longitude", "total_pfam_hits", "final_peptide_records",
        "is_uae",
    ]
    merged = cohort[metadata_columns].merge(
        gee,
        left_on="Genome",
        right_on="genome_id",
        how="inner",
        suffixes=("", "_archived_gee"),
        validate="one_to_one",
    )
    excluded = cohort.loc[~cohort["Genome"].isin(merged["Genome"]), ["Genome", "Species", "Phylum"]].copy()
    excluded["exclusion_reason"] = "not_exactly_present_in_archived_GEE_cohort"
    merged["geo_block_20deg"] = (
        np.floor((merged["DD latitude"] + 90.0) / GEO_BLOCK_DEGREES).astype(int).astype(str)
        + "_"
        + np.floor((merged["DD longitude"] + 180.0) / GEO_BLOCK_DEGREES).astype(int).astype(str)
    )
    return merged.set_index("Genome", drop=False), excluded


def residual_rank_correlations(
    matrix: np.ndarray,
    exposure: np.ndarray,
    covariates: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Vectorized Spearman/partial-Spearman correlations across matrix columns."""
    matrix = np.asarray(matrix, dtype=float)
    exposure = np.asarray(exposure, dtype=float)
    n, p = matrix.shape
    ranked_matrix = rankdata(matrix, axis=0, method="average")
    ranked_exposure = rankdata(exposure, method="average")

    if covariates is None:
        ranked_matrix -= ranked_matrix.mean(axis=0)
        ranked_exposure -= ranked_exposure.mean()
        df = n - 2
    else:
        c = np.asarray(covariates, dtype=float)
        if c.ndim != 2 or c.shape[0] != n:
            raise ValueError("Covariate matrix shape mismatch")
        c = np.column_stack([np.ones(n), c])
        rank_c = np.linalg.matrix_rank(c)
        projection = c @ np.linalg.pinv(c)
        ranked_matrix = ranked_matrix - projection @ ranked_matrix
        ranked_exposure = ranked_exposure - projection @ ranked_exposure
        df = n - rank_c - 1
    numerator = ranked_exposure @ ranked_matrix
    denominator = np.linalg.norm(ranked_exposure) * np.linalg.norm(ranked_matrix, axis=0)
    effects = np.divide(
        numerator,
        denominator,
        out=np.full(p, np.nan, dtype=float),
        where=denominator > 0,
    )
    effects = np.clip(effects, -1.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        t_stat = effects * np.sqrt(df / np.maximum(1e-15, 1.0 - effects ** 2))
    p_values = 2.0 * stats.t.sf(np.abs(t_stat), df=df)
    p_values[~np.isfinite(effects)] = np.nan
    return effects, p_values, int(df)


def screen_genomewide_combo(
    representation: pd.DataFrame,
    metadata: pd.DataFrame,
    env_columns: Sequence[str],
    eligible_pfams: Sequence[str],
    method: str,
) -> pd.DataFrame:
    records: list[pd.DataFrame] = []
    for env in env_columns:
        required = [env]
        include_cov = method == "quality_phylum_adjusted"
        include_peptide = include_cov
        if include_peptide:
            required.append("final_peptide_records")
        valid = metadata[required].replace([np.inf, -np.inf], np.nan).notna().all(axis=1)
        if method == "busco_ge50":
            valid &= metadata["BUSCOs-%present"] >= PRIMARY_BUSCO_THRESHOLD
        elif method == "busco_ge70":
            valid &= metadata["BUSCOs-%present"] >= STRICT_BUSCO_THRESHOLD
        submeta = metadata.loc[valid]
        if len(submeta) < 10:
            continue
        values = representation.loc[submeta.index, eligible_pfams].to_numpy(float)
        exposure = pd.to_numeric(submeta[env], errors="coerce").to_numpy(float)
        cov_array = None
        if include_cov:
            cov = covariate_frame(
                submeta,
                include_phylum=True,
                include_peptide=True,
                extra=None,
            )
            if cov.isna().any().any():
                raise ValueError(f"Unexpected missing covariates for {method}/{env}")
            cov_array = cov.to_numpy(float)
        effects, p_values, df = residual_rank_correlations(values, exposure, cov_array)
        records.append(pd.DataFrame({
            "method": method,
            "environment": env,
            "pfam": list(eligible_pfams),
            "n": len(submeta),
            "df": df,
            "effect": effects,
            "p_value": p_values,
        }))
    if not records:
        return pd.DataFrame()
    result = pd.concat(records, ignore_index=True)
    result["genomewide_bh_q"] = np.nan
    valid_p = result["p_value"].notna()
    if valid_p.any():
        result.loc[valid_p, "genomewide_bh_q"] = multipletests(
            result.loc[valid_p, "p_value"], method="fdr_bh"
        )[1]
    result["n_finite_tests_in_family"] = int(valid_p.sum())
    return result


def genomewide_screens(
    representations: Mapping[str, pd.DataFrame],
    gee_meta: pd.DataFrame,
    env_columns: Sequence[str],
    selected_pairs: pd.DataFrame,
    output_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    summary_rows: list[dict] = []
    candidate_frames: list[pd.DataFrame] = []
    lookup: dict[tuple[str, str, str, str], dict] = {}
    selected_keys = set(zip(selected_pairs["pfam"], selected_pairs["env_var"]))
    write_header = True
    methods = ["unadjusted", "quality_phylum_adjusted", "busco_ge50", "busco_ge70"]

    with gzip.open(output_path, "wt", encoding="utf-8", newline="") as handle:
        for representation_name, full_rep in representations.items():
            available = full_rep.loc[gee_meta.index]
            complete_rows = available.notna().all(axis=1)
            base_meta = gee_meta.loc[complete_rows].copy()
            base_rep = available.loc[complete_rows]
            prevalence = (base_rep > 0).sum(axis=0)
            eligible = prevalence[prevalence >= 10].index.tolist()
            if not eligible:
                raise ValueError(f"No prevalence-eligible Pfams for {representation_name}")
            print(
                f"Genome-wide screens: {representation_name}, n={len(base_meta)}, "
                f"eligible Pfams={len(eligible)}",
                flush=True,
            )
            for method in methods:
                result = screen_genomewide_combo(
                    base_rep,
                    base_meta,
                    env_columns,
                    eligible,
                    method,
                )
                if result.empty:
                    continue
                result.insert(0, "representation", representation_name)
                result.to_csv(handle, index=False, header=write_header)
                write_header = False

                for env, group in result.groupby("environment"):
                    finite = group["p_value"].notna()
                    top = group.loc[finite].sort_values("p_value").head(1)
                    summary_rows.append({
                        "representation": representation_name,
                        "method": method,
                        "environment": env,
                        "n_samples": int(group["n"].iloc[0]),
                        "finite_tests": int(finite.sum()),
                        "p_lt_0.001": int((group["p_value"] < 0.001).sum()),
                        "bh_q_lt_0.05": int((group["genomewide_bh_q"] < 0.05).sum()),
                        "top_pfam": top["pfam"].iloc[0] if not top.empty else None,
                        "top_effect": float(top["effect"].iloc[0]) if not top.empty else np.nan,
                        "min_p": float(top["p_value"].iloc[0]) if not top.empty else np.nan,
                        "min_q": float(top["genomewide_bh_q"].iloc[0]) if not top.empty else np.nan,
                    })

                selected_mask = [
                    (pfam, env) in selected_keys
                    for pfam, env in zip(result["pfam"], result["environment"])
                ]
                keep = (
                    (result["genomewide_bh_q"] < 0.05)
                    | (result["p_value"] < 0.001)
                    | np.asarray(selected_mask)
                    | (result["pfam"] == "PF00092")
                )
                candidate_frames.append(result.loc[keep].copy())
                for row in result.loc[np.asarray(selected_mask)].to_dict("records"):
                    lookup[(method, representation_name, row["pfam"], row["environment"])] = row
                del result
    candidates = pd.concat(candidate_frames, ignore_index=True) if candidate_frames else pd.DataFrame()
    return pd.DataFrame(summary_rows), candidates, lookup


def build_pev_from_reviewed_mapping(
    tree_path: Path,
    mapping_path: Path,
    canonical_genomes: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Build patristic PCoA axes from an integrity-reviewed tip mapping.

    The mapping is an input, not inferred here, because source phylogeny tables
    contain documented row/ID inconsistencies. Only rows explicitly marked as
    accepted/exact are used.
    """
    mapping = pd.read_csv(mapping_path)
    if {"tip", "Genome", "status"}.issubset(mapping.columns):
        accepted_words = {"accepted", "exact", "traceable", "ok"}
        accepted = mapping[
            mapping["status"].astype(str).str.lower().isin(accepted_words)
            & mapping["Genome"].isin(canonical_genomes)
        ].copy()
    elif {"rbcl_tip", "Genome", "safe_for_rbcl_pfam_analysis"}.issubset(mapping.columns):
        mapping = mapping.rename(columns={"rbcl_tip": "tip"}).copy()
        mapping["status"] = np.where(
            mapping["safe_for_rbcl_pfam_analysis"].fillna(False).astype(bool),
            "accepted",
            "excluded",
        )
        accepted = mapping[
            (mapping["status"] == "accepted")
            & mapping["tip"].notna()
            & mapping["Genome"].isin(canonical_genomes)
        ].copy()
    else:
        raise ValueError(
            "Tree mapping lacks either reviewed tip/Genome/status or integrity-manifest columns: "
            f"{mapping.columns.tolist()}"
        )
    mapping["tip"] = mapping["tip"].astype(str).str.replace("\xa0", "_", regex=False).str.replace(" ", "_", regex=False)
    accepted["tip"] = accepted["tip"].astype(str).str.replace("\xa0", "_", regex=False).str.replace(" ", "_", regex=False)
    if accepted["tip"].duplicated().any() or accepted["Genome"].duplicated().any():
        raise ValueError("Reviewed tree mapping is not one-tip-to-one-genome")
    raw_tree = tree_path.read_text().replace("\xa0", "_").replace(" ", "_")
    tree = Phylo.read(StringIO(raw_tree), "newick")
    original_terminal_lengths = np.asarray([
        float(t.branch_length) if t.branch_length is not None else np.nan
        for t in tree.get_terminals()
    ])
    # The deposited FastTree has many zero/near-zero terminal branches and is
    # unsuitable for a nonsingular Brownian covariance/PGLS. We therefore use
    # unit-edge topological distances for a transparent eigenvector sensitivity.
    # This controls hierarchical relatedness without pretending branch lengths
    # are calibrated evolutionary time.
    for clade in tree.find_clades(order="level"):
        if clade is not tree.root:
            clade.branch_length = 1.0
    terminals = {str(t.name): t for t in tree.get_terminals()}

    missing_tips = sorted(set(accepted["tip"]) - set(terminals))
    if missing_tips:
        raise ValueError(f"Reviewed mapping tips absent after Newick normalization: {missing_tips[:10]}")
    accepted = accepted[accepted["tip"].isin(terminals)].copy()
    terminal_list = [terminals[t] for t in accepted["tip"]]
    genomes = accepted["Genome"].tolist()
    n = len(terminal_list)
    if n < 30:
        raise ValueError(f"Too few accepted canonical tips for PEV analysis: {n}")
    distances = np.zeros((n, n), dtype=float)
    for i, tip_i in enumerate(terminal_list):
        for j in range(i, n):
            d = float(tree.distance(tip_i, terminal_list[j]))
            distances[i, j] = distances[j, i] = d
    jmat = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * jmat @ (distances ** 2) @ jmat
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    positive = eigenvalues > max(1e-12, eigenvalues[0] * 1e-10)
    positive_values = eigenvalues[positive]
    positive_vectors = eigenvectors[:, positive]
    fractions = positive_values / positive_values.sum()
    cumulative = np.cumsum(fractions)
    n_axes = len(fractions)
    broken_stick = np.asarray([
        sum(1.0 / j for j in range(i, n_axes + 1)) / n_axes
        for i in range(1, n_axes + 1)
    ])
    # Outcome-blind broken-stick selection: retain the leading consecutive axes
    # whose inertia exceeds the corresponding broken-stick expectation.
    leading_above = fractions > broken_stick
    k_broken_stick = next((i for i, value in enumerate(leading_above) if not value), n_axes)
    k_broken_stick = max(1, k_broken_stick)
    k5 = min(5, positive_vectors.shape[1])
    k10 = min(10, positive_vectors.shape[1])
    columns = [f"PEV{i+1}" for i in range(k10)]
    pev = pd.DataFrame(positive_vectors[:, :k10], index=genomes, columns=columns)
    pev.index.name = "Genome"
    eigen_audit = pd.DataFrame({
        "axis": [f"PEV{i+1}" for i in range(len(positive_values))],
        "eigenvalue": positive_values,
        "positive_inertia_fraction": fractions,
        "cumulative_positive_inertia": cumulative,
        "broken_stick_expectation": broken_stick,
        "selected_broken_stick": np.arange(len(positive_values)) < k_broken_stick,
        "selected_pev5": np.arange(len(positive_values)) < k5,
        "selected_pev10": np.arange(len(positive_values)) < k10,
    })
    details = {
        "n_mapped_canonical_tips": n,
        "n_positive_axes": len(positive_values),
        "k_broken_stick": k_broken_stick,
        "k5": k5,
        "k10": k10,
        "distance_basis": "unit_edge_topological_distance",
        "original_zero_or_nearzero_terminal_lengths": int(
            np.sum(np.isfinite(original_terminal_lengths) & (original_terminal_lengths <= 5e-9))
        ),
        "original_missing_terminal_lengths": int(np.sum(~np.isfinite(original_terminal_lengths))),
    }
    mapping_audit = mapping.copy()
    mapping_audit["used_in_canonical_pev"] = mapping_audit["Genome"].isin(pev.index)
    return pev, pd.concat([mapping_audit, eigen_audit], ignore_index=True, sort=False), details


def headline_gee_models(
    selected_pairs: pd.DataFrame,
    representations: Mapping[str, pd.DataFrame],
    gee_meta: pd.DataFrame,
    pev: pd.DataFrame | None,
    pev_details: Mapping[str, int] | None,
    n_boot: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    rows: list[dict] = []
    for pair in selected_pairs.itertuples(index=False):
        pfam = pair.pfam
        env = pair.env_var
        for rep_name, rep in representations.items():
            if pfam not in rep.columns:
                continue
            base = gee_meta.copy()
            outcome = rep.loc[base.index, pfam]
            exposure = pd.to_numeric(base[env], errors="coerce")
            complete = outcome.notna() & exposure.notna()
            submeta = base.loc[complete]
            x = outcome.loc[complete].to_numpy(float)
            y = exposure.loc[complete].to_numpy(float)
            effect, p_value = spearman_effect(x, y)
            ci_low, ci_high, kept = bootstrap_spearman_ci(x, y, n_boot, rng)
            common = {
                "pfam": pfam,
                "environment": env,
                "representation": rep_name,
                "selection_source": pair.selection_source,
                "archived_spearman_r": pair.archived_spearman_r,
                "archived_p_value": pair.archived_p_value,
                "archived_p_fdr": pair.archived_p_fdr,
            }
            rows.append({
                **common,
                "method": "spearman_full",
                "effect_type": "spearman_rho",
                "n": int(complete.sum()),
                "effect": effect,
                "se": np.nan,
                "p_value": p_value,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "bootstrap_replicates_retained": kept,
                "df_resid": max(0, int(complete.sum()) - 2),
                "condition_number": np.nan,
                "status": "ok" if np.isfinite(effect) else "not_estimable",
            })

            quality = partial_rank_hc3(
                outcome.loc[complete], exposure.loc[complete], submeta,
                include_phylum=True,
                include_peptide=True,
            )
            rows.append({
                **common,
                "method": "quality_phylum_adjusted_hc3",
                "effect_type": "standardized_partial_rank_beta",
                **quality,
                "bootstrap_replicates_retained": np.nan,
            })

            for threshold, method in [
                (PRIMARY_BUSCO_THRESHOLD, "busco_ge50_spearman"),
                (STRICT_BUSCO_THRESHOLD, "busco_ge70_spearman"),
            ]:
                filt = complete & (base["BUSCOs-%present"] >= threshold)
                fx = outcome.loc[filt].to_numpy(float)
                fy = exposure.loc[filt].to_numpy(float)
                fe, fp = spearman_effect(fx, fy)
                flo, fhi, fkept = bootstrap_spearman_ci(fx, fy, n_boot, rng)
                rows.append({
                    **common,
                    "method": method,
                    "effect_type": "spearman_rho",
                    "n": int(filt.sum()),
                    "effect": fe,
                    "se": np.nan,
                    "p_value": fp,
                    "ci_low": flo,
                    "ci_high": fhi,
                    "bootstrap_replicates_retained": fkept,
                    "df_resid": max(0, int(filt.sum()) - 2),
                    "condition_number": np.nan,
                    "status": "ok" if np.isfinite(fe) else "not_estimable",
                })

            if pev is not None and pev_details is not None:
                tree_ids = base.index.intersection(pev.index)
                tree_mask = complete & base.index.isin(tree_ids)
                tree_meta = base.loc[tree_mask]
                for key, method in [
                    ("k_broken_stick", "phylo_pev_brokenstick_quality_hc3"),
                    ("k5", "phylo_pev5_quality_hc3"),
                    ("k10", "phylo_pev10_quality_hc3"),
                ]:
                    k = int(pev_details[key])
                    extra = pev.loc[tree_meta.index, pev.columns[:k]]
                    model = partial_rank_hc3(
                        outcome.loc[tree_mask],
                        exposure.loc[tree_mask],
                        tree_meta,
                        include_phylum=False,
                        include_peptide=True,
                        extra_covariates=extra,
                    )
                    rows.append({
                        **common,
                        "method": method,
                        "effect_type": "standardized_partial_rank_beta",
                        **model,
                        "bootstrap_replicates_retained": np.nan,
                    })
            else:
                rows.append({
                    **common,
                    "method": "phylo_pev_status_gate",
                    "effect_type": "not_estimated",
                    "n": 0,
                    "effect": np.nan,
                    "se": np.nan,
                    "p_value": np.nan,
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "bootstrap_replicates_retained": np.nan,
                    "df_resid": np.nan,
                    "condition_number": np.nan,
                    "status": "blocked_without_integrity_reviewed_tip_mapping",
                })
    result = pd.DataFrame(rows)
    result = bh_adjust_by_group(result, ["method", "representation"])
    return result


def grouped_permutation_matrix(
    values: np.ndarray,
    groups: Sequence[object],
    n_perm: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, dict]:
    """Shuffle observed values only within supplied structural groups."""
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups, dtype=object)
    group_indices = [np.flatnonzero(groups == group) for group in pd.unique(groups)]
    result = np.empty((n_perm, len(values)), dtype=float)
    for b in range(n_perm):
        permuted = values.copy()
        for idx in group_indices:
            if len(idx) > 1:
                permuted[idx] = values[rng.permutation(idx)]
        result[b] = permuted
    sizes = np.asarray([len(idx) for idx in group_indices])
    audit = {
        "n_groups": len(group_indices),
        "n_multisample_groups": int((sizes > 1).sum()),
        "n_samples_in_multisample_groups": int(sizes[sizes > 1].sum()),
        "max_group_size": int(sizes.max()) if len(sizes) else 0,
    }
    return result, audit


def structured_permutation_tests(
    selected_pairs: pd.DataFrame,
    representations: Mapping[str, pd.DataFrame],
    gee_meta: pd.DataFrame,
    n_perm: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pair_rows: list[dict] = []
    null_by_key: dict[tuple[str, str, str, str], np.ndarray] = {}
    observed_by_key: dict[tuple[str, str, str, str], tuple[float, int]] = {}
    structures = {"within_phylum": "Phylum", "within_geo20": "geo_block_20deg"}

    for env, env_pairs in selected_pairs.groupby("env_var"):
        for rep_name, rep in representations.items():
            pfams = [pf for pf in env_pairs["pfam"] if pf in rep.columns]
            if not pfams:
                continue
            available = rep.loc[gee_meta.index, pfams]
            valid = pd.to_numeric(gee_meta[env], errors="coerce").notna() & available.notna().all(axis=1)
            meta = gee_meta.loc[valid]
            y_rank = rankdata(pd.to_numeric(meta[env], errors="coerce").to_numpy(float), method="average")
            for structure_name, group_col in structures.items():
                permuted_y, audit = grouped_permutation_matrix(
                    y_rank, meta[group_col].astype(str).to_numpy(), n_perm, rng
                )
                permuted_y -= permuted_y.mean(axis=1, keepdims=True)
                perm_norm = np.linalg.norm(permuted_y, axis=1)
                for pfam in pfams:
                    x_rank = rankdata(rep.loc[meta.index, pfam].to_numpy(float), method="average")
                    x_rank -= x_rank.mean()
                    denominator = perm_norm * np.linalg.norm(x_rank)
                    null = np.divide(
                        permuted_y @ x_rank,
                        denominator,
                        out=np.full(n_perm, np.nan),
                        where=denominator > 0,
                    )
                    observed = float(np.dot(y_rank - y_rank.mean(), x_rank) / (
                        np.linalg.norm(y_rank - y_rank.mean()) * np.linalg.norm(x_rank)
                    )) if np.linalg.norm(x_rank) > 0 else np.nan
                    finite_null = null[np.isfinite(null)]
                    empirical_p = (
                        (1 + np.sum(np.abs(finite_null) >= abs(observed))) / (1 + len(finite_null))
                        if np.isfinite(observed) and len(finite_null) else np.nan
                    )
                    q025, q975 = (
                        np.percentile(finite_null, [2.5, 97.5]) if len(finite_null) else (np.nan, np.nan)
                    )
                    row = {
                        "pfam": pfam,
                        "environment": env,
                        "representation": rep_name,
                        "structure": structure_name,
                        "n": len(meta),
                        "observed_spearman_rho": observed,
                        "empirical_two_sided_p": empirical_p,
                        "null_mean": float(np.mean(finite_null)) if len(finite_null) else np.nan,
                        "null_sd": float(np.std(finite_null, ddof=1)) if len(finite_null) > 1 else np.nan,
                        "null_q025": q025,
                        "null_q975": q975,
                        "permutations": len(finite_null),
                        **audit,
                    }
                    pair_rows.append(row)
                    key = (rep_name, structure_name, pfam, env)
                    null_by_key[key] = null
                    observed_by_key[key] = (observed, len(meta))

    pair_result = pd.DataFrame(pair_rows)
    pair_result = bh_adjust_by_group(
        pair_result,
        ["representation", "structure"],
        p_column="empirical_two_sided_p",
        output_column="selected_set_empirical_bh_q",
    )

    count_rows: list[dict] = []
    for rep_name in representations:
        for structure_name in structures:
            keys = [k for k in null_by_key if k[0] == rep_name and k[1] == structure_name]
            if not keys:
                continue
            observed_ps = []
            null_ps = []
            for key in keys:
                observed, n = observed_by_key[key]
                if not np.isfinite(observed) or n <= 2:
                    continue
                observed_t = observed * math.sqrt((n - 2) / max(1e-15, 1 - observed ** 2))
                observed_ps.append(2 * stats.t.sf(abs(observed_t), df=n - 2))
                null = null_by_key[key]
                null_t = null * np.sqrt((n - 2) / np.maximum(1e-15, 1 - null ** 2))
                null_ps.append(2 * stats.t.sf(np.abs(null_t), df=n - 2))
            if not null_ps:
                continue
            null_count = (np.vstack(null_ps) < 0.05).sum(axis=0)
            observed_count = int((np.asarray(observed_ps) < 0.05).sum())
            empirical_count_p = (1 + np.sum(null_count >= observed_count)) / (1 + len(null_count))
            count_rows.append({
                "representation": rep_name,
                "structure": structure_name,
                "selected_pairs_tested": len(observed_ps),
                "observed_nominal_p_lt_0.05": observed_count,
                "null_mean_count": float(null_count.mean()),
                "null_sd_count": float(null_count.std(ddof=1)),
                "null_q025_count": float(np.percentile(null_count, 2.5)),
                "null_q975_count": float(np.percentile(null_count, 97.5)),
                "empirical_count_p": float(empirical_count_p),
                "permutations": len(null_count),
            })
    return pair_result, pd.DataFrame(count_rows)


def cliffs_delta(group1: np.ndarray, group0: np.ndarray) -> float:
    group1 = np.asarray(group1, dtype=float)
    group0 = np.asarray(group0, dtype=float)
    if len(group1) == 0 or len(group0) == 0:
        return np.nan
    differences = group1[:, None] - group0[None, :]
    return float((np.sum(differences > 0) - np.sum(differences < 0)) / differences.size)


def bootstrap_mean_ratio_ci(
    group1: np.ndarray,
    group0: np.ndarray,
    n_boot: int,
    rng: np.random.Generator,
) -> tuple[float, float, int]:
    group1 = np.asarray(group1, dtype=float)
    group0 = np.asarray(group0, dtype=float)
    if len(group1) < 2 or len(group0) < 2 or n_boot <= 0 or group0.mean() <= 0:
        return np.nan, np.nan, 0
    ratios = np.empty(n_boot, dtype=float)
    kept = 0
    for _ in range(n_boot):
        a = group1[rng.integers(0, len(group1), len(group1))].mean()
        b = group0[rng.integers(0, len(group0), len(group0))].mean()
        if b > 0:
            ratios[kept] = a / b
            kept += 1
    if kept < max(50, int(0.8 * n_boot)):
        return np.nan, np.nan, kept
    low, high = np.percentile(ratios[:kept], [2.5, 97.5])
    return float(low), float(high), kept


def model_sample_summary(
    model: Mapping[str, object],
    metadata: pd.DataFrame,
    outcome: pd.Series,
    exposure: pd.Series,
) -> dict:
    """Summarize exactly the complete cases used by a fitted partial-rank model."""
    index = pd.Index(model["complete_case_index"])
    if len(index) != int(model["n"]):
        raise AssertionError("Model n does not match its complete-case index")
    sample_meta = metadata.loc[index]
    sample_outcome = outcome.loc[index]
    sample_exposure = exposure.loc[index].astype(int)
    group_1 = sample_outcome[sample_exposure == 1]
    group_0 = sample_outcome[sample_exposure == 0]
    return {
        "n_total": int(len(index)),
        "n_sites": int(sample_meta["site_id"].nunique()),
        "n_uae": int((sample_exposure == 1).sum()),
        "n_non_uae": int((sample_exposure == 0).sum()),
        "uae_mean": float(group_1.mean()) if len(group_1) else np.nan,
        "non_uae_mean": float(group_0.mean()) if len(group_0) else np.nan,
        "uae_median": float(group_1.median()) if len(group_1) else np.nan,
        "non_uae_median": float(group_0.median()) if len(group_0) else np.nan,
    }


def pf00092_uae_analysis(
    representations: Mapping[str, pd.DataFrame],
    cohort_meta: pd.DataFrame,
    pev: pd.DataFrame | None,
    pev_details: Mapping[str, int] | None,
    n_boot: int,
    n_perm: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    rows: list[dict] = []
    metadata = cohort_meta.set_index("Genome", drop=False)
    for rep_name, rep in representations.items():
        outcome = rep.loc[metadata.index, "PF00092"]
        valid = outcome.notna()
        meta = metadata.loc[valid]
        values = outcome.loc[valid]
        uae = meta["is_uae"].astype(int)
        g1 = values[uae == 1].to_numpy(float)
        g0 = values[uae == 0].to_numpy(float)
        mw = stats.mannwhitneyu(g1, g0, alternative="two-sided") if len(g1) and len(g0) else None
        fold = float(g1.mean() / g0.mean()) if len(g0) and g0.mean() > 0 else np.nan
        fold_low, fold_high, kept = bootstrap_mean_ratio_ci(g1, g0, n_boot, rng)
        common = {
            "representation": rep_name,
            "n_total": len(meta),
            "n_sites": int(meta["site_id"].nunique()),
            "n_uae": int((uae == 1).sum()),
            "n_non_uae": int((uae == 0).sum()),
        }
        rows.append({
            **common,
            "method": "genome_level_group_descriptive",
            "effect_type": "mean_ratio_UAE_over_nonUAE",
            "effect": fold,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p_value": np.nan,
            "naive_unclustered_mannwhitney_p": float(mw.pvalue) if mw is not None else np.nan,
            "secondary_effect_cliffs_delta": cliffs_delta(g1, g0),
            "uae_mean": float(g1.mean()) if len(g1) else np.nan,
            "non_uae_mean": float(g0.mean()) if len(g0) else np.nan,
            "uae_median": float(np.median(g1)) if len(g1) else np.nan,
            "non_uae_median": float(np.median(g0)) if len(g0) else np.nan,
            "bootstrap_replicates_retained": np.nan,
            "status": "descriptive_only_repeated_coordinate_sites",
        })

        site = site_level_frame(values, uae, meta)
        site_label = site["exposure"].astype(int)
        s1 = site.loc[site_label == 1, "outcome"].to_numpy(float)
        s0 = site.loc[site_label == 0, "outcome"].to_numpy(float)
        site_test = stats.mannwhitneyu(s1, s0, alternative="two-sided") if len(s1) and len(s0) else None
        site_ratio = float(s1.mean() / s0.mean()) if len(s0) and s0.mean() > 0 else np.nan
        slo, shi, ske = bootstrap_mean_ratio_ci(s1, s0, n_boot, rng)
        rows.append({
            "representation": rep_name,
            "method": "site_mean_group_comparison",
            "effect_type": "site_mean_ratio_UAE_over_nonUAE",
            "n_total": len(site),
            "n_sites": len(site),
            "n_uae": int((site_label == 1).sum()),
            "n_non_uae": int((site_label == 0).sum()),
            "effect": site_ratio,
            "ci_low": slo,
            "ci_high": shi,
            "p_value": float(site_test.pvalue) if site_test is not None else np.nan,
            "secondary_effect_cliffs_delta": cliffs_delta(s1, s0),
            "uae_mean": float(s1.mean()),
            "non_uae_mean": float(s0.mean()),
            "uae_median": float(np.median(s1)),
            "non_uae_median": float(np.median(s0)),
            "bootstrap_replicates_retained": ske,
            "status": "unique_coordinate_site_inference",
        })

        quality = partial_rank_hc3(
            values,
            uae,
            meta,
            include_phylum=True,
            include_peptide=True,
        )
        quality_sample = model_sample_summary(quality, meta, values, uae)
        rows.append({
            "representation": rep_name,
            **quality_sample,
            "method": "genome_quality_phylum_hc3_noncluster",
            "effect_type": "standardized_partial_rank_beta",
            "effect": quality["effect"],
            "ci_low": quality["ci_low"],
            "ci_high": quality["ci_high"],
            "p_value": quality["p_value"],
            "secondary_effect_cliffs_delta": np.nan,
            "bootstrap_replicates_retained": np.nan,
            "status": quality["status"],
        })
        site_cluster = partial_rank_cluster(
            values, uae, meta, cluster=meta["site_id"],
            include_phylum=True, include_peptide=True,
        )
        site_cluster_sample = model_sample_summary(site_cluster, meta, values, uae)
        if site_cluster_sample["n_sites"] != int(site_cluster["n_clusters"]):
            raise AssertionError("Clustered-model site count does not match complete cases")
        rows.append({
            "representation": rep_name,
            **site_cluster_sample,
            "method": "quality_phylum_sitecluster",
            "effect_type": "standardized_partial_rank_beta",
            "effect": site_cluster["effect"],
            "ci_low": site_cluster["ci_low"],
            "ci_high": site_cluster["ci_high"],
            "p_value": site_cluster["p_value"],
            "secondary_effect_cliffs_delta": np.nan,
            "bootstrap_replicates_retained": np.nan,
            "condition_number": site_cluster["condition_number"],
            "status": site_cluster["status"],
        })

        for threshold, method in [
            (PRIMARY_BUSCO_THRESHOLD, "busco_ge50_group_comparison"),
            (STRICT_BUSCO_THRESHOLD, "busco_ge70_group_comparison"),
        ]:
            keep = meta["BUSCOs-%present"] >= threshold
            fm = meta.loc[keep]
            fv = values.loc[keep]
            fu = fm["is_uae"].astype(int)
            f1 = fv[fu == 1].to_numpy(float)
            f0 = fv[fu == 0].to_numpy(float)
            test = stats.mannwhitneyu(f1, f0, alternative="two-sided") if len(f1) and len(f0) else None
            fratio = float(f1.mean() / f0.mean()) if len(f0) and f0.mean() > 0 else np.nan
            flo, fhi, fkept = bootstrap_mean_ratio_ci(f1, f0, n_boot, rng)
            rows.append({
                "representation": rep_name,
                "method": f"genome_{method}_descriptive",
                "effect_type": "mean_ratio_UAE_over_nonUAE",
                "n_total": len(fm),
                "n_uae": int((fu == 1).sum()),
                "n_non_uae": int((fu == 0).sum()),
                "effect": fratio,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p_value": np.nan,
                "naive_unclustered_mannwhitney_p": float(test.pvalue) if test is not None else np.nan,
                "secondary_effect_cliffs_delta": cliffs_delta(f1, f0),
                "uae_mean": float(f1.mean()) if len(f1) else np.nan,
                "non_uae_mean": float(f0.mean()) if len(f0) else np.nan,
                "uae_median": float(np.median(f1)) if len(f1) else np.nan,
                "non_uae_median": float(np.median(f0)) if len(f0) else np.nan,
                "bootstrap_replicates_retained": np.nan,
                "status": "descriptive_only_repeated_coordinate_sites",
            })

            site_filtered = site_level_frame(fv, fu, fm)
            slabel = site_filtered["exposure"].astype(int)
            sf1 = site_filtered.loc[slabel == 1, "outcome"].to_numpy(float)
            sf0 = site_filtered.loc[slabel == 0, "outcome"].to_numpy(float)
            stest = stats.mannwhitneyu(sf1, sf0, alternative="two-sided") if len(sf1) and len(sf0) else None
            sratio = float(sf1.mean() / sf0.mean()) if len(sf0) and sf0.mean() > 0 else np.nan
            sflo, sfhi, sfkept = bootstrap_mean_ratio_ci(sf1, sf0, n_boot, rng)
            rows.append({
                "representation": rep_name,
                "method": method.replace("group_comparison", "site_mean_group_comparison"),
                "effect_type": "site_mean_ratio_UAE_over_nonUAE",
                "n_total": len(site_filtered),
                "n_sites": len(site_filtered),
                "n_uae": int((slabel == 1).sum()),
                "n_non_uae": int((slabel == 0).sum()),
                "effect": sratio,
                "ci_low": sflo,
                "ci_high": sfhi,
                "p_value": float(stest.pvalue) if stest is not None else np.nan,
                "secondary_effect_cliffs_delta": cliffs_delta(sf1, sf0),
                "uae_mean": float(sf1.mean()),
                "non_uae_mean": float(sf0.mean()),
                "uae_median": float(np.median(sf1)),
                "non_uae_median": float(np.median(sf0)),
                "bootstrap_replicates_retained": sfkept,
                "status": "unique_coordinate_site_inference",
            })

        # UAE labels are shuffled within trusted phyla. This preserves the
        # observed lineage composition of the nine UAE genomes.
        uae_values = site_label.to_numpy(float)
        permuted, audit = grouped_permutation_matrix(
            uae_values, site["site_phyla"].to_numpy(), n_perm, rng
        )
        observed_diff = float(s1.mean() - s0.mean())
        null_diff = np.full(n_perm, np.nan)
        val_array = site["outcome"].to_numpy(float)
        for b in range(n_perm):
            label = permuted[b].astype(int)
            if (label == 1).any() and (label == 0).any():
                null_diff[b] = val_array[label == 1].mean() - val_array[label == 0].mean()
        finite = null_diff[np.isfinite(null_diff)]
        empirical_p = (
            (1 + np.sum(np.abs(finite) >= abs(observed_diff))) / (1 + len(finite))
            if len(finite) else np.nan
        )
        rows.append({
            "representation": rep_name,
            "method": "site_phylum_composition_label_permutation",
            "effect_type": "mean_difference_UAE_minus_nonUAE",
            "n_total": len(site),
            "n_sites": len(site),
            "n_uae": int((site_label == 1).sum()),
            "n_non_uae": int((site_label == 0).sum()),
            "effect": observed_diff,
            "ci_low": float(np.percentile(finite, 2.5)) if len(finite) else np.nan,
            "ci_high": float(np.percentile(finite, 97.5)) if len(finite) else np.nan,
            "p_value": empirical_p,
            "secondary_effect_cliffs_delta": np.nan,
            "uae_mean": float(s1.mean()),
            "non_uae_mean": float(s0.mean()),
            "uae_median": float(np.median(s1)),
            "non_uae_median": float(np.median(s0)),
            "bootstrap_replicates_retained": np.nan,
            "status": "ok" if len(finite) else "not_estimable",
            **audit,
        })

        if pev is not None and pev_details is not None:
            tree_ids = meta.index.intersection(pev.index)
            tm = meta.loc[tree_ids]
            tv = values.loc[tree_ids]
            tu = tm["is_uae"].astype(int)
            for key, method in [
                ("k_broken_stick", "phylo_pev_brokenstick_quality_hc3"),
                ("k5", "phylo_pev5_quality_hc3"),
                ("k10", "phylo_pev10_quality_hc3"),
            ]:
                extra = pev.loc[tree_ids, pev.columns[:int(pev_details[key])]]
                model = partial_rank_hc3(
                    tv,
                    tu,
                    tm,
                    include_phylum=False,
                    include_peptide=True,
                    extra_covariates=extra,
                )
                model_sample = model_sample_summary(model, tm, tv, tu)
                rows.append({
                    "representation": rep_name,
                    **model_sample,
                    "method": f"genome_{method}_noncluster",
                    "effect_type": "standardized_partial_rank_beta",
                    "effect": model["effect"],
                    "ci_low": model["ci_low"],
                    "ci_high": model["ci_high"],
                    "p_value": model["p_value"],
                    "secondary_effect_cliffs_delta": np.nan,
                    "bootstrap_replicates_retained": np.nan,
                    "condition_number": model["condition_number"],
                    "status": model["status"],
                })
                cluster_model = partial_rank_cluster(
                    tv, tu, tm, cluster=tm["site_id"], include_phylum=False,
                    include_peptide=True, extra_covariates=extra,
                )
                cluster_sample = model_sample_summary(cluster_model, tm, tv, tu)
                if cluster_sample["n_sites"] != int(cluster_model["n_clusters"]):
                    raise AssertionError("Tree clustered-model site count does not match complete cases")
                rows.append({
                    "representation": rep_name,
                    **cluster_sample,
                    "method": method.replace("quality_hc3", "quality_sitecluster"),
                    "effect_type": "standardized_partial_rank_beta",
                    "effect": cluster_model["effect"],
                    "ci_low": cluster_model["ci_low"],
                    "ci_high": cluster_model["ci_high"],
                    "p_value": cluster_model["p_value"],
                    "secondary_effect_cliffs_delta": np.nan,
                    "bootstrap_replicates_retained": np.nan,
                    "condition_number": cluster_model["condition_number"],
                    "status": cluster_model["status"],
                })
    result = pd.DataFrame(rows)
    result = bh_adjust_by_group(
        result, ["method"], p_column="p_value", output_column="representation_set_bh_q"
    )
    return result


def archived_aef_sensitivity(
    aef: pd.DataFrame,
    cohort_meta: pd.DataFrame,
    representations: Mapping[str, pd.DataFrame],
    pev: pd.DataFrame | None,
    pev_details: Mapping[str, int] | None,
    n_boot: int,
    n_perm: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sensitivity only; latent axes are not assigned physical meanings."""
    overlap = cohort_meta.index.intersection(aef["Genome"])
    embedding = aef.set_index("Genome").loc[overlap]
    meta = cohort_meta.loc[overlap]
    explicit_pairs = [
        ("PF13411", "A18"),
        ("PF01638", "A52"),
        ("PF01638", "A53"),
        ("PF10988", "A36"),
    ]
    axis_columns = [f"A{i:02d}" for i in range(64) if f"A{i:02d}" in embedding.columns]
    pairs = explicit_pairs + [("PF00092", axis) for axis in axis_columns]
    rows: list[dict] = []
    permutation_rows: list[dict] = []
    for pfam, axis in pairs:
        if axis not in embedding or axis not in axis_columns:
            continue
        for rep_name, rep in representations.items():
            if pfam not in rep:
                continue
            outcome = rep.loc[overlap, pfam]
            exposure = pd.to_numeric(embedding.loc[overlap, axis], errors="coerce")
            valid = outcome.notna() & exposure.notna()
            m = meta.loc[valid]
            x = outcome.loc[valid].to_numpy(float)
            y = exposure.loc[valid].to_numpy(float)
            rho, p = spearman_effect(x, y)
            rows.append({
                "pfam": pfam,
                "latent_axis": axis,
                "representation": rep_name,
                "method": "genome_level_spearman_descriptive",
                "n": len(m),
                "n_sites": int(m["site_id"].nunique()),
                "effect": rho,
                "p_value": p,
                "status": "descriptive_only_repeated_site_exposures",
            })
            site = site_level_frame(outcome, exposure, meta)
            site_rho, site_p = spearman_effect(site["outcome"].to_numpy(), site["exposure"].to_numpy())
            site_low, site_high, site_boot_kept = bootstrap_spearman_ci(
                site["outcome"].to_numpy(), site["exposure"].to_numpy(), n_boot, rng
            )
            rows.append({
                "pfam": pfam,
                "latent_axis": axis,
                "representation": rep_name,
                "method": "site_mean_spearman",
                "n": len(site),
                "n_sites": len(site),
                "effect": site_rho,
                "se": np.nan,
                "p_value": site_p,
                "ci_low": site_low,
                "ci_high": site_high,
                "bootstrap_replicates_retained": site_boot_kept,
                "status": "unique_coordinate_site_inference",
            })
            quality = partial_rank_hc3(
                outcome.loc[valid], exposure.loc[valid], m,
                include_phylum=True, include_peptide=True,
            )
            rows.append({
                "pfam": pfam,
                "latent_axis": axis,
                "representation": rep_name,
                "method": "genome_quality_phylum_hc3_noncluster",
                "n": quality["n"],
                "effect": quality["effect"],
                "se": quality["se"],
                "p_value": quality["p_value"],
                "ci_low": quality["ci_low"],
                "ci_high": quality["ci_high"],
                "df_resid": quality["df_resid"],
                "condition_number": quality["condition_number"],
                "status": quality["status"],
            })
            clustered = partial_rank_cluster(
                outcome.loc[valid], exposure.loc[valid], m,
                cluster=m["site_id"], include_phylum=True, include_peptide=True,
            )
            rows.append({
                "pfam": pfam,
                "latent_axis": axis,
                "representation": rep_name,
                "method": "quality_phylum_sitecluster",
                "n": clustered["n"],
                "n_sites": clustered["n_clusters"],
                "effect": clustered["effect"],
                "se": clustered["se"],
                "p_value": clustered["p_value"],
                "ci_low": clustered["ci_low"],
                "ci_high": clustered["ci_high"],
                "df_resid": clustered["df_resid"],
                "condition_number": clustered["condition_number"],
                "status": clustered["status"],
            })
            for threshold, method in [
                (PRIMARY_BUSCO_THRESHOLD, "busco_ge50_spearman"),
                (STRICT_BUSCO_THRESHOLD, "busco_ge70_spearman"),
            ]:
                filt = valid & (meta.reindex(valid.index)["BUSCOs-%present"] >= threshold)
                # `valid` and cohort metadata share the same canonical index.
                ids = valid.index[filt]
                fe, fp = spearman_effect(
                    outcome.loc[ids].to_numpy(float),
                    exposure.loc[ids].to_numpy(float),
                )
                rows.append({
                    "pfam": pfam,
                    "latent_axis": axis,
                    "representation": rep_name,
                    "method": f"genome_{method}_descriptive",
                    "n": len(ids),
                    "effect": fe,
                    "p_value": fp,
                    "status": "descriptive_only_repeated_site_exposures" if np.isfinite(fe) else "not_estimable",
                })
                filtered_outcome = outcome.loc[ids]
                filtered_exposure = exposure.loc[ids]
                filtered_meta = meta.loc[ids]
                site_filtered = site_level_frame(filtered_outcome, filtered_exposure, filtered_meta)
                sfe, sfp = spearman_effect(
                    site_filtered["outcome"].to_numpy(), site_filtered["exposure"].to_numpy()
                )
                sflo, sfhi, sfkept = bootstrap_spearman_ci(
                    site_filtered["outcome"].to_numpy(),
                    site_filtered["exposure"].to_numpy(),
                    n_boot, rng,
                )
                rows.append({
                    "pfam": pfam,
                    "latent_axis": axis,
                    "representation": rep_name,
                    "method": f"site_mean_{method}",
                    "n": len(site_filtered),
                    "n_sites": len(site_filtered),
                    "effect": sfe,
                    "se": np.nan,
                    "p_value": sfp,
                    "ci_low": sflo,
                    "ci_high": sfhi,
                    "bootstrap_replicates_retained": sfkept,
                    "status": "unique_coordinate_site_inference" if np.isfinite(sfe) else "not_estimable",
                })
            if pev is not None and pev_details is not None:
                tree_ids = m.index.intersection(pev.index)
                tm = m.loc[tree_ids]
                for key, method in [
                    ("k_broken_stick", "phylo_pev_brokenstick_quality_hc3"),
                    ("k5", "phylo_pev5_quality_hc3"),
                    ("k10", "phylo_pev10_quality_hc3"),
                ]:
                    extra = pev.loc[tree_ids, pev.columns[:int(pev_details[key])]]
                    model = partial_rank_hc3(
                        outcome.loc[tree_ids], exposure.loc[tree_ids], tm,
                        include_phylum=False, include_peptide=True, extra_covariates=extra,
                    )
                    rows.append({
                        "pfam": pfam,
                        "latent_axis": axis,
                        "representation": rep_name,
                        "method": f"genome_{method}_noncluster",
                        "n": model["n"],
                        "effect": model["effect"],
                        "se": model["se"],
                        "p_value": model["p_value"],
                        "ci_low": model["ci_low"],
                        "ci_high": model["ci_high"],
                        "df_resid": model["df_resid"],
                        "condition_number": model["condition_number"],
                        "status": model["status"],
                    })
                    cluster_model = partial_rank_cluster(
                        outcome.loc[tree_ids], exposure.loc[tree_ids], tm,
                        cluster=tm["site_id"], include_phylum=False,
                        include_peptide=True, extra_covariates=extra,
                    )
                    rows.append({
                        "pfam": pfam,
                        "latent_axis": axis,
                        "representation": rep_name,
                        "method": method.replace("quality_hc3", "quality_sitecluster"),
                        "n": cluster_model["n"],
                        "n_sites": cluster_model["n_clusters"],
                        "effect": cluster_model["effect"],
                        "se": cluster_model["se"],
                        "p_value": cluster_model["p_value"],
                        "ci_low": cluster_model["ci_low"],
                        "ci_high": cluster_model["ci_high"],
                        "df_resid": cluster_model["df_resid"],
                        "condition_number": cluster_model["condition_number"],
                        "status": cluster_model["status"],
                    })

            # Unique-site empirical null. AEF vectors are permuted as intact
            # site units within the exact set of phyla represented at each site.
            y_rank = rankdata(site["exposure"].to_numpy(float), method="average")
            x_rank = rankdata(site["outcome"].to_numpy(float), method="average")
            perm_y, audit = grouped_permutation_matrix(
                y_rank, site["site_phyla"].to_numpy(), n_perm, rng
            )
            y0 = y_rank - y_rank.mean()
            x0 = x_rank - x_rank.mean()
            obs = float(np.dot(x0, y0) / (np.linalg.norm(x0) * np.linalg.norm(y0)))
            perm_y -= perm_y.mean(axis=1, keepdims=True)
            den = np.linalg.norm(perm_y, axis=1) * np.linalg.norm(x0)
            null = np.divide(perm_y @ x0, den, out=np.full(n_perm, np.nan), where=den > 0)
            finite = null[np.isfinite(null)]
            empirical = (1 + np.sum(np.abs(finite) >= abs(obs))) / (1 + len(finite))
            permutation_rows.append({
                "pfam": pfam,
                "latent_axis": axis,
                "representation": rep_name,
                "analysis_level": "unique_coordinate_site_mean",
                "n": len(site),
                "n_sites": len(site),
                "n_genomes_contributing": int(site["n_genomes"].sum()),
                "observed_spearman_rho": obs,
                "empirical_two_sided_p": empirical,
                "null_mean": float(finite.mean()),
                "null_sd": float(finite.std(ddof=1)),
                "null_q025": float(np.percentile(finite, 2.5)),
                "null_q975": float(np.percentile(finite, 97.5)),
                "permutations": len(finite),
                **audit,
            })
    result = pd.DataFrame(rows)
    result = bh_adjust_by_group(
        result, ["representation", "method"], p_column="p_value",
        output_column="selected_AEF_set_bh_q",
    )
    perm = pd.DataFrame(permutation_rows)
    perm = bh_adjust_by_group(
        perm, ["representation"], p_column="empirical_two_sided_p",
        output_column="selected_AEF_empirical_bh_q",
    )
    return result, perm


def build_robust_candidate_table(
    aef_results: pd.DataFrame,
    aef_permutations: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict] = []
    for (pfam, axis), group in aef_results.groupby(["pfam", "latent_axis"]):
        effects = group.loc[group["effect"].notna(), "effect"]
        raw_full = group[
            (group["method"] == "site_mean_spearman")
            & (group["representation"] == "raw_count")
        ]
        reference_sign = np.sign(raw_full["effect"].iloc[0]) if not raw_full.empty else np.nan
        sign_consistent = bool(
            len(effects)
            and np.isfinite(reference_sign)
            and reference_sign != 0
            and (np.sign(effects) == reference_sign).all()
        )
        checks = {
            "raw_full": group[(group.method == "site_mean_spearman") & (group.representation == "raw_count")],
            "total_full": group[(group.method == "site_mean_spearman") & (group.representation == "per_total_pfam_hit")],
            "peptide_full": group[(group.method == "site_mean_spearman") & (group.representation == "per_final_peptide_record")],
            "total_quality": group[(group.method == "quality_phylum_sitecluster") & (group.representation == "per_total_pfam_hit")],
            "total_busco50": group[(group.method == "site_mean_busco_ge50_spearman") & (group.representation == "per_total_pfam_hit")],
            "total_tree": group[(group.method == "phylo_pev_brokenstick_quality_sitecluster") & (group.representation == "per_total_pfam_hit")],
        }
        check_values = {
            name: bool(not frame.empty and frame["selected_AEF_set_bh_q"].iloc[0] < 0.05)
            for name, frame in checks.items()
        }
        perm = aef_permutations[
            (aef_permutations["pfam"] == pfam)
            & (aef_permutations["latent_axis"] == axis)
            & (aef_permutations["representation"] == "per_total_pfam_hit")
        ]
        phylum_empirical = bool(
            not perm.empty and perm["selected_AEF_empirical_bh_q"].iloc[0] < 0.05
        )
        required = [
            check_values["raw_full"], check_values["total_full"],
            check_values["peptide_full"], check_values["total_quality"],
            check_values["total_busco50"], check_values["total_tree"],
            phylum_empirical,
        ]
        rows.append({
            "pfam": pfam,
            "latent_axis": axis,
            "reference_raw_effect": float(raw_full["effect"].iloc[0]) if not raw_full.empty else np.nan,
            "direction_consistent_all_estimable_checks": sign_consistent,
            **{f"selected_set_q_lt_0.05_{k}": v for k, v in check_values.items()},
            "empirical_selected_set_q_lt_0.05_total_hit_within_phylum": phylum_empirical,
            "checks_passed": int(sum(required)),
            "checks_required": len(required),
            "robust_candidate_all_required_checks": bool(sign_consistent and all(required)),
            "interpretation_boundary": "archived latent-feature association; no physical meaning or adaptation inferred",
        })
    return pd.DataFrame(rows).sort_values(
        ["robust_candidate_all_required_checks", "checks_passed", "reference_raw_effect"],
        ascending=[False, False, False],
    )


def create_robustness_figure(
    aef_results: pd.DataFrame,
    aef_permutations: pd.DataFrame,
    output_pdf: Path,
    output_svg: Path,
) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

    mpl.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica"],
        "font.size": 6,
        "axes.labelsize": 6,
        "axes.titlesize": 6,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.linewidth": 0.25,
        "xtick.major.width": 0.25,
        "ytick.major.width": 0.25,
        "xtick.major.size": 2,
        "ytick.major.size": 2,
        "axes.labelpad": 1,
        "xtick.major.pad": 1,
        "ytick.major.pad": 1,
    })
    colors = [
        (0.0, 0.3, 0.7), (0.3, 0.5, 0.9), (0.7, 0.8, 0.95),
        (0.95, 0.95, 0.95), (0.95, 0.8, 0.7), (0.9, 0.4, 0.3), (0.7, 0.1, 0.1),
    ]
    cmap = LinearSegmentedColormap.from_list("journal_diverging", colors)

    raw = aef_results[
        (aef_results["representation"] == "raw_count")
        & (aef_results["method"] == "site_mean_spearman")
    ].sort_values("p_value")
    explicit = raw[raw["pfam"] != "PF00092"][["pfam", "latent_axis"]].drop_duplicates()
    pf00092 = raw[raw["pfam"] == "PF00092"].head(8)[["pfam", "latent_axis"]]
    selected = pd.concat([explicit, pf00092], ignore_index=True).drop_duplicates().head(12)
    selected_keys = set(map(tuple, selected.to_numpy()))
    plot_data = aef_results[
        [(p, a) in selected_keys for p, a in zip(aef_results.pfam, aef_results.latent_axis)]
    ].copy()
    method_order = [
        "site_mean_spearman", "quality_phylum_sitecluster",
        "site_mean_busco_ge50_spearman", "phylo_pev_brokenstick_quality_sitecluster",
    ]
    rep_order = ["raw_count", "per_total_pfam_hit", "per_final_peptide_record"]
    columns = [(m, r) for m in method_order for r in rep_order]
    row_order = [tuple(x) for x in selected.to_numpy()]
    matrix = np.full((len(row_order), len(columns)), np.nan)
    qmatrix = np.full_like(matrix, np.nan)
    for i, key in enumerate(row_order):
        for j, (method, rep) in enumerate(columns):
            hit = plot_data[
                (plot_data.pfam == key[0]) & (plot_data.latent_axis == key[1])
                & (plot_data.method == method) & (plot_data.representation == rep)
            ]
            if not hit.empty:
                matrix[i, j] = hit.effect.iloc[0]
                qmatrix[i, j] = hit.selected_AEF_set_bh_q.iloc[0]

    fig = plt.figure(figsize=(7.0, 5.7))
    gs = fig.add_gridspec(
        2, 2, height_ratios=[3.2, 1.4], width_ratios=[4.2, 1.8],
        hspace=0.58, wspace=0.30,
    )
    ax_heat = fig.add_subplot(gs[0, :])
    vmax = max(0.4, float(np.nanmax(np.abs(matrix))))
    image = ax_heat.imshow(matrix, aspect="auto", cmap=cmap, norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
    short_method = {
        "site_mean_spearman": "Site", "quality_phylum_sitecluster": "Qual",
        "site_mean_busco_ge50_spearman": "B50",
        "phylo_pev_brokenstick_quality_sitecluster": "Phy",
    }
    short_rep = {
        "raw_count": "R", "per_total_pfam_hit": "T",
        "per_final_peptide_record": "P",
    }
    ax_heat.set_xticks(
        range(len(columns)),
        [f"{short_method[m]}\n{short_rep[r]}" for m, r in columns],
        rotation=0, ha="center",
    )
    ax_heat.set_yticks(range(len(row_order)), [f"{p}/{a}" for p, a in row_order])
    ax_heat.set_title("A  Priority archived latent-feature associations across robustness checks", loc="left", fontweight="bold")
    ax_heat.set_xlabel(
        "Qual, coordinate-clustered quality model; B50, BUSCO >=50% site mean; "
        "Phy, broken-stick topology + site cluster\nR, raw; T, total-hit ratio; P, peptide-record ratio"
    )
    for boundary in [2.5, 5.5, 8.5]:
        ax_heat.axvline(boundary, color="black", linewidth=0.25)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if np.isfinite(qmatrix[i, j]) and qmatrix[i, j] < 0.05:
                ax_heat.text(j, i, "*", ha="center", va="center", color="black")
    cbar = fig.colorbar(image, ax=ax_heat, fraction=0.018, pad=0.01)
    cbar.set_label("Rank association effect")
    for spine in ax_heat.spines.values():
        spine.set_linewidth(0.25)

    ax_null = fig.add_subplot(gs[1, 0])
    selected_pair_mask = pd.Series(
        [
            (pfam, axis) in selected_keys
            for pfam, axis in zip(aef_permutations.pfam, aef_permutations.latent_axis)
        ],
        index=aef_permutations.index,
    )
    null = aef_permutations[
        (aef_permutations.representation == "per_total_pfam_hit")
        & selected_pair_mask
    ].copy()
    null["key"] = list(zip(null.pfam, null.latent_axis))
    null = null.set_index("key").reindex(row_order).dropna(subset=["observed_spearman_rho"])
    ypos = np.arange(len(null))
    ax_null.hlines(ypos, null.null_q025, null.null_q975, color="#7f7f7f", linewidth=0.25)
    sig = null.empirical_two_sided_p < 0.05
    ax_null.scatter(null.observed_spearman_rho[~sig], ypos[~sig], s=9, facecolor="white", edgecolor="#1f4e79", linewidth=0.25, zorder=3)
    ax_null.scatter(null.observed_spearman_rho[sig], ypos[sig], s=9, color="#b22222", linewidth=0, zorder=3)
    ax_null.axvline(0, color="black", linewidth=0.25)
    ax_null.set_yticks(ypos, [f"{p}/{a}" for p, a in null.index])
    ax_null.invert_yaxis()
    ax_null.set_xlabel("Site-mean rho (point); site-block null 95% interval")
    ax_null.set_title("B  Structured empirical null", loc="left", fontweight="bold")

    ax_count = fig.add_subplot(gs[1, 1])
    count = (
        aef_results[aef_results[["pfam", "latent_axis"]].apply(tuple, axis=1).isin(selected_keys)]
        .groupby(["method", "representation"])["p_value"]
        .apply(lambda x: int((x < 0.05).sum()))
        .reindex(pd.MultiIndex.from_tuples(columns), fill_value=0)
    )
    bar_colors = ["#355c7d", "#6c8ebf", "#b7c9e2"] * len(method_order)
    ax_count.bar(np.arange(len(count)), count.to_numpy(), color=bar_colors, linewidth=0)
    ax_count.set_xticks(np.arange(len(count)), [short_rep[r] for _, r in columns], rotation=0)
    ax_count.set_ylabel("Nominal p<0.05 pairs")
    ax_count.set_title("C  Retention by check", loc="left", fontweight="bold")
    ax_count.set_ylim(0, max(1, int(count.max()) + 1))
    ax_count.grid(axis="y", color="#cccccc", linewidth=0.25)

    for ax in [ax_heat, ax_null, ax_count]:
        ax.tick_params(width=0.25)
    fig.savefig(output_pdf, format="pdf", bbox_inches="tight", transparent=True, edgecolor="none")
    fig.savefig(output_svg, format="svg", bbox_inches="tight", transparent=True, edgecolor="none")
    plt.close(fig)


def write_manuscript_summary(
    path: Path,
    cohort: pd.DataFrame,
    raw_manifest: pd.DataFrame,
    peptide_manifest: pd.DataFrame,
    consolidated_audit: pd.DataFrame,
    pev_details: Mapping[str, object],
    robust_candidates: pd.DataFrame,
    uae_results: pd.DataFrame,
    aef_results: pd.DataFrame,
    gee_audit: Mapping[str, object],
    args: argparse.Namespace,
) -> None:
    phylum_counts = cohort["Phylum"].value_counts().to_dict()
    peptide_n = int(peptide_manifest["final_peptide_records"].notna().sum())
    zero_rows = int(consolidated_audit["consolidated_all_zero"].sum())
    robust = robust_candidates[robust_candidates["robust_candidate_all_required_checks"]]
    pfam_annotations = (
        robust_candidates[
            ["pfam", "raw_hmm_query_name", "verified_interpro_short_name", "verified_interpro_name"]
        ]
        .drop_duplicates("pfam")
        .set_index("pfam")
    )
    uae_focus = uae_results[
        (uae_results.representation == "per_total_pfam_hit")
        & (uae_results.method.isin(["site_mean_group_comparison", "quality_phylum_sitecluster", "phylo_pev_brokenstick_quality_sitecluster", "site_phylum_composition_label_permutation"]))
    ]
    priority = aef_results[
        (aef_results.pfam.isin(["PF13411", "PF01638", "PF10988"]))
        & (aef_results.representation == "per_total_pfam_hit")
        & (aef_results.method.isin([
            "site_mean_spearman", "quality_phylum_sitecluster",
            "site_mean_busco_ge50_spearman",
            "phylo_pev_brokenstick_quality_sitecluster",
        ]))
    ]
    tree_condition = pd.to_numeric(
        aef_results.loc[
            aef_results.method.str.contains("phylo_pev", na=False),
            "condition_number",
        ],
        errors="coerce",
    )

    if str(gee_audit.get("status", "")).startswith("VERIFIED_REEXTRACTION"):
        gee_status_line = (
            f"- A verified GEE baseline containing {int(gee_audit['baseline_rows'])} genomes "
            "was ingested and audited; GEE association estimates are generated by the separate "
            "exact-ID workflow."
        )
    else:
        gee_status_line = (
            "- GEE association estimates are generated by the separate exact-ID workflow; this "
            "robustness workflow does not read archived derived GEE values."
        )

    lines = [
        "# Manuscript-ready quantitative robustness summary",
        "",
        "## Data and inferential status",
        "",
        f"- The canonical AEF cohort contains {len(cohort)} genomes: {phylum_counts}.",
        f"- These genomes occupy {cohort['site_id'].nunique()} exact coordinate sites; {(cohort['site_id'].value_counts() > 1).sum()} sites are repeated (maximum {cohort['site_id'].value_counts().max()} genomes). Genome-level correlations are descriptive, while inference uses site means or coordinate-clustered covariance.",
        f"- Raw HMM-search outputs yielded {int(raw_manifest.parsed_pfam_hits.sum()):,} sequence-level Pfam-hit records across {raw_manifest.shape[0]} genomes; every raw total passed the authoritative metadata checksum.",
        f"- The consolidated S1E matrix contains {zero_rows} false all-zero rows in this cohort. All analyses use raw reconstructed counts, not those zeros.",
        f"- Final-BLEACH amino-acid FASTAs provide a directly counted peptide/translated-ORF denominator for {peptide_n}/{len(cohort)} genomes. This is a search-space proxy, not a documented predicted-gene count; exact total predicted proteins remain unavailable.",
        f"- The explicit tree sensitivity uses {pev_details.get('n_mapped_canonical_tips', 0)} accession-gated AEF genomes and unit-edge topological distances. The outcome-blind broken-stick rule retained {pev_details.get('k_broken_stick', 'NA')} axes; 5- and 10-axis models are reported as sensitivity checks. It is a phylogenetic-eigenvector analysis, not PGLS. The deposited tree has {pev_details.get('original_zero_or_nearzero_terminal_lengths', 'NA')} zero/near-zero terminal branches. Across fitted Pfam–AEF tree models, the maximum design condition number was {tree_condition.max():.2f}.",
        gee_status_line,
        "",
        "## Chosen robustness remedies",
        "",
        "- Raw count models remain visible. Total-Pfam-hit ratios are a compositional sensitivity, not a replacement primary outcome; denominator-adjusted rank models provide the complementary check.",
        "- Continuous BUSCO and assembly/search-space covariates are primary. BUSCO >=50% and >=70% exclusions are threshold sensitivities because no universal biological cutoff exists for these divergent algal lineages.",
        "- Rank-based effects are used because Pfam counts are sparse and overdispersed. Primary genome-level regression sensitivities use covariance clustered by exact coordinate site; HC3 fits are retained only as non-clustered descriptive checks. No Gaussian count PGLS is claimed.",
        f"- Structured empirical p-values use {args.permutations:,} shuffles of intact coordinate-site units within site phylum-composition strata; bootstrap intervals use {args.bootstraps:,} resamples of unique site means.",
        "",
        "## Archived latent-feature sensitivity",
        "",
    ]
    if robust.empty:
        lines.append("No tested Pfam–latent-axis pair passed every prespecified normalization, quality, BUSCO, topology, and within-phylum empirical-null check. The revised manuscript should not call any AEF pair uniquely robust.")
    else:
        pairs = ", ".join(
            f"{r.pfam} ({r.raw_hmm_query_name}; {r.verified_interpro_name})/{r.latent_axis}"
            for r in robust.itertuples()
        )
        lines.append(f"Pairs passing every prespecified check: {pairs}. These remain latent-feature associations; no physical meaning, adaptation, or mechanism is inferred.")
    lines.extend(["", "Priority numerical results (total-hit normalization):", ""])
    if priority.empty:
        lines.append("- No priority result was estimable.")
    else:
        for row in priority.itertuples():
            q = row.selected_AEF_set_bh_q
            annotation = pfam_annotations.loc[row.pfam]
            lines.append(
                f"- {row.pfam} ({annotation.raw_hmm_query_name}; {annotation.verified_interpro_name})/"
                f"{row.latent_axis}, {row.method}: n={int(row.n)}, effect={row.effect:.3f}, "
                f"95% CI [{row.ci_low:.3f}, {row.ci_high:.3f}], p={row.p_value:.3g}, "
                f"selected-set q={q:.3g}."
            )
    lines.extend(["", "## PF00092 Arabian Gulf comparison", ""])
    for row in uae_focus.itertuples():
        interval = (
            f", 95% CI [{row.ci_low:.3f}, {row.ci_high:.3f}]"
            if np.isfinite(row.ci_low) and np.isfinite(row.ci_high)
            else ""
        )
        lines.append(
            f"- {row.method}: n={int(row.n_total)} ({int(row.n_uae)} UAE), "
            f"{row.effect_type}={row.effect:.3f}{interval}, p={row.p_value:.3g}, status={row.status}."
        )
    corrected = pfam_annotations.loc[["PF10988", "PF13411"]]
    corrected_labels = "; ".join(
        f"{accession} as {row.raw_hmm_query_name} ({row.verified_interpro_name})"
        for accession, row in corrected.iterrows()
    )
    lines.append(
        "Accession-level raw-HMM and verified InterPro records identify "
        f"{corrected_labels}. The prior NAD(P)-binding and peroxidase labels do not map "
        "to these accessions and should be replaced."
    )
    lines.extend([
        "",
        "## Text suitable for the Results/response letter",
        "",
        "We reconstructed Pfam abundances directly from per-genome HMM-search outputs and removed false zero profiles introduced during table consolidation. We evaluated raw counts, counts per total Pfam hit, and—where a final amino-acid FASTA was uniquely traceable—counts per final peptide record. Because 126 genomes represented 90 coordinate sites, genome-level correlations are reported descriptively and inference is based on site-mean analyses or coordinate-clustered covariance. Continuous assembly size, BUSCO completeness, total Pfam hits, peptide-search-space size, and phylum were included in rank-based models; BUSCO thresholds were retained as sensitivity checks rather than treated as universal quality cutoffs. A topology-based phylogenetic eigenvector analysis was applied to accession-gated tree matches because zero-length branches precluded a defensible Gaussian PGLS covariance. Structured null tests permuted intact site units within site phylum-composition strata. These analyses prioritize statistical stability and do not convert association into evidence of selection, adaptation, or mechanism.",
        "",
        "The archived AEF dimensions are treated strictly as latent geospatial descriptors. Any retained Pfam–dimension result is evidence of covariance with a latent feature, not decoding of an environmental variable. GEE association estimates are produced by the separate exact-ID extraction and validation workflow.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    if args.permutations < 99 or args.bootstraps < 99:
        raise ValueError("At least 99 real-data permutations/bootstraps are required")
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    output_dir = script_path.parent
    run_id = now_stamp()
    rng = np.random.default_rng(args.seed)  # resampling observed data only
    inputs = load_inputs(root)

    cohort, aef, legacy = load_canonical_cohort(inputs)
    verified_gee, verified_gee_audit = validate_verified_gee_ingest(
        args.verified_gee_csv,
        args.verified_gee_provenance,
        set(cohort["Genome"]),
    )
    raw_pfams, raw_manifest = reconstruct_raw_pfam_matrix(cohort, inputs, root)
    priority_annotations = pd.read_csv(inputs["priority_annotations"])
    raw_query_names = extract_raw_hmm_query_names(raw_manifest, ANALYSIS_PFAMS)
    annotation_audit = add_verified_pfam_annotations(
        pd.DataFrame({"pfam": list(ANALYSIS_PFAMS)}),
        raw_query_names,
        priority_annotations,
    )
    annotation_audit["raw_hmm_interpro_short_name_exact_match"] = (
        annotation_audit["raw_hmm_query_name"]
        == annotation_audit["verified_interpro_short_name"]
    )
    peptide_manifest = derive_peptide_denominators(cohort, inputs, root)
    cohort = cohort.merge(
        peptide_manifest[["Genome", "final_peptide_records", "status"]].rename(columns={"status": "peptide_status"}),
        on="Genome", how="left", validate="one_to_one",
    )
    cohort["total_pfam_hits"] = cohort["Genome"].map(raw_pfams.sum(axis=1))
    cohort["site_id"] = (
        cohort["DD latitude"].round(10).map(lambda x: f"{x:.10f}")
        + "|"
        + cohort["DD longitude"].round(10).map(lambda x: f"{x:.10f}")
    )
    uae_meta = pd.read_excel(inputs["final_s1"], sheet_name="Table S1A New Species Meta")
    uae_ids = set(uae_meta["Genome ID"].dropna().astype(str))
    cohort["is_uae"] = cohort["Genome"].isin(uae_ids).astype(int)
    if cohort["is_uae"].sum() != len(uae_ids):
        raise ValueError("UAE identifiers do not map one-to-one to canonical cohort")

    consolidated_audit = audit_consolidated_table(cohort, raw_pfams, inputs)
    tree_map_path = args.tree_map.resolve() if args.tree_map else inputs["integrity_manifest"]
    pev, phylo_audit, pev_details = build_pev_from_reviewed_mapping(
        inputs["tree"], tree_map_path, set(cohort["Genome"])
    )

    # Timestamped validation/provenance outputs are always written, including
    # validate-only runs, so every scientific matrix remains traceable.
    paths = {
        "raw_manifest": output_dir / f"raw_pfam_source_manifest_{run_id}.csv",
        "raw_matrix": output_dir / f"reconstructed_raw_pfam_counts_{run_id}.csv.gz",
        "peptide_manifest": output_dir / f"final_peptide_denominator_manifest_{run_id}.csv",
        "consolidated_audit": output_dir / f"consolidated_pfam_integrity_audit_{run_id}.csv",
        "phylo_audit": output_dir / f"phylogenetic_mapping_and_eigen_audit_{run_id}.csv",
        "pev": output_dir / f"topology_phylogenetic_eigenvectors_{run_id}.csv",
        "gee_gate": output_dir / f"GEE_STATUS_GATE_{run_id}.csv",
        "annotation_audit": output_dir / f"priority_pfam_annotation_audit_{run_id}.csv",
        "aef_results": output_dir / f"archived_AEF_priority_robustness_{run_id}.csv",
        "aef_perm": output_dir / f"archived_AEF_within_phylum_null_{run_id}.csv",
        "uae": output_dir / f"PF00092_UAE_robustness_{run_id}.csv",
        "candidates": output_dir / f"robust_candidate_table_{run_id}.csv",
        "figure_pdf": output_dir / f"Figure_Robustness_Sensitivity_{run_id}.pdf",
        "figure_svg": output_dir / f"Figure_Robustness_Sensitivity_{run_id}.svg",
        "summary": output_dir / f"manuscript_ready_quantitative_summary_{run_id}.md",
        "manifest": output_dir / f"run_manifest_{run_id}.json",
    }
    raw_manifest.to_csv(paths["raw_manifest"], index=False)
    raw_pfams.reset_index().to_csv(paths["raw_matrix"], index=False, compression="gzip")
    peptide_manifest.to_csv(paths["peptide_manifest"], index=False)
    consolidated_audit.to_csv(paths["consolidated_audit"], index=False)
    phylo_audit.to_csv(paths["phylo_audit"], index=False)
    pev.reset_index().to_csv(paths["pev"], index=False)
    annotation_audit.to_csv(paths["annotation_audit"], index=False)

    gate_rows = []
    for variable in GEE_ENV_COLUMNS:
        gate_rows.append({
            "variable": variable,
            "status": verified_gee_audit["status"],
            "source_sha256": verified_gee_audit.get("csv_sha256", "not supplied"),
            "required_remedy": verified_gee_audit.get(
                "follow_up",
                "run corrected Earth Engine extraction, then pass --verified-gee-csv and --verified-gee-provenance",
            ),
            "baseline_rows": verified_gee_audit.get("baseline_rows", pd.NA),
            "nonmissing_rows": verified_gee_audit.get("complete_counts", {}).get(variable, pd.NA),
        })
    pd.DataFrame(gate_rows).to_csv(paths["gee_gate"], index=False)

    run_manifest = {
        "run_id": run_id,
        "script_version": SCRIPT_VERSION,
        "script": str(script_path),
        "command": sys.argv,
        "parameters": vars(args) | {
            "primary_busco_threshold": PRIMARY_BUSCO_THRESHOLD,
            "strict_busco_threshold": STRICT_BUSCO_THRESHOLD,
            "geo_block_degrees": GEO_BLOCK_DEGREES,
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "statsmodels": statsmodels.__version__,
            "biopython": Bio.__version__,
        },
        "git": git_metadata(root),
        "inputs": {k: input_record(v, root, with_hash=(v.is_file() and v.stat().st_size < 100_000_000)) for k, v in inputs.items()},
        "data_counts": {
            "canonical_genomes": len(cohort),
            "raw_pfam_accessions": raw_pfams.shape[1],
            "raw_pfam_hits": int(raw_pfams.to_numpy().sum()),
            "traceable_peptide_denominators": int(peptide_manifest.final_peptide_records.notna().sum()),
            "false_zero_consolidated_rows": int(consolidated_audit.consolidated_all_zero.sum()),
            "tree_mapped_AEF_genomes": int(pev_details["n_mapped_canonical_tips"]),
            "UAE_genomes": int(cohort.is_uae.sum()),
            "unique_coordinate_sites": int(cohort.site_id.nunique()),
            "repeated_coordinate_sites": int((cohort.site_id.value_counts() > 1).sum()),
            "max_genomes_per_coordinate_site": int(cohort.site_id.value_counts().max()),
            "priority_pfam_annotations_verified": int(len(annotation_audit)),
        },
        "phylogenetic_method": pev_details,
        "gee_status": verified_gee_audit,
        "outputs": {k: str(v.resolve()) for k, v in paths.items()},
    }

    if args.validate_only:
        run_manifest["status"] = "validation_complete_statistics_skipped"
        paths["manifest"].write_text(json.dumps(run_manifest, indent=2, default=str), encoding="utf-8")
        print(json.dumps(run_manifest["data_counts"], indent=2), flush=True)
        return 0

    representations = make_representations(raw_pfams, cohort)
    indexed_meta = cohort.set_index("Genome", drop=False)
    uae_results = pf00092_uae_analysis(
        representations, cohort, pev, pev_details,
        args.bootstraps, args.permutations, rng,
    )
    aef_results, aef_perm = archived_aef_sensitivity(
        aef, indexed_meta, representations, pev, pev_details,
        args.bootstraps, args.permutations, rng,
    )
    candidates = build_robust_candidate_table(aef_results, aef_perm)

    uae_results = add_verified_pfam_annotations(
        uae_results, raw_query_names, priority_annotations, fixed_pfam="PF00092"
    )
    aef_results = add_verified_pfam_annotations(
        aef_results, raw_query_names, priority_annotations
    )
    aef_perm = add_verified_pfam_annotations(
        aef_perm, raw_query_names, priority_annotations
    )
    candidates = add_verified_pfam_annotations(
        candidates, raw_query_names, priority_annotations
    )

    uae_results.to_csv(paths["uae"], index=False)
    aef_results.to_csv(paths["aef_results"], index=False)
    aef_perm.to_csv(paths["aef_perm"], index=False)
    candidates.to_csv(paths["candidates"], index=False)
    create_robustness_figure(aef_results, aef_perm, paths["figure_pdf"], paths["figure_svg"])
    write_manuscript_summary(
        paths["summary"], cohort, raw_manifest, peptide_manifest,
        consolidated_audit, pev_details, candidates, uae_results,
        aef_results, verified_gee_audit, args,
    )
    run_manifest["status"] = "complete"
    paths["manifest"].write_text(json.dumps(run_manifest, indent=2, default=str), encoding="utf-8")
    print("Analysis complete. Outputs:", flush=True)
    for key, value in paths.items():
        if value.exists():
            print(f"  {key}: {value.resolve()}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

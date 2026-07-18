#!/usr/bin/env python3
"""Test canonical MIDAS coordinating residues in retained PF00092 hit proteins.

This is a fresh, sequence-level reannotation of real retained data. It does not
reconstruct the missing original PF00092.31 domtblout files. The program:

1. Selects the authenticated 126-genome cohort from the reconciliation manifest.
2. Reads only the per-genome HMMER tblout file named in that manifest.
3. Extracts target protein IDs reported for PF00092.31 (VWA).
4. Recovers those exact proteins from the retained peptide FASTA named in the
   peptide-denominator manifest.
5. Searches them with the current official Pfam PF00092 HMM downloaded from the
   InterPro API and applies that model's gathering threshold.
6. Reads residues aligned to five structurally homologous HMM states:
   D1/S2/S3/T4/D5 = states 7/9/11/79/112 in PF00092.35.

The state mapping is checked at runtime against two experimentally resolved
controls downloaded from RCSB PDB: 1SHU has a complete canonical MIDAS sequence
(D,S,S,T,D), whereas 1AO3 (human VWF A3) has the vestigial pattern
(D,S,S,S,T). A complete sequence signature supports MIDAS compatibility only;
it does not demonstrate metal binding, adhesion, tissue cohesion, substrate
attachment, or any macroalgal biological function.

No synthetic, simulated, randomly generated, or hardcoded result values are
used. Every reported count derives from retained input files or computations on
their sequences. Output directories and filenames are timestamped and created
with overwrite protection.

Created: 2026-07-15
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import platform
import re
import sys
import time
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Sequence, Tuple

import pandas as pd
import pyhmmer


SCRIPT_PATH = Path(__file__).resolve()
ROOT = SCRIPT_PATH.parents[2]
REVISION = ROOT / "ISCIENCE_REVISION_20260711"

RECONCILIATION_MANIFEST = (
    REVISION / "integrity/reconciled_analysis_manifest_20260711_110650.csv"
)
PEPTIDE_MANIFEST = (
    REVISION
    / "analysis_stats/final_peptide_denominator_manifest_20260711_131706.csv"
)
RAW_COUNT_MATRIX = (
    REVISION
    / "analysis_stats/reconstructed_raw_pfam_counts_20260711_131706.csv.gz"
)

PFAM_HMM_URL = (
    "https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00092?annotation=hmm"
)
PFAM_SEED_URL = (
    "https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00092?annotation=alignment:seed"
)
PDB_FASTA_URLS = {
    "1SHU": "https://www.rcsb.org/fasta/entry/1SHU/display",
    "1AO3": "https://www.rcsb.org/fasta/entry/1AO3/display",
}

EXPECTED_COHORT_SIZE = 126
EXPECTED_PHYLA = {"Rhodophyta": 70, "Ochrophyta": 43, "Chlorophyta": 13}
ORIGINAL_ACCESSION = "PF00092.31"
EXPECTED_CURRENT_ACCESSION = "PF00092.35"
EXPECTED_CURRENT_HMM_SHA256 = (
    "8489d4407aa1db1b0c94fa424c5f09fbd190670f967c0485b7e100b3dd936c64"
)
MIDAS_STATES = (7, 9, 11, 79, 112)
MIDAS_LABELS = ("D1", "S2", "S3", "T4", "D5")
MIDAS_EXPECTED = ("D", "S", "S", "T", "D")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, timeout: int = 60) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Macroalgae-PF00092-MIDAS-audit/20260715"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        payload = response.read()
    if not payload:
        raise RuntimeError(f"Empty response from {url}")
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help=(
            "Run the complete 126-cohort extraction, current-HMM search, and "
            "control checks in memory without creating result files."
        ),
    )
    parser.add_argument(
        "--cpus",
        type=int,
        default=4,
        help="Threads for pyhmmer.hmmsearch (default: 4).",
    )
    return parser.parse_args()


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def load_cohort() -> pd.DataFrame:
    manifest = pd.read_csv(RECONCILIATION_MANIFEST, low_memory=False)
    required = {
        "Genome",
        "Species",
        "Phylum",
        "pfam_source_path",
        "pfam_source_sha256",
        "safe_for_aef_pfam_analysis",
    }
    missing = required.difference(manifest.columns)
    if missing:
        raise ValueError(f"Reconciliation manifest lacks columns: {sorted(missing)}")

    cohort = manifest.loc[as_bool(manifest["safe_for_aef_pfam_analysis"])].copy()
    if len(cohort) != EXPECTED_COHORT_SIZE:
        raise ValueError(f"Expected 126 cohort rows; observed {len(cohort)}")
    if cohort["Genome"].nunique() != EXPECTED_COHORT_SIZE:
        raise ValueError("Cohort Genome IDs are not unique")
    observed_phyla = cohort["Phylum"].value_counts().to_dict()
    if observed_phyla != EXPECTED_PHYLA:
        raise ValueError(
            f"Unexpected phylum composition: {observed_phyla}; expected {EXPECTED_PHYLA}"
        )

    for path_text, expected_hash in cohort[
        ["pfam_source_path", "pfam_source_sha256"]
    ].itertuples(index=False, name=None):
        path = ROOT / str(path_text)
        if not path.is_file():
            raise FileNotFoundError(f"Missing retained HMM tblout: {path}")
        observed_hash = sha256_file(path)
        if observed_hash != expected_hash:
            raise ValueError(f"HMM tblout hash mismatch: {path}")
    return cohort


def parse_pf00092_tblout(path: Path) -> List[dict]:
    """Parse tagged or untagged HMMER tblout using accession-relative fields."""
    records: List[dict] = []
    with path.open("rt", encoding="utf-8", errors="strict") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields or "PF00092.31" not in fields:
                continue
            accession_index = fields.index("PF00092.31")
            if accession_index < 3 or len(fields) <= accession_index + 14:
                raise ValueError(f"Malformed PF00092 row at {path}:{line_number}")
            target_id = fields[accession_index - 3]
            records.append(
                {
                    "target_id": target_id,
                    "source_line": line_number,
                    "original_query_name": fields[accession_index - 1],
                    "original_accession": fields[accession_index],
                    "original_full_evalue": float(fields[accession_index + 1]),
                    "original_full_score": float(fields[accession_index + 2]),
                    "original_best_domain_evalue": float(fields[accession_index + 4]),
                    "original_best_domain_score": float(fields[accession_index + 5]),
                    "original_expected_domains": float(fields[accession_index + 7]),
                    "original_reported_domains": int(fields[accession_index + 13]),
                    "original_included_domains": int(fields[accession_index + 14]),
                }
            )
    targets = [record["target_id"] for record in records]
    if len(targets) != len(set(targets)):
        raise ValueError(f"Duplicate PF00092 target IDs in {path}")
    return records


def collect_original_hits(cohort: pd.DataFrame) -> pd.DataFrame:
    records: List[dict] = []
    for row in cohort.itertuples(index=False):
        source_path = ROOT / str(row.pfam_source_path)
        for hit in parse_pf00092_tblout(source_path):
            records.append(
                {
                    "Genome": row.Genome,
                    "Species": row.Species,
                    "Phylum": row.Phylum,
                    "tblout_path": str(source_path),
                    "tblout_sha256": row.pfam_source_sha256,
                    **hit,
                }
            )
    hits = pd.DataFrame.from_records(records)
    if hits.empty:
        raise ValueError("No PF00092.31 hit records were found")
    if hits.duplicated(["Genome", "target_id"]).any():
        raise ValueError("Genome/target PF00092 keys are not unique")

    raw_counts = pd.read_csv(
        RAW_COUNT_MATRIX, usecols=["Genome", "PF00092"], low_memory=False
    )
    if raw_counts["Genome"].duplicated().any():
        raise ValueError("Raw count matrix Genome IDs are not unique")
    if set(raw_counts["Genome"]) != set(cohort["Genome"]):
        raise ValueError("Raw count matrix Genome set differs from the 126-genome cohort")
    expected_total = int(raw_counts["PF00092"].sum())
    if len(hits) != expected_total:
        raise ValueError(
            f"Parsed PF00092 records ({len(hits)}) != raw matrix total ({expected_total})"
        )
    per_genome = hits.groupby("Genome").size()
    observed = raw_counts.set_index("Genome")["PF00092"].astype(int)
    parsed = per_genome.reindex(observed.index, fill_value=0).astype(int)
    if not parsed.equals(observed):
        mismatches = pd.DataFrame({"parsed": parsed, "matrix": observed})
        mismatches = mismatches[mismatches["parsed"] != mismatches["matrix"]]
        raise ValueError(f"Per-genome PF00092 mismatch:\n{mismatches}")
    return hits


def fasta_records(path: Path) -> Iterator[Tuple[str, str]]:
    name: str | None = None
    chunks: List[str] = []
    with path.open("rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"Sequence before first FASTA header in {path}")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def extract_hit_sequences(
    hits: pd.DataFrame, cohort: pd.DataFrame,
) -> Tuple[pd.DataFrame, Dict[str, str], dict]:
    peptide_manifest = pd.read_csv(PEPTIDE_MANIFEST, low_memory=False)
    if peptide_manifest["Genome"].duplicated().any():
        raise ValueError("Peptide manifest Genome IDs are not unique")
    expected_genomes = set(cohort["Genome"])
    if set(peptide_manifest["Genome"]) != expected_genomes:
        raise ValueError("Peptide manifest Genome set differs from the 126-genome cohort")
    peptide = peptide_manifest.set_index("Genome", drop=False)

    accounting = hits.copy()
    accounting["peptide_fasta"] = None
    accounting["peptide_fasta_sha256"] = None
    accounting["sequence_status"] = "not_processed"
    accounting["sequence_length"] = pd.NA
    accounting["sequence_sha256"] = None
    accounting["search_key"] = None

    sequences: Dict[str, str] = {}
    terminal_stops_removed = 0
    internal_stop_records: List[Tuple[str, str]] = []

    for genome, index_values in accounting.groupby("Genome").groups.items():
        genome_indices = list(index_values)
        peptide_row = peptide.loc[genome]
        fasta_text = peptide_row.get("peptide_fasta_realpath")
        expected_fasta_hash = peptide_row.get("peptide_fasta_sha256")
        if pd.isna(fasta_text) or not str(fasta_text).strip():
            accounting.loc[genome_indices, "sequence_status"] = (
                "unevaluable_missing_retained_peptide_fasta"
            )
            continue

        fasta_path = Path(str(fasta_text))
        if not fasta_path.is_file():
            raise FileNotFoundError(f"Peptide manifest path missing: {fasta_path}")
        observed_fasta_hash = sha256_file(fasta_path)
        if observed_fasta_hash != expected_fasta_hash:
            raise ValueError(f"Peptide FASTA hash mismatch: {fasta_path}")

        wanted_to_index = {
            accounting.at[index, "target_id"]: index for index in genome_indices
        }
        found: Dict[str, str] = {}
        for name, sequence in fasta_records(fasta_path):
            if name in wanted_to_index:
                if name in found:
                    raise ValueError(f"Duplicate FASTA header {name} in {fasta_path}")
                found[name] = sequence

        missing = sorted(set(wanted_to_index).difference(found))
        if missing:
            raise ValueError(
                f"{len(missing)} PF00092 target IDs absent from {fasta_path}: {missing[:5]}"
            )

        for target_id, sequence in found.items():
            index = wanted_to_index[target_id]
            raw_sequence = sequence.upper()
            cleaned = raw_sequence.rstrip("*")
            terminal_stops_removed += len(raw_sequence) - len(cleaned)
            if "*" in cleaned:
                accounting.at[index, "sequence_status"] = (
                    "unevaluable_internal_stop_character"
                )
                internal_stop_records.append((genome, target_id))
                continue
            search_key = f"hit_{index:06d}"
            sequences[search_key] = cleaned
            accounting.at[index, "peptide_fasta"] = str(fasta_path)
            accounting.at[index, "peptide_fasta_sha256"] = observed_fasta_hash
            accounting.at[index, "sequence_status"] = "extracted_for_current_hmm_search"
            accounting.at[index, "sequence_length"] = len(cleaned)
            accounting.at[index, "sequence_sha256"] = sha256_bytes(cleaned.encode("ascii"))
            accounting.at[index, "search_key"] = search_key

    diagnostics = {
        "terminal_stop_characters_removed": terminal_stops_removed,
        "internal_stop_records": internal_stop_records,
        "total_extracted_sequences": len(sequences),
        "total_extracted_residues": sum(map(len, sequences.values())),
    }
    return accounting, sequences, diagnostics


def load_hmm(payload_gzip: bytes) -> pyhmmer.plan7.HMM:
    try:
        payload = gzip.decompress(payload_gzip)
    except gzip.BadGzipFile as exc:
        raise ValueError("Official Pfam HMM response was not valid gzip") from exc
    with pyhmmer.plan7.HMMFile(io.BytesIO(payload)) as handle:
        hmm = handle.read()
        if handle.read() is not None:
            raise ValueError("Expected exactly one HMM in PF00092 payload")
    if hmm is None:
        raise ValueError("PF00092 HMM payload contained no model")
    accession = hmm.accession.decode() if hmm.accession else ""
    if accession != EXPECTED_CURRENT_ACCESSION:
        raise ValueError(
            f"Expected pinned HMM {EXPECTED_CURRENT_ACCESSION}; downloaded {accession}"
        )
    if hmm.M < max(MIDAS_STATES):
        raise ValueError(f"PF00092 model length {hmm.M} does not cover MIDAS states")
    return hmm


def parse_rcsb_fasta(payload: bytes, pdb_id: str) -> str:
    text = payload.decode("utf-8")
    sequences = [sequence for _, sequence in _fasta_text_records(text)]
    if not sequences:
        raise ValueError(f"No sequence in RCSB FASTA for {pdb_id}")
    unique = list(dict.fromkeys(sequences))
    if len(unique) != 1:
        raise ValueError(f"Expected one unique sequence for {pdb_id}; got {len(unique)}")
    return unique[0].upper().rstrip("*")


def _fasta_text_records(text: str) -> Iterator[Tuple[str, str]]:
    name: str | None = None
    chunks: List[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        else:
            if name is None:
                raise ValueError("Malformed FASTA payload")
            chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def aligned_state_residues(
    alignment: pyhmmer.plan7.Alignment,
    states: Sequence[int] = MIDAS_STATES,
) -> Dict[int, str | None]:
    state = int(alignment.hmm_from)
    residues: Dict[int, str | None] = {position: None for position in states}
    for hmm_character, target_character in zip(
        alignment.hmm_sequence, alignment.target_sequence
    ):
        if hmm_character in ".-":
            continue
        if state in residues:
            residue = target_character.upper()
            residues[state] = None if residue in {"-", "."} else residue
        state += 1
    return residues


def classify_residues(residue_map: Mapping[int, str | None]) -> Tuple[str, str, int]:
    observed = tuple(residue_map[position] for position in MIDAS_STATES)
    observed_text = "".join(residue if residue is not None else "-" for residue in observed)
    matches = sum(
        residue == expected for residue, expected in zip(observed, MIDAS_EXPECTED)
    )
    if observed == MIDAS_EXPECTED:
        classification = "complete_canonical_MIDAS_signature"
    elif any(residue is None for residue in observed):
        classification = "incomplete_alignment_at_MIDAS_states"
    elif observed[:3] == MIDAS_EXPECTED[:3]:
        classification = "noncanonical_complete_region1_DxSxS_retained"
    else:
        classification = "noncanonical_complete_signature"
    return classification, observed_text, matches


def search_controls(
    hmm: pyhmmer.plan7.HMM, control_sequences: Mapping[str, str]
) -> dict:
    alphabet = pyhmmer.easel.Alphabet.amino()
    digital = [
        pyhmmer.easel.TextSequence(name=name.encode(), sequence=sequence).digitize(
            alphabet
        )
        for name, sequence in control_sequences.items()
    ]
    top_hits = next(
        pyhmmer.hmmsearch(hmm, digital, cpus=1, bit_cutoffs="gathering")
    )
    results = {}
    for hit in top_hits:
        name = hit.name.decode()
        included_domains = [domain for domain in hit.domains if domain.included]
        if not hit.included or not included_domains:
            continue
        domain = max(included_domains, key=lambda item: item.score)
        residues = aligned_state_residues(domain.alignment)
        classification, observed, matches = classify_residues(residues)
        results[name] = {
            "observed_states": observed,
            "classification": classification,
            "matching_states": matches,
            "domain_score": float(domain.score),
            "target_from": int(domain.alignment.target_from),
            "target_to": int(domain.alignment.target_to),
        }
    if results.get("1SHU", {}).get("observed_states") != "DSSTD":
        raise ValueError(f"Positive structural control failed: {results.get('1SHU')}")
    if results.get("1AO3", {}).get("observed_states") != "DSSST":
        raise ValueError(f"Vestigial structural control failed: {results.get('1AO3')}")
    return results


def search_macroalgal_sequences(
    hmm: pyhmmer.plan7.HMM,
    sequences: Mapping[str, str],
    accounting: pd.DataFrame,
    cpus: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    alphabet = pyhmmer.easel.Alphabet.amino()
    digital_sequences = []
    invalid: Dict[str, str] = {}
    for key, sequence in sequences.items():
        try:
            digital_sequences.append(
                pyhmmer.easel.TextSequence(
                    name=key.encode("ascii"), sequence=sequence
                ).digitize(alphabet)
            )
        except (ValueError, UnicodeEncodeError) as exc:
            invalid[key] = str(exc)

    if invalid:
        for key, reason in invalid.items():
            index = int(key.rsplit("_", 1)[1])
            accounting.at[index, "sequence_status"] = (
                f"unevaluable_pyhmmer_sequence_error:{reason}"
            )
        digital_sequences = [
            sequence
            for sequence in digital_sequences
            if sequence.name.decode("ascii") not in invalid
        ]

    top_hits = next(
        pyhmmer.hmmsearch(
            hmm,
            digital_sequences,
            cpus=max(1, cpus),
            bit_cutoffs="gathering",
        )
    )

    result_rows: List[dict] = []
    included_keys: set[str] = set()
    for hit in top_hits:
        key = hit.name.decode("ascii")
        if not hit.included:
            continue
        included_domains = [domain for domain in hit.domains if domain.included]
        if not included_domains:
            continue
        included_keys.add(key)
        index = int(key.rsplit("_", 1)[1])
        source = accounting.loc[index]
        for domain_number, domain in enumerate(included_domains, start=1):
            alignment = domain.alignment
            residues = aligned_state_residues(alignment)
            classification, observed, matches = classify_residues(residues)
            aligned_target = alignment.target_sequence
            domain_sequence = "".join(
                character.upper()
                for character in aligned_target
                if character not in {"-", "."}
            )
            row = {
                "Genome": source["Genome"],
                "Species": source["Species"],
                "Phylum": source["Phylum"],
                "target_id": source["target_id"],
                "sequence_length": int(source["sequence_length"]),
                "sequence_sha256": source["sequence_sha256"],
                "original_accession": source["original_accession"],
                "original_full_score": source["original_full_score"],
                "original_best_domain_score": source["original_best_domain_score"],
                "original_included_domains": source["original_included_domains"],
                "current_hmm_accession": hmm.accession.decode(),
                "current_sequence_score": float(hit.score),
                "current_domain_number": domain_number,
                "current_domain_score": float(domain.score),
                "current_domain_i_evalue": float(domain.i_evalue),
                "hmm_from": int(alignment.hmm_from),
                "hmm_to": int(alignment.hmm_to),
                "target_from": int(alignment.target_from),
                "target_to": int(alignment.target_to),
                "D1_hmm_state_7": residues[7],
                "S2_hmm_state_9": residues[9],
                "S3_hmm_state_11": residues[11],
                "T4_hmm_state_79": residues[79],
                "D5_hmm_state_112": residues[112],
                "observed_MIDAS_states": observed,
                "matching_canonical_states": matches,
                "MIDAS_sequence_classification": classification,
                "complete_canonical_MIDAS_signature": (
                    classification == "complete_canonical_MIDAS_signature"
                ),
                "aligned_domain_sequence": domain_sequence,
                "tblout_path": source["tblout_path"],
                "peptide_fasta": source["peptide_fasta"],
            }
            result_rows.append(row)

    extracted_mask = accounting["sequence_status"].eq(
        "extracted_for_current_hmm_search"
    )
    for index in accounting.index[extracted_mask]:
        key = accounting.at[index, "search_key"]
        accounting.at[index, "sequence_status"] = (
            "current_PF00092_GA_hit"
            if key in included_keys
            else "not_retained_by_current_PF00092_GA"
        )

    domains = pd.DataFrame.from_records(result_rows)
    if domains.empty:
        raise ValueError("Current PF00092 HMM produced no included domains")
    return domains, accounting


def summarize(
    cohort: pd.DataFrame,
    accounting: pd.DataFrame,
    domains: pd.DataFrame,
    extraction_diagnostics: dict,
    controls: dict,
    elapsed_seconds: float,
) -> dict:
    class_counts = domains["MIDAS_sequence_classification"].value_counts().to_dict()
    status_counts = accounting["sequence_status"].value_counts().to_dict()
    canonical = domains["complete_canonical_MIDAS_signature"].astype(bool)
    summary = {
        "cohort_genomes": int(cohort["Genome"].nunique()),
        "cohort_phyla": {
            key: int(value) for key, value in cohort["Phylum"].value_counts().items()
        },
        "original_pf00092_protein_hit_records": int(len(accounting)),
        "original_pf00092_positive_genomes": int(accounting["Genome"].nunique()),
        "original_tblout_included_domain_estimate": int(
            accounting["original_included_domains"].sum()
        ),
        "sequence_status_counts": {key: int(value) for key, value in status_counts.items()},
        "extracted_sequences": int(extraction_diagnostics["total_extracted_sequences"]),
        "extracted_residues": int(extraction_diagnostics["total_extracted_residues"]),
        "current_hmm_included_proteins": int(
            accounting["sequence_status"].eq("current_PF00092_GA_hit").sum()
        ),
        "current_hmm_included_domains": int(len(domains)),
        "domain_classification_counts": {
            key: int(value) for key, value in class_counts.items()
        },
        "fully_evaluable_domains": int(
            len(domains)
            - class_counts.get("incomplete_alignment_at_MIDAS_states", 0)
        ),
        "fully_evaluable_noncanonical_domains": int(
            class_counts.get("noncanonical_complete_region1_DxSxS_retained", 0)
            + class_counts.get("noncanonical_complete_signature", 0)
        ),
        "incomplete_alignment_domains": int(
            class_counts.get("incomplete_alignment_at_MIDAS_states", 0)
        ),
        "canonical_signature_domains": int(canonical.sum()),
        "canonical_signature_proteins": int(
            domains.loc[canonical, ["Genome", "target_id"]].drop_duplicates().shape[0]
        ),
        "canonical_signature_genomes": int(
            domains.loc[canonical, "Genome"].nunique()
        ),
        "structural_controls": controls,
        "extraction_diagnostics": {
            "terminal_stop_characters_removed": int(
                extraction_diagnostics["terminal_stop_characters_removed"]
            ),
            "internal_stop_record_count": int(
                len(extraction_diagnostics["internal_stop_records"])
            ),
        },
        "analysis_elapsed_seconds": elapsed_seconds,
    }
    return summary


def build_genome_summary(
    cohort: pd.DataFrame, accounting: pd.DataFrame, domains: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    for row in cohort.itertuples(index=False):
        acc = accounting[accounting["Genome"] == row.Genome]
        dom = domains[domains["Genome"] == row.Genome]
        canonical = dom["complete_canonical_MIDAS_signature"].astype(bool)
        rows.append(
            {
                "Genome": row.Genome,
                "Species": row.Species,
                "Phylum": row.Phylum,
                "original_PF00092_protein_hits": int(len(acc)),
                "original_tblout_included_domain_estimate": int(
                    acc["original_included_domains"].sum()
                ),
                "extracted_proteins": int(
                    acc["sequence_sha256"].notna().sum()
                ),
                "current_hmm_included_proteins": int(
                    acc["sequence_status"].eq("current_PF00092_GA_hit").sum()
                ),
                "current_hmm_included_domains": int(len(dom)),
                "canonical_signature_domains": int(canonical.sum()),
                "canonical_signature_proteins": int(
                    dom.loc[canonical, "target_id"].nunique()
                ),
                "fully_evaluable_noncanonical_domains": int(
                    dom["MIDAS_sequence_classification"]
                    .isin(
                        {
                            "noncanonical_complete_region1_DxSxS_retained",
                            "noncanonical_complete_signature",
                        }
                    )
                    .sum()
                ),
                "incomplete_alignment_domains": int(
                    dom["MIDAS_sequence_classification"]
                    .eq("incomplete_alignment_at_MIDAS_states")
                    .sum()
                ),
            }
        )
    return pd.DataFrame(rows)


def write_csv(path: Path, frame: pd.DataFrame) -> None:
    if path.exists():
        raise FileExistsError(path)
    frame.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)


def output_summary_text(summary: dict, hmm_accession: str) -> str:
    classes = summary["domain_classification_counts"]
    statuses = summary["sequence_status_counts"]
    return "\n".join(
        [
            "PF00092 MIDAS SEQUENCE AUDIT",
            "=" * 31,
            "",
            "Exact computational rule:",
            f"  A complete canonical signature requires D/S/S/T/D at {EXPECTED_CURRENT_ACCESSION}",
            "  HMM match states 7/9/11/79/112, respectively.",
            "  States are read from current-Pfam domain alignments passing the",
            "  model-specific gathering threshold.",
            "",
            f"Cohort: {summary['cohort_genomes']} genomes",
            f"Original PF00092.31 protein-hit records: {summary['original_pf00092_protein_hit_records']}",
            f"Original tblout included-domain estimate: {summary['original_tblout_included_domain_estimate']}",
            f"Extracted exact retained proteins: {summary['extracted_sequences']}",
            f"Current model: {hmm_accession}",
            f"Current-model included proteins: {summary['current_hmm_included_proteins']}",
            f"Current-model included domains: {summary['current_hmm_included_domains']}",
            f"Complete canonical signature domains: {summary['canonical_signature_domains']}",
            f"Fully evaluable noncanonical domains: {summary['fully_evaluable_noncanonical_domains']}",
            f"Incomplete alignments at >=1 MIDAS state: {summary['incomplete_alignment_domains']}",
            f"Complete canonical signature proteins: {summary['canonical_signature_proteins']}",
            f"Genomes with >=1 complete canonical signature: {summary['canonical_signature_genomes']}",
            "",
            "Domain classifications:",
            *[f"  {key}: {value}" for key, value in sorted(classes.items())],
            "",
            "Sequence accounting:",
            *[f"  {key}: {value}" for key, value in sorted(statuses.items())],
            "",
            "Structural alignment controls:",
            f"  1SHU: {summary['structural_controls']['1SHU']['observed_states']}",
            f"  1AO3: {summary['structural_controls']['1AO3']['observed_states']}",
            "",
            "Classification scope:",
            "  'Complete canonical signature' is a sequence-level classification at",
            "  five topologically homologous positions. Fully evaluable noncanonical",
            "  domains lack the complete five-residue signature; incomplete alignments",
            "  remain unevaluable.",
        ]
    )


def run(args: argparse.Namespace) -> int:
    started = time.perf_counter()
    cohort = load_cohort()
    hits = collect_original_hits(cohort)
    accounting, sequences, extraction_diagnostics = extract_hit_sequences(hits, cohort)

    hmm_payload = fetch(PFAM_HMM_URL)
    hmm_payload_sha256 = sha256_bytes(hmm_payload)
    if hmm_payload_sha256 != EXPECTED_CURRENT_HMM_SHA256:
        raise ValueError(
            "Downloaded PF00092 HMM checksum differs from the pinned PF00092.35 "
            f"model: observed {hmm_payload_sha256}"
        )
    seed_payload = fetch(PFAM_SEED_URL)
    pdb_payloads = {name: fetch(url) for name, url in PDB_FASTA_URLS.items()}
    hmm = load_hmm(hmm_payload)
    controls = search_controls(
        hmm,
        {name: parse_rcsb_fasta(payload, name) for name, payload in pdb_payloads.items()},
    )
    domains, accounting = search_macroalgal_sequences(
        hmm, sequences, accounting, cpus=args.cpus
    )
    elapsed = time.perf_counter() - started
    summary = summarize(
        cohort,
        accounting,
        domains,
        extraction_diagnostics,
        controls,
        elapsed,
    )
    hmm_accession = hmm.accession.decode()

    print(output_summary_text(summary, hmm_accession))
    print(f"\nMeasured complete-run time: {elapsed:.2f} seconds")
    approx_sequence_bytes = extraction_diagnostics["total_extracted_residues"]
    print(
        "Extracted sequence characters processed: "
        f"{approx_sequence_bytes:,} (~{approx_sequence_bytes / 1024**2:.2f} MiB raw text)"
    )

    if args.validate_only:
        print("VALIDATION PASSED: full cohort and full sequence set; no files created.")
        return 0

    timestamp = datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")
    output_dir = SCRIPT_PATH.parent / f"PF00092_MIDAS_run_{timestamp}"
    output_dir.mkdir(parents=False, exist_ok=False)

    paths = {
        "domain_calls": output_dir / f"PF00092_MIDAS_domain_calls_{timestamp}.csv",
        "sequence_accounting": output_dir
        / f"PF00092_MIDAS_sequence_accounting_{timestamp}.csv",
        "genome_summary": output_dir / f"PF00092_MIDAS_genome_summary_{timestamp}.csv",
        "canonical_fasta": output_dir
        / f"PF00092_complete_canonical_domains_{timestamp}.fasta",
        "summary": output_dir / f"PF00092_MIDAS_summary_{timestamp}.txt",
        "pfam_hmm": output_dir / f"PF00092_official_HMM_{timestamp}.hmm.gz",
        "pfam_seed": output_dir / f"PF00092_official_seed_{timestamp}.sto.gz",
        "pdb_1shu": output_dir / f"PDB_1SHU_control_{timestamp}.fasta",
        "pdb_1ao3": output_dir / f"PDB_1AO3_control_{timestamp}.fasta",
        "provenance": output_dir / f"PF00092_MIDAS_provenance_{timestamp}.json",
    }

    write_csv(paths["domain_calls"], domains)
    write_csv(paths["sequence_accounting"], accounting)
    write_csv(paths["genome_summary"], build_genome_summary(cohort, accounting, domains))

    canonical_domains = domains[
        domains["complete_canonical_MIDAS_signature"].astype(bool)
    ]
    with paths["canonical_fasta"].open("x", encoding="utf-8") as handle:
        for row in canonical_domains.itertuples(index=False):
            header = (
                f">{row.Genome}|{row.target_id}|domain{row.current_domain_number}|"
                f"target_{row.target_from}_{row.target_to}|{hmm_accession}"
            )
            handle.write(header + "\n")
            sequence = row.aligned_domain_sequence
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    paths["summary"].write_text(
        output_summary_text(summary, hmm_accession) + "\n", encoding="utf-8"
    )
    paths["pfam_hmm"].write_bytes(hmm_payload)
    paths["pfam_seed"].write_bytes(seed_payload)
    paths["pdb_1shu"].write_bytes(pdb_payloads["1SHU"])
    paths["pdb_1ao3"].write_bytes(pdb_payloads["1AO3"])

    provenance = {
        "created_at": datetime.now().astimezone().isoformat(),
        "command_parameters": {
            "cpus": int(args.cpus),
            "validate_only": bool(args.validate_only),
        },
        "script": {
            "path": str(SCRIPT_PATH),
            "sha256": sha256_file(SCRIPT_PATH),
        },
        "inputs": {
            "reconciliation_manifest": {
                "path": str(RECONCILIATION_MANIFEST),
                "sha256": sha256_file(RECONCILIATION_MANIFEST),
            },
            "peptide_manifest": {
                "path": str(PEPTIDE_MANIFEST),
                "sha256": sha256_file(PEPTIDE_MANIFEST),
            },
            "raw_count_matrix": {
                "path": str(RAW_COUNT_MATRIX),
                "sha256": sha256_file(RAW_COUNT_MATRIX),
            },
            "pfam_hmm": {
                "url": PFAM_HMM_URL,
                "download_sha256": hmm_payload_sha256,
                "expected_sha256": EXPECTED_CURRENT_HMM_SHA256,
                "saved_path": str(paths["pfam_hmm"]),
                "accession": hmm_accession,
                "model_length": int(hmm.M),
                "gathering_cutoffs": [float(value) for value in hmm.cutoffs.gathering],
            },
            "pfam_seed": {
                "url": PFAM_SEED_URL,
                "download_sha256": sha256_bytes(seed_payload),
                "saved_path": str(paths["pfam_seed"]),
            },
            "structural_controls": {
                name: {
                    "url": PDB_FASTA_URLS[name],
                    "download_sha256": sha256_bytes(payload),
                    "saved_path": str(paths[f"pdb_{name.lower()}"]),
                }
                for name, payload in pdb_payloads.items()
            },
        },
        "method": {
            "original_hit_definition": "PF00092.31 tblout protein-hit row",
            "fresh_reannotation_model": hmm_accession,
            "fresh_reannotation_threshold": "model-specific Pfam gathering cutoff",
            "midas_hmm_states": list(MIDAS_STATES),
            "midas_labels": list(MIDAS_LABELS),
            "canonical_expected_residues": list(MIDAS_EXPECTED),
            "canonical_rule": (
                "Complete canonical MIDAS signature iff residues aligned to "
                f"{EXPECTED_CURRENT_ACCESSION} HMM states 7,9,11,79,112 equal D,S,S,T,D"
            ),
            "classification_scope": (
                "Sequence-level classification at five homologous HMM states; "
                "incomplete alignments remain unevaluable."
            ),
        },
        "summary": summary,
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "pandas": pd.__version__,
            "pyhmmer": pyhmmer.__version__,
        },
        "outputs": {},
    }
    for key, path in paths.items():
        if key == "provenance":
            continue
        provenance["outputs"][key] = {
            "path": str(path),
            "sha256": sha256_file(path),
            "bytes": path.stat().st_size,
        }
    with paths["provenance"].open("x", encoding="utf-8") as handle:
        json.dump(provenance, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"\nOutputs: {output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run(parse_args()))

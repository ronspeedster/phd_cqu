# %% Cell 1 - Imports and configuration
from __future__ import annotations

from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import product
import os
from pathlib import Path
import argparse
import csv
from decimal import Decimal, InvalidOperation
import logging
import re
import shutil
import subprocess
import sys


LOGGER = logging.getLogger("hypermute3_pairwise")
DNA_STRICT_ALLOWED = set("ACGT-")

LANL_HEADER_FIELDS: tuple[str, ...] = (
    "Se ID", "Patient Code", "PAT id(SSAM)", "Accession", "Name",
    "Subtype", "Country", "Sampling Year", "Problematic Sequence",
    "HXB2/MAC239 start", "HXB2/MAC239 stop", "Sequence Length", "Organism",
)


@dataclass
class HypermutRunConfig:
    name: str
    match: str  # strict | partial
    keepgaps: bool


@dataclass
class PipelineConfig:
    group_manifest_csv: Path
    sequence_to_group_csv: Path
    hypermut_script: Path
    output_dir: Path
    python_executable: str = sys.executable
    mafft_binary: str = "mafft"
    mafft_threads: int = -1
    mutation_from: str = "G"
    mutation_to: str = "A"
    upstream_context: str = ""
    downstream_context: str = "RD"
    enforce: str = "D"
    begin: int = 0
    finish: int | None = None
    number_of_sequence: int | None = None
    run_all_mutation_directions: bool = False
    reuse_existing_pairwise_alignments: bool = False
    parallel_workers: int | None = None
    hxb2_reference_fasta: Path | None = None
    log_level: str = "INFO"
    full_run: bool = False
    problematic_only: bool = False


def configure_logging(log_level: str) -> None:
    level_name = (log_level or "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def get_mutation_directions(cfg: PipelineConfig) -> list[tuple[str, str]]:
    if cfg.run_all_mutation_directions:
        bases = ["A", "C", "G", "T"]
        return [(base_from, base_to) for base_from in bases for base_to in bases if base_from != base_to]
    return [(cfg.mutation_from, cfg.mutation_to)]


def ensure_hxb2_reference(data_dir: Path) -> Path:
    hxb2_path = data_dir / "hxb2.fasta"
    if hxb2_path.exists():
        LOGGER.info("HXB2 reference already exists: %s", hxb2_path)
        return hxb2_path
    LOGGER.info("Downloading HXB2 reference from LANL...")
    import urllib.request
    url = (
        "https://www.hiv.lanl.gov/cgi-bin/NEWALIGN/saveAlign.cgi"
        "?format=fasta&alignID=HXB2&bg=genomicDNA"
    )
    try:
        urllib.request.urlretrieve(url, hxb2_path)
        LOGGER.info("HXB2 reference saved: %s", hxb2_path)
    except Exception as exc:
        raise RuntimeError(
            f"Failed to download HXB2 reference: {exc}. "
            f"Please download manually and save to {hxb2_path}"
        ) from exc
    return hxb2_path


# %% Cell 2 - FASTA and formatting helpers
def read_fasta_records(file_path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = None
    seq_parts: list[str] = []

    with file_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        records.append((header, "".join(seq_parts)))

    return records


def write_fasta_records(
    records: list[tuple[str, str]],
    file_path: Path,
) -> None:
    with file_path.open("w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for idx in range(0, len(seq), 80):
                handle.write(seq[idx : idx + 80] + "\n")


def _parse_lanl_header(header: str) -> dict[str, str]:
    parts = header.split(",")
    parsed = {}
    for i, field_name in enumerate(LANL_HEADER_FIELDS):
        parsed[field_name] = parts[i].strip() if i < len(parts) else ""
    return parsed


def _is_problematic(seq_header: str) -> bool:
    """Return True if the sequence is marked as problematic (value 3)."""
    fields = _parse_lanl_header(seq_header)
    return fields.get("Problematic Sequence", "").strip() == "3"


MONOMER_MUTATION_COLUMNS: list[str] = (
    ["consensus_A_count", "consensus_G_count", "consensus_C_count", "consensus_T_count"]
    + [f"{a}_to_{b}" for a in "ACGT" for b in "ACGT" if a != b]
)


def generate_nmer_mutation_columns(n: int) -> list[str]:
    bases = "ACGT"
    nmers = ["".join(combo) for combo in product(bases, repeat=n)]
    columns = []
    for src in nmers:
        for dst in nmers:
            if src != dst:
                columns.append(f"{src}_to_{dst}")
    return columns


def count_monomer_mutations(consensus: str, query: str) -> dict[str, int]:
    counts: dict[str, int] = {col: 0 for col in MONOMER_MUTATION_COLUMNS}
    for c_base, q_base in zip(consensus.upper(), query.upper()):
        if c_base == "-" and q_base == "-":
            continue
        if c_base in "ACGT":
            counts[f"consensus_{c_base}_count"] += 1
            if q_base in "ACGT" and c_base != q_base:
                counts[f"{c_base}_to_{q_base}"] += 1
    return counts


def count_nmer_mutations(consensus: str, query: str, n: int) -> dict[str, int]:
    columns = generate_nmer_mutation_columns(n)
    counts: dict[str, int] = {col: 0 for col in columns}
    con_upper = consensus.upper()
    que_upper = query.upper()
    for i in range(len(con_upper) - n + 1):
        c_kmer = con_upper[i : i + n]
        q_kmer = que_upper[i : i + n]
        if "-" in c_kmer or "-" in q_kmer:
            continue
        if all(b in "ACGT" for b in c_kmer) and all(b in "ACGT" for b in q_kmer):
            if c_kmer != q_kmer:
                key = f"{c_kmer}_to_{q_kmer}"
                if key in counts:
                    counts[key] += 1
    return counts


def _strip_flanking_gaps(consensus: str, query: str) -> tuple[str, str]:
    start = 0
    while start < len(consensus) and (consensus[start] == "-" or query[start] == "-"):
        start += 1
    end = len(consensus)
    while end > start and (consensus[end - 1] == "-" or query[end - 1] == "-"):
        end -= 1
    return consensus[start:end], query[start:end]


def strip_pair_fasta_flanks(fasta_path: Path) -> None:
    records = read_fasta_records(fasta_path)
    if len(records) < 2:
        return
    con_trimmed, que_trimmed = _strip_flanking_gaps(records[0][1], records[1][1])
    write_fasta_records(
        [(records[0][0], con_trimmed), (records[1][0], que_trimmed)],
        fasta_path,
    )


@dataclass
class AlignmentQC:
    aligned_positions: int
    total_positions: int
    gap_fraction_consensus: float
    gap_fraction_query: float
    identity: float
    passed: bool


def validate_alignment_quality(
    consensus: str,
    query: str,
    min_identity: float = 0.3,
    max_gap_fraction: float = 0.5,
) -> AlignmentQC:
    total = len(consensus)
    if total == 0:
        return AlignmentQC(0, 0, 1.0, 1.0, 0.0, False)

    con_upper = consensus.upper()
    que_upper = query.upper()

    gaps_c = sum(1 for b in con_upper if b == "-")
    gaps_q = sum(1 for b in que_upper if b == "-")
    gap_frac_c = gaps_c / total
    gap_frac_q = gaps_q / total

    aligned = 0
    matches = 0
    for c_base, q_base in zip(con_upper, que_upper):
        if c_base == "-" or q_base == "-":
            continue
        aligned += 1
        if c_base == q_base:
            matches += 1

    identity = matches / aligned if aligned > 0 else 0.0
    passed = identity >= min_identity and gap_frac_q <= max_gap_fraction

    return AlignmentQC(
        aligned_positions=aligned,
        total_positions=total,
        gap_fraction_consensus=gap_frac_c,
        gap_fraction_query=gap_frac_q,
        identity=identity,
        passed=passed,
    )


def normalize_sequence(seq: str) -> str:
    seq = seq.upper().replace("U", "T")
    return "".join(base if base in DNA_STRICT_ALLOWED else "-" for base in seq)


def to_non_scientific(value: str) -> str:
    text = (value or "").strip()
    if not text:
        return text
    if "e" not in text.lower():
        return text
    try:
        return format(Decimal(text), "f")
    except (InvalidOperation, ValueError):
        return text


def sanitize_file_component(text: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    safe = safe.strip("_")
    return safe[:90] if safe else "seq"


def resolve_manifest_path(path_text: str, group_manifest_csv: Path) -> Path:
    """Resolve FASTA paths from manifest across Windows/macOS environments."""
    text = (path_text or "").strip()
    if not text:
        return Path("")

    # Normalize slashes first so we can parse path segments consistently.
    normalized = text.replace("\\", "/")

    # 1) Direct path as written in manifest.
    direct = Path(normalized)
    if direct.exists():
        return direct

    # 2) Relative to project root (script lives in project root).
    project_root = group_manifest_csv.resolve().parents[2]
    relative_candidate = project_root / normalized
    if relative_candidate.exists():
        return relative_candidate

    # 3) Convert absolute Windows paths containing /data/... to local project path.
    lower = normalized.lower()
    marker = "/data/"
    idx = lower.find(marker)
    if idx != -1:
        data_tail = normalized[idx + len(marker):]
        data_candidate = project_root / "data" / data_tail
        if data_candidate.exists():
            return data_candidate

    return direct


# %% Cell 3 - Pairwise alignment and Hypermut execution
def align_pair_with_mafft(
    consensus_header: str,
    consensus_seq: str,
    query_header: str,
    query_seq: str,
    cfg: PipelineConfig,
    pair_input_path: Path,
    pair_aligned_path: Path,
) -> list[tuple[str, str]]:
    mafft_exe = shutil.which(cfg.mafft_binary)
    if not mafft_exe:
        raise RuntimeError(
            "MAFFT was not found on PATH. Install MAFFT to run pairwise alignment."
        )

    input_records = [
        (consensus_header, normalize_sequence(consensus_seq)),
        (query_header, normalize_sequence(query_seq)),
    ]
    write_fasta_records(input_records, pair_input_path)

    cmd = [
        mafft_exe,
        "--thread",
        str(cfg.mafft_threads),
        "--auto",
        "--inputorder",
        str(pair_input_path),
    ]
    proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"MAFFT failed ({proc.returncode}) for query '{query_header}':\n"
            f"{proc.stderr}"
        )

    pair_aligned_path.write_text(proc.stdout, encoding="utf-8")
    aligned_records = read_fasta_records(pair_aligned_path)
    if len(aligned_records) < 2:
        raise RuntimeError(
            f"Aligned FASTA malformed for query '{query_header}': "
            f"expected at least 2 records, got {len(aligned_records)}"
        )

    # Keep the first two records in the expected order for Hypermut3.
    return [
        (aligned_records[0][0], normalize_sequence(aligned_records[0][1])),
        (aligned_records[1][0], normalize_sequence(aligned_records[1][1])),
    ]



def run_hypermut_once(
    cfg: PipelineConfig,
    run_cfg: HypermutRunConfig,
    aligned_pair_fasta: Path,
    seq_output_dir: Path,
    mutation_from: str,
    mutation_to: str,
    direction_tag: str,
) -> dict[str, str]:
    prefix = str(seq_output_dir / f"{run_cfg.name}-{direction_tag}-")
    run_input_fasta = aligned_pair_fasta

    cmd = [
        cfg.python_executable,
        str(cfg.hypermut_script),
        str(run_input_fasta),
        mutation_from,
        mutation_to,
        "-m",
        run_cfg.match,
        "-e",
        cfg.enforce,
        "-p",
        prefix,
    ]

    if cfg.upstream_context:
        cmd.extend(["-u", cfg.upstream_context])
    if cfg.downstream_context:
        cmd.extend(["-d", cfg.downstream_context])
    if run_cfg.keepgaps:
        cmd.append("-k")
    if cfg.begin > 0:
        cmd.extend(["-b", str(cfg.begin)])
    if cfg.finish is not None:
        cmd.extend(["-f", str(cfg.finish)])

    proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            "Hypermut3 failed for mode "
            f"'{run_cfg.name}' (code {proc.returncode})\n"
            f"STDERR:\n{proc.stderr}\nSTDOUT:\n{proc.stdout}"
        )

    summary_csv = seq_output_dir / f"{run_cfg.name}-{direction_tag}-summary.csv"
    positions_csv = seq_output_dir / f"{run_cfg.name}-{direction_tag}-positions.csv"
    args_csv = seq_output_dir / f"{run_cfg.name}-{direction_tag}-args.csv"

    return {
        "run_name": run_cfg.name,
        "match": run_cfg.match,
        "keepgaps": str(run_cfg.keepgaps),
        "mutation_from": mutation_from,
        "mutation_to": mutation_to,
        "input_fasta": str(run_input_fasta),
        "summary_csv": str(summary_csv),
        "positions_csv": str(positions_csv),
        "args_csv": str(args_csv),
    }


def read_single_hypermut_summary(summary_csv: Path) -> dict[str, str]:
    empty = {
        "seq_name": "",
        "primary_matches": "",
        "potential_primaries": "",
        "control_matches": "",
        "potential_controls": "",
        "rate_ratio": "",
        "fisher_p": "",
    }

    if not summary_csv.exists():
        return empty

    # Hypermut3 may write unquoted seq_name values containing commas.
    # Parse from the right: last 6 columns are numeric metrics, the left side is seq_name.
    lines = summary_csv.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        return empty

    parts = [p.strip() for p in lines[1].split(",")]
    if len(parts) < 7:
        return empty

    metrics = parts[-6:]
    seq_name = ",".join(parts[:-6]).strip()
    return {
        "seq_name": seq_name,
        "primary_matches": to_non_scientific(metrics[0]),
        "potential_primaries": to_non_scientific(metrics[1]),
        "control_matches": to_non_scientific(metrics[2]),
        "potential_controls": to_non_scientific(metrics[3]),
        "rate_ratio": to_non_scientific(metrics[4]),
        "fisher_p": to_non_scientific(metrics[5]),
    }


# %% Cell 4 - Group mapping and CSV helpers
def load_group_consensus_map(group_manifest_csv: Path) -> dict[str, tuple[str, str]]:
    if not group_manifest_csv.exists():
        raise FileNotFoundError(f"Group manifest not found: {group_manifest_csv}")

    out: dict[str, tuple[str, str]] = {}
    with group_manifest_csv.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            full_key = row.get("full_key", "").strip()
            sampling_year = row.get("sampling_year", "").strip() or "UNKNOWN"
            hxb2_start = row.get("hxb2_start", "").strip() or "UNKNOWN"
            cluster_id = row.get("cluster_id", "").strip()
            bin_id = row.get("bin_id", "").strip()
            reduced_key = f"SamplingYear={sampling_year}|HXB2Start={hxb2_start}"
            consensus_fasta = resolve_manifest_path(
                row.get("consensus_fasta", ""),
                group_manifest_csv,
            )
            if not consensus_fasta.exists():
                continue
            records = read_fasta_records(consensus_fasta)
            if not records:
                continue

            # Legacy key used by older sequence_to_group_mapping.csv outputs.
            out[reduced_key] = (records[0][0], records[0][1])

            # New clustering key used by sampling_year_clustering pipeline.
            if cluster_id and bin_id:
                clustered_key = (
                    f"SamplingYear={sampling_year}|ClusterId={cluster_id}|BinId={bin_id}"
                )
                out[clustered_key] = (records[0][0], records[0][1])

            # Keep compatibility with full-key lookups when available.
            if full_key:
                out[full_key] = (records[0][0], records[0][1])
                if "|n=" in full_key:
                    out[full_key.split("|n=", 1)[0]] = (records[0][0], records[0][1])

    return out


def load_sequence_to_group_map(sequence_to_group_csv: Path) -> dict[str, str]:
    if not sequence_to_group_csv.exists():
        raise FileNotFoundError(
            f"Sequence-to-group CSV not found: {sequence_to_group_csv}"
        )

    out: dict[str, str] = {}
    with sequence_to_group_csv.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            header = row.get("header", "").strip()
            final_key = row.get("final_group_key", "").strip()
            if header:
                out[header] = final_key
    return out


def write_csv(rows: list[dict[str, str]], out_file: Path) -> None:
    if not rows:
        LOGGER.warning("No rows to write for %s", out_file)
        return
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    LOGGER.info("Wrote %d rows: %s", len(rows), out_file)


def append_row_to_csv(row: dict[str, str], out_file: Path, fieldnames: list[str] | None = None) -> None:
    """Append a single row to a CSV file. Create file with headers if it doesn't exist."""
    out_file.parent.mkdir(parents=True, exist_ok=True)
    file_exists = out_file.exists()

    if fieldnames is None:
        fieldnames = list(row.keys())

    with out_file.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)
        handle.flush()
        os.fsync(handle.fileno())


def initialize_csv_file(out_file: Path, fieldnames: list[str]) -> None:
    """Create or reset a CSV file with its header row."""
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        handle.flush()
        os.fsync(handle.fileno())


def create_csv_backup(source: Path, backup: Path) -> None:
    """Create a backup copy of a CSV file. If source doesn't exist, create empty backup."""
    backup.parent.mkdir(parents=True, exist_ok=True)
    if source.exists():
        shutil.copy2(source, backup)
        LOGGER.info("Created backup: %s", backup)
    else:
        backup.touch()
        LOGGER.info("Created empty backup (source does not exist): %s", backup)


def append_row_to_per_donor_csv(
    row: dict[str, str],
    per_donor_dir: Path,
    fieldnames: list[str],
) -> None:
    donor_id = row.get("PAT id(SSAM)", "").strip()
    if not donor_id or donor_id == "-":
        donor_id = row.get("Patient Code", "").strip()
    if not donor_id or donor_id == "-":
        donor_id = "unknown_donor"
    donor_id = re.sub(r"[^A-Za-z0-9._-]+", "_", donor_id).strip("_")[:90]
    donor_csv = per_donor_dir / f"{donor_id}.csv"
    append_row_to_csv(row, donor_csv, fieldnames)


def build_merged_row(
    *,
    lanl_fields: dict[str, str],
    consensus_header: str,
    consensus_seq: str,
    query_seq: str,
    run_info: dict[str, str],
    sum_row: dict[str, str],
    hxb2_p_value: str = "",
    fasta_input_local: str = "",
    fasta_input_hxb2: str = "",
) -> dict[str, str]:
    con_trimmed, que_trimmed = _strip_flanking_gaps(consensus_seq, query_seq)
    mutation_counts = count_monomer_mutations(con_trimmed, que_trimmed)
    row: dict[str, str] = {}
    for field_name in LANL_HEADER_FIELDS:
        row[field_name] = lanl_fields.get(field_name, "")
    row["consensus_header"] = consensus_header
    for col, val in mutation_counts.items():
        row[col] = str(val)
    row.update({
        "mutation_from": run_info["mutation_from"],
        "mutation_to": run_info["mutation_to"],
        "match_mode": run_info["match"],
        "keepgaps": run_info["keepgaps"],
        "hypermut_p_value_local": sum_row["fisher_p"],
        "fasta_input_local": fasta_input_local,
        "hypermut_p_value_hxb2": hxb2_p_value,
        "fasta_input_hxb2": fasta_input_hxb2,
        "rate_ratio": sum_row["rate_ratio"],
        "primary_matches": sum_row["primary_matches"],
        "potential_primaries": sum_row["potential_primaries"],
        "control_matches": sum_row["control_matches"],
        "potential_controls": sum_row["potential_controls"],
    })
    return row


def process_existing_pair_task(
    *,
    cfg: PipelineConfig,
    pair_aligned: Path,
    seq_to_group: dict[str, str],
    mutation_directions: list[tuple[str, str]],
    run_modes: list[HypermutRunConfig],
    per_seq_dir: Path,
    hxb2_header: str = "",
    hxb2_seq: str = "",
    pair_input_hxb2_dir: Path | None = None,
    pair_aligned_hxb2_dir: Path | None = None,
) -> tuple[list[dict[str, str]], dict[str, str] | None]:
    rows: list[dict[str, str]] = []
    pair_records = read_fasta_records(pair_aligned)
    if len(pair_records) < 2:
        return rows, {
            "seq_header": pair_aligned.name,
            "final_group_key": "",
            "error": "Aligned pair FASTA has fewer than 2 records",
        }

    consensus_header = pair_records[0][0]
    consensus_seq = pair_records[0][1]
    seq_header = pair_records[1][0]
    query_seq = pair_records[1][1]
    LOGGER.info("Processing sequence (reuse mode): %s", seq_header)
    final_key = seq_to_group.get(seq_header, "")
    seq_slug = pair_aligned.name.replace("_aligned.fasta", "")
    seq_run_dir = per_seq_dir / seq_slug
    seq_run_dir.mkdir(parents=True, exist_ok=True)

    lanl_fields = _parse_lanl_header(seq_header)

    con_trimmed, que_trimmed = _strip_flanking_gaps(consensus_seq, query_seq)
    consensus_seq = con_trimmed
    query_seq = que_trimmed

    stripped_path = seq_run_dir / "stripped_aligned.fasta"
    write_fasta_records(
        [(consensus_header, consensus_seq), (seq_header, query_seq)],
        stripped_path,
    )

    qc = validate_alignment_quality(consensus_seq, query_seq)
    if not qc.passed:
        LOGGER.warning(
            "QC FAILED for %s: identity=%.3f, gap_frac=%.3f — skipping Hypermut3",
            seq_header, qc.identity, qc.gap_fraction_query,
        )
        return [], {
            "seq_header": seq_header,
            "final_group_key": final_key,
            "error": f"Alignment QC failed: identity={qc.identity:.3f}, gap_frac={qc.gap_fraction_query:.3f}",
        }

    try:
        hxb2_info: dict[tuple[str, str, str], tuple[str, str]] = {}
        if hxb2_seq and pair_input_hxb2_dir and pair_aligned_hxb2_dir:
            hxb2_start_str = lanl_fields.get("HXB2/MAC239 start", "").strip()
            hxb2_stop_str = lanl_fields.get("HXB2/MAC239 stop", "").strip()
            hxb2_start = int(hxb2_start_str) if hxb2_start_str not in ("", "-", "0") else 0
            hxb2_stop = int(hxb2_stop_str) if hxb2_stop_str not in ("", "-", "0") else 0

            if hxb2_start > 0 and hxb2_stop > hxb2_start:
                hxb2_region = hxb2_seq[hxb2_start - 1 : hxb2_stop]
            else:
                LOGGER.warning("Missing HXB2 coordinates for %s — skipping HXB2 p-value", seq_header)
                hxb2_region = ""

            if hxb2_region:
                pair_input_hxb2 = pair_input_hxb2_dir / f"{seq_slug}.fasta"
                pair_aligned_hxb2 = pair_aligned_hxb2_dir / f"{seq_slug}_aligned.fasta"
                hxb2_run_dir = seq_run_dir / "hxb2"
                hxb2_run_dir.mkdir(parents=True, exist_ok=True)
                try:
                    query_body = query_seq.replace("-", "")
                    aligned_hxb2 = align_pair_with_mafft(
                        consensus_header=hxb2_header,
                        consensus_seq=hxb2_region,
                        query_header=seq_header,
                        query_seq=query_body,
                        cfg=cfg,
                        pair_input_path=pair_input_hxb2,
                        pair_aligned_path=pair_aligned_hxb2,
                    )
                    write_fasta_records(aligned_hxb2, pair_aligned_hxb2)
                    strip_pair_fasta_flanks(pair_aligned_hxb2)
                    for mutation_from, mutation_to in mutation_directions:
                        direction_tag = f"{mutation_from}to{mutation_to}"
                        for mode in run_modes:
                            hxb2_run_info = run_hypermut_once(
                                cfg=cfg,
                                run_cfg=mode,
                                aligned_pair_fasta=pair_aligned_hxb2,
                                seq_output_dir=hxb2_run_dir,
                                mutation_from=mutation_from,
                                mutation_to=mutation_to,
                                direction_tag=f"hxb2_{direction_tag}",
                            )
                            hxb2_sum = read_single_hypermut_summary(Path(hxb2_run_info["summary_csv"]))
                            hxb2_info[(mutation_from, mutation_to, mode.name)] = (
                                hxb2_sum["fisher_p"], hxb2_run_info["input_fasta"]
                            )
                except Exception as exc:
                    LOGGER.warning("HXB2 Hypermut3 failed for %s: %s", seq_header, exc)

        for mutation_from, mutation_to in mutation_directions:
            direction_tag = f"{mutation_from}to{mutation_to}"
            for mode in run_modes:
                run_info = run_hypermut_once(
                    cfg=cfg,
                    run_cfg=mode,
                    aligned_pair_fasta=stripped_path,
                    seq_output_dir=seq_run_dir,
                    mutation_from=mutation_from,
                    mutation_to=mutation_to,
                    direction_tag=direction_tag,
                )
                sum_row = read_single_hypermut_summary(Path(run_info["summary_csv"]))
                hxb2_p, hxb2_fasta = hxb2_info.get((mutation_from, mutation_to, mode.name), ("", ""))
                rows.append(
                    build_merged_row(
                        lanl_fields=lanl_fields,
                        consensus_header=consensus_header,
                        consensus_seq=consensus_seq,
                        query_seq=query_seq,
                        run_info=run_info,
                        sum_row=sum_row,
                        hxb2_p_value=hxb2_p,
                        fasta_input_local=run_info["input_fasta"],
                        fasta_input_hxb2=hxb2_fasta,
                    )
                )
    except Exception as exc:
        return [], {
            "seq_header": seq_header,
            "final_group_key": final_key,
            "error": str(exc),
        }

    return rows, None


def process_sequence_task(
    *,
    cfg: PipelineConfig,
    idx: int,
    seq_header: str,
    final_key: str,
    seq_body: str,
    consensus_header: str,
    consensus_seq: str,
    mutation_directions: list[tuple[str, str]],
    run_modes: list[HypermutRunConfig],
    pair_input_dir: Path,
    pair_aligned_dir: Path,
    per_seq_dir: Path,
    hxb2_header: str = "",
    hxb2_seq: str = "",
    pair_input_hxb2_dir: Path | None = None,
    pair_aligned_hxb2_dir: Path | None = None,
) -> tuple[list[dict[str, str]], dict[str, str] | None]:
    rows: list[dict[str, str]] = []
    seq_slug = f"{idx:04d}_{sanitize_file_component(seq_header)}"
    LOGGER.info("Processing sequence %d: %s", idx, seq_header)
    pair_input = pair_input_dir / f"{seq_slug}.fasta"
    pair_aligned = pair_aligned_dir / f"{seq_slug}_aligned.fasta"
    seq_run_dir = per_seq_dir / seq_slug
    seq_run_dir.mkdir(parents=True, exist_ok=True)

    lanl_fields = _parse_lanl_header(seq_header)

    try:
        aligned_pair = align_pair_with_mafft(
            consensus_header=consensus_header,
            consensus_seq=consensus_seq,
            query_header=seq_header,
            query_seq=seq_body,
            cfg=cfg,
            pair_input_path=pair_input,
            pair_aligned_path=pair_aligned,
        )
        write_fasta_records(aligned_pair, pair_aligned)
        strip_pair_fasta_flanks(pair_aligned)
        aligned_pair = read_fasta_records(pair_aligned)

        aligned_consensus_seq = aligned_pair[0][1]
        aligned_query_seq = aligned_pair[1][1]

        qc = validate_alignment_quality(aligned_consensus_seq, aligned_query_seq)
        if not qc.passed:
            LOGGER.warning(
                "QC FAILED for %s: identity=%.3f, gap_frac=%.3f — skipping Hypermut3",
                seq_header, qc.identity, qc.gap_fraction_query,
            )
            return [], {
                "seq_header": seq_header,
                "final_group_key": final_key,
                "error": f"Alignment QC failed: identity={qc.identity:.3f}, gap_frac={qc.gap_fraction_query:.3f}",
            }

        hxb2_info: dict[tuple[str, str, str], tuple[str, str]] = {}
        if hxb2_seq and pair_input_hxb2_dir and pair_aligned_hxb2_dir:
            hxb2_start_str = lanl_fields.get("HXB2/MAC239 start", "").strip()
            hxb2_stop_str = lanl_fields.get("HXB2/MAC239 stop", "").strip()
            hxb2_start = int(hxb2_start_str) if hxb2_start_str not in ("", "-", "0") else 0
            hxb2_stop = int(hxb2_stop_str) if hxb2_stop_str not in ("", "-", "0") else 0

            if hxb2_start > 0 and hxb2_stop > hxb2_start:
                hxb2_region = hxb2_seq[hxb2_start - 1 : hxb2_stop]
            else:
                LOGGER.warning("Missing HXB2 coordinates for %s — skipping HXB2 p-value", seq_header)
                hxb2_region = ""

            if hxb2_region:
                pair_input_hxb2 = pair_input_hxb2_dir / f"{seq_slug}.fasta"
                pair_aligned_hxb2 = pair_aligned_hxb2_dir / f"{seq_slug}_aligned.fasta"
                hxb2_run_dir = seq_run_dir / "hxb2"
                hxb2_run_dir.mkdir(parents=True, exist_ok=True)
                try:
                    query_body = seq_body.replace("-", "")
                    aligned_hxb2 = align_pair_with_mafft(
                        consensus_header=hxb2_header,
                        consensus_seq=hxb2_region,
                        query_header=seq_header,
                        query_seq=query_body,
                        cfg=cfg,
                        pair_input_path=pair_input_hxb2,
                        pair_aligned_path=pair_aligned_hxb2,
                    )
                    write_fasta_records(aligned_hxb2, pair_aligned_hxb2)
                    strip_pair_fasta_flanks(pair_aligned_hxb2)
                    for mutation_from, mutation_to in mutation_directions:
                        direction_tag = f"{mutation_from}to{mutation_to}"
                        for mode in run_modes:
                            hxb2_run_info = run_hypermut_once(
                                cfg=cfg,
                                run_cfg=mode,
                                aligned_pair_fasta=pair_aligned_hxb2,
                                seq_output_dir=hxb2_run_dir,
                                mutation_from=mutation_from,
                                mutation_to=mutation_to,
                                direction_tag=f"hxb2_{direction_tag}",
                            )
                            hxb2_sum = read_single_hypermut_summary(Path(hxb2_run_info["summary_csv"]))
                            hxb2_info[(mutation_from, mutation_to, mode.name)] = (
                                hxb2_sum["fisher_p"], hxb2_run_info["input_fasta"]
                            )
                except Exception as exc:
                    LOGGER.warning("HXB2 Hypermut3 failed for %s: %s", seq_header, exc)

        for mutation_from, mutation_to in mutation_directions:
            direction_tag = f"{mutation_from}to{mutation_to}"
            for mode in run_modes:
                run_info = run_hypermut_once(
                    cfg=cfg,
                    run_cfg=mode,
                    aligned_pair_fasta=pair_aligned,
                    seq_output_dir=seq_run_dir,
                    mutation_from=mutation_from,
                    mutation_to=mutation_to,
                    direction_tag=direction_tag,
                )
                sum_row = read_single_hypermut_summary(Path(run_info["summary_csv"]))
                hxb2_p, hxb2_fasta = hxb2_info.get((mutation_from, mutation_to, mode.name), ("", ""))
                rows.append(
                    build_merged_row(
                        lanl_fields=lanl_fields,
                        consensus_header=consensus_header,
                        consensus_seq=aligned_consensus_seq,
                        query_seq=aligned_query_seq,
                        run_info=run_info,
                        sum_row=sum_row,
                        hxb2_p_value=hxb2_p,
                        fasta_input_local=run_info["input_fasta"],
                        fasta_input_hxb2=hxb2_fasta,
                    )
                )
    except Exception as exc:
        return [], {
            "seq_header": seq_header,
            "final_group_key": final_key,
            "error": str(exc),
        }

    return rows, None


# %% Cell 5 - Pairwise Hypermut3 pipeline
def run_pipeline(cfg: PipelineConfig) -> dict[str, str | int]:
    LOGGER.info("Pairwise Hypermut3 pipeline started")
    LOGGER.info("Group manifest CSV: %s", cfg.group_manifest_csv)
    LOGGER.info("Sequence mapping CSV: %s", cfg.sequence_to_group_csv)
    LOGGER.info("Output directory: %s", cfg.output_dir)

    if not cfg.hypermut_script.exists():
        raise FileNotFoundError(f"Hypermut script not found: {cfg.hypermut_script}")
    if not shutil.which(cfg.mafft_binary):
        raise RuntimeError(f"MAFFT binary not found: {cfg.mafft_binary}")

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    pair_input_dir = cfg.output_dir / "pairwise_inputs"
    pair_aligned_dir = cfg.output_dir / "pairwise_aligned"
    pair_input_hxb2_dir = cfg.output_dir / "pairwise_inputs_hxb2"
    pair_aligned_hxb2_dir = cfg.output_dir / "pairwise_aligned_hxb2"
    per_seq_dir = cfg.output_dir / "per_sequence_hypermut"
    per_donor_dir = cfg.output_dir / "per_donor"
    for path in [pair_input_dir, pair_aligned_dir, pair_input_hxb2_dir,
                 pair_aligned_hxb2_dir, per_seq_dir, per_donor_dir]:
        path.mkdir(parents=True, exist_ok=True)

    hxb2_header = ""
    hxb2_seq = ""
    if cfg.hxb2_reference_fasta and cfg.hxb2_reference_fasta.exists():
        hxb2_records = read_fasta_records(cfg.hxb2_reference_fasta)
        if hxb2_records:
            hxb2_header = hxb2_records[0][0]
            hxb2_seq = normalize_sequence(hxb2_records[0][1])
            LOGGER.info("Loaded HXB2 reference: %s (%d bp)", hxb2_header, len(hxb2_seq))
        else:
            LOGGER.warning("HXB2 FASTA is empty: %s", cfg.hxb2_reference_fasta)
    else:
        LOGGER.warning("HXB2 reference not provided or not found — skipping HXB2 p-values")

    merged_csv = cfg.output_dir / "per_sequence_hypermut_merged.csv"
    live_progress_csv = (
        cfg.output_dir / "per_sequence_hypermut_merged_for_initial_showing.csv"
    )
    create_csv_backup(merged_csv, live_progress_csv)

    seq_to_group = load_sequence_to_group_map(cfg.sequence_to_group_csv)
    group_consensus: dict[str, tuple[str, str]] = {}
    if cfg.reuse_existing_pairwise_alignments:
        LOGGER.info(
            "Loaded %d sequence mappings (reuse mode: pairwise alignments already exist)",
            len(seq_to_group),
        )
    else:
        group_consensus = load_group_consensus_map(cfg.group_manifest_csv)
        LOGGER.info(
            "Loaded %d group consensus entries and %d sequence mappings",
            len(group_consensus),
            len(seq_to_group),
        )

    if cfg.full_run:
        run_modes = [
            HypermutRunConfig(name="strict-keepgaps", match="strict", keepgaps=True),
            HypermutRunConfig(name="strict-skipgaps", match="strict", keepgaps=False),
            HypermutRunConfig(name="partial-keepgaps", match="partial", keepgaps=True),
            HypermutRunConfig(name="partial-skipgaps", match="partial", keepgaps=False),
        ]
        LOGGER.info("Full run mode: all 4 Hypermut3 configurations enabled")
    else:
        run_modes = [
            HypermutRunConfig(name="strict-keepgaps", match="strict", keepgaps=True),
        ]
        LOGGER.info("Default run mode: strict-keepgaps only (use --full-run true for all 4 modes + all mutation directions)")
    mutation_directions = get_mutation_directions(cfg)
    LOGGER.info("Mutation directions configured: %d", len(mutation_directions))
    if cfg.run_all_mutation_directions:
        LOGGER.info("All-direction mode enabled (12 directions expected)")

    workers = cfg.parallel_workers or max(1, (os.cpu_count() or 1))
    workers = max(1, workers)
    if workers > 1 and cfg.mafft_threads == -1 and not cfg.reuse_existing_pairwise_alignments:
        LOGGER.warning(
            "parallel_workers=%d with mafft_threads=-1 can oversubscribe CPUs; "
            "forcing mafft_threads=1 for balanced parallel execution",
            workers,
        )
        cfg.mafft_threads = 1
    LOGGER.info("Parallel workers: %d", workers)

    sequence_records: list[tuple[str, str, str]] = []
    for seq_header, final_key in seq_to_group.items():
        sequence_records.append((seq_header, "", final_key))
    if cfg.number_of_sequence is not None:
        sequence_records = sequence_records[: cfg.number_of_sequence]
        LOGGER.info(
            "Sequence limit enabled: running first %d sequence records",
            len(sequence_records),
        )

    query_map: dict[str, str] = {}
    if not cfg.reuse_existing_pairwise_alignments:
        # Reconstruct query sequences from group aligned FASTA files listed
        # in group manifest. Fallback to raw FASTA if aligned file is missing.
        # Remove alignment gaps before running pairwise MAFFT.
        with cfg.group_manifest_csv.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                full_key = row.get("full_key", "").strip()
                aligned_fasta = resolve_manifest_path(
                    row.get("aligned_fasta", ""),
                    cfg.group_manifest_csv,
                )
                raw_fasta = resolve_manifest_path(
                    row.get("raw_fasta", ""),
                    cfg.group_manifest_csv,
                )
                if not full_key:
                    continue

                source_fasta = aligned_fasta if aligned_fasta.exists() else raw_fasta
                if not source_fasta.exists():
                    continue

                for header, seq in read_fasta_records(source_fasta):
                    query_map[header] = seq.replace("-", "")

    rows_merged: list[dict[str, str]] = []
    rows_failures: list[dict[str, str]] = []

    merged_csv_fieldnames = (
        list(LANL_HEADER_FIELDS)
        + ["consensus_header"]
        + MONOMER_MUTATION_COLUMNS
        + ["mutation_from", "mutation_to", "match_mode", "keepgaps"]
        + ["hypermut_p_value_local", "fasta_input_local",
           "hypermut_p_value_hxb2", "fasta_input_hxb2",
           "rate_ratio", "primary_matches",
           "potential_primaries", "control_matches", "potential_controls"]
    )
    initialize_csv_file(merged_csv, merged_csv_fieldnames)
    initialize_csv_file(live_progress_csv, merged_csv_fieldnames)

    total = len(sequence_records)
    if cfg.reuse_existing_pairwise_alignments:
        LOGGER.info(
            "Starting Hypermut3-only execution using existing pairwise aligned files"
        )
    else:
        LOGGER.info(
            "Starting per-sequence pairwise alignment + Hypermut3 for %d sequences",
            total,
        )

    if cfg.reuse_existing_pairwise_alignments:
        aligned_files = sorted(pair_aligned_dir.glob("*_aligned.fasta"))
        if cfg.number_of_sequence is not None:
            aligned_files = aligned_files[: cfg.number_of_sequence]
            LOGGER.info(
                "Sequence limit enabled: running first %d aligned FASTA files",
                len(aligned_files),
            )

        if cfg.problematic_only:
            filtered_files: list[Path] = []
            for fasta_path in aligned_files:
                records = read_fasta_records(fasta_path)
                if len(records) < 2:
                    continue
                seq_header = records[1][0]
                if _is_problematic(seq_header):
                    filtered_files.append(fasta_path)
                else:
                    LOGGER.info(
                        "Skipping non-problematic sequence: %s",
                        seq_header[:120],
                    )
            LOGGER.info(
                "Problematic-only filter: %d / %d sequences retained",
                len(filtered_files), len(aligned_files),
            )
            aligned_files = filtered_files

        total = len(aligned_files)
        LOGGER.info(
            "Hypermut-only mode: reusing %d existing pairwise aligned FASTA files (no MAFFT run)",
            total,
        )

        completed = 0
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(
                    process_existing_pair_task,
                    cfg=cfg,
                    pair_aligned=pair_aligned,
                    seq_to_group=seq_to_group,
                    mutation_directions=mutation_directions,
                    run_modes=run_modes,
                    per_seq_dir=per_seq_dir,
                    hxb2_header=hxb2_header,
                    hxb2_seq=hxb2_seq,
                    pair_input_hxb2_dir=pair_input_hxb2_dir,
                    pair_aligned_hxb2_dir=pair_aligned_hxb2_dir,
                )
                for pair_aligned in aligned_files
            ]
            for future in as_completed(futures):
                rows_batch, failure = future.result()
                completed += 1
                if failure is not None:
                    rows_failures.append(failure)
                    LOGGER.info(
                        "Completed %d/%d (FAILED): %s | error=%s",
                        completed,
                        total,
                        failure.get("seq_header", ""),
                        failure.get("error", ""),
                    )
                    continue
                for row_data in rows_batch:
                    append_row_to_csv(row_data, merged_csv, merged_csv_fieldnames)
                    append_row_to_csv(row_data, live_progress_csv, merged_csv_fieldnames)
                    append_row_to_per_donor_csv(row_data, per_donor_dir, merged_csv_fieldnames)
                rows_merged.extend(rows_batch)
                seq_name = rows_batch[0].get("Accession", "") if rows_batch else ""
                LOGGER.info(
                    "Completed %d/%d: %s | runs=%d",
                    completed,
                    total,
                    seq_name,
                    len(rows_batch),
                )

    else:
        valid_tasks: list[tuple[int, str, str, str, str, str]] = []
        for idx, (seq_header, _, final_key) in enumerate(sequence_records, start=1):
            if cfg.problematic_only and not _is_problematic(seq_header):
                LOGGER.info(
                    "Skipping non-problematic sequence: %s",
                    seq_header[:120],
                )
                continue

            seq_body = query_map.get(seq_header)
            if seq_body is None:
                rows_failures.append(
                    {
                        "seq_header": seq_header,
                        "final_group_key": final_key,
                        "error": "Sequence not found in aligned/raw group FASTA map",
                    }
                )
                continue

            consensus_entry = group_consensus.get(final_key)
            if consensus_entry is None:
                rows_failures.append(
                    {
                        "seq_header": seq_header,
                        "final_group_key": final_key,
                        "error": "Consensus not found for final group key",
                    }
                )
                continue

            consensus_header, consensus_seq = consensus_entry
            valid_tasks.append(
                (idx, seq_header, final_key, seq_body, consensus_header, consensus_seq)
            )

        total_valid = len(valid_tasks)
        LOGGER.info(
            "Prepared %d valid sequence tasks (pre-check failures: %d)",
            total_valid,
            len(rows_failures),
        )

        completed = 0
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(
                    process_sequence_task,
                    cfg=cfg,
                    idx=idx,
                    seq_header=seq_header,
                    final_key=final_key,
                    seq_body=seq_body,
                    consensus_header=consensus_header,
                    consensus_seq=consensus_seq,
                    mutation_directions=mutation_directions,
                    run_modes=run_modes,
                    pair_input_dir=pair_input_dir,
                    pair_aligned_dir=pair_aligned_dir,
                    per_seq_dir=per_seq_dir,
                    hxb2_header=hxb2_header,
                    hxb2_seq=hxb2_seq,
                    pair_input_hxb2_dir=pair_input_hxb2_dir,
                    pair_aligned_hxb2_dir=pair_aligned_hxb2_dir,
                )
                for idx, seq_header, final_key, seq_body, consensus_header, consensus_seq in valid_tasks
            ]
            for future in as_completed(futures):
                rows_batch, failure = future.result()
                completed += 1
                if failure is not None:
                    rows_failures.append(failure)
                    LOGGER.info(
                        "Completed %d/%d (FAILED): %s | error=%s",
                        completed,
                        total_valid,
                        failure.get("seq_header", ""),
                        failure.get("error", ""),
                    )
                    continue
                for row_data in rows_batch:
                    append_row_to_csv(row_data, merged_csv, merged_csv_fieldnames)
                    append_row_to_csv(row_data, live_progress_csv, merged_csv_fieldnames)
                    append_row_to_per_donor_csv(row_data, per_donor_dir, merged_csv_fieldnames)
                rows_merged.extend(rows_batch)
                seq_name = rows_batch[0].get("Accession", "") if rows_batch else ""
                LOGGER.info(
                    "Completed %d/%d: %s | runs=%d",
                    completed,
                    total_valid,
                    seq_name,
                    len(rows_batch),
                )

    LOGGER.info("Merged CSV was updated incrementally during processing: %s", merged_csv)

    failed_csv = cfg.output_dir / "per_sequence_hypermut_failures.csv"
    write_csv(rows_failures, failed_csv)

    LOGGER.info("Pairwise Hypermut3 pipeline finished")
    LOGGER.info("Success rows: %d", len(rows_merged))
    LOGGER.info("Failure rows: %d", len(rows_failures))

    return {
        "group_manifest_csv": str(cfg.group_manifest_csv),
        "sequence_to_group_csv": str(cfg.sequence_to_group_csv),
        "output_dir": str(cfg.output_dir),
        "merged_csv": str(merged_csv),
        "live_progress_csv": str(live_progress_csv),
        "failures_csv": str(failed_csv),
        "n_sequences": total,
        "n_success_rows": len(rows_merged),
        "n_failure_rows": len(rows_failures),
    }


# %% Cell 6 - Entry point
def parse_bool(value: str) -> bool:
    text = (value or "").strip().lower()
    if text in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "f", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(
        "Expected a boolean value: true/false, yes/no, 1/0"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run pairwise Hypermut3 pipeline with optional sequence limit"
    )
    parser.add_argument(
        "--number-of-sequence",
        type=int,
        default=None,
        help="Run only the first N sequences (or first N aligned files in reuse mode)",
    )
    parser.add_argument(
        "--number-of-sequences",
        dest="number_of_sequence",
        type=int,
        default=None,
        help="Alias for --number-of-sequence",
    )
    parser.add_argument(
        "--reuse-existing-pairwise-alignments",
        type=parse_bool,
        default=False,
        help=(
            "Reuse precomputed pairwise aligned FASTA files from "
            "data/hypermute3_pairwise_results/pairwise_aligned. "
            "Default is false (fresh pairwise MAFFT + Hypermut run)."
        ),
    )
    parser.add_argument(
        "--parallel-workers",
        type=int,
        default=max(1, (os.cpu_count() or 1)),
        help="Number of parallel sequence workers (default: CPU core count)",
    )
    parser.add_argument(
        "--hxb2-reference",
        type=Path,
        default=None,
        help="Path to HXB2 reference FASTA file. If not provided, auto-downloads from LANL.",
    )
    parser.add_argument(
        "--full-run",
        type=parse_bool,
        default=False,
        help=(
            "Run all 4 Hypermut3 modes (strict/partial × keepgaps/skipgaps) "
            "and all 12 mutation directions. Default is false "
            "(G→A, strict, keepgaps only)."
        ),
    )
    parser.add_argument(
        "--problematic-only",
        type=parse_bool,
        default=False,
        help=(
            "Only process sequences marked as problematic (value 3). "
            "Default is false (process all sequences)."
        ),
    )
    args = parser.parse_args()
    if args.number_of_sequence is not None and args.number_of_sequence <= 0:
        parser.error("--number-of-sequence must be a positive integer")
    if args.parallel_workers is not None and args.parallel_workers <= 0:
        parser.error("--parallel-workers must be a positive integer")
    return args


def main() -> None:
    args = parse_args()
    project_root = Path(__file__).resolve().parent

    if args.hxb2_reference:
        hxb2_fasta = args.hxb2_reference
    else:
        try:
            hxb2_fasta = ensure_hxb2_reference(project_root / "data")
        except RuntimeError as exc:
            LOGGER.warning("%s — HXB2 p-values will be empty", exc)
            hxb2_fasta = None

    cfg = PipelineConfig(
        group_manifest_csv=project_root
        / "data"
        / "processed_grouped"
        / "group_manifest.csv",
        sequence_to_group_csv=project_root
        / "data"
        / "processed_grouped"
        / "sequence_to_group_mapping.csv",
        hypermut_script=project_root / "Hypermut3" / "hypermut.py",
        output_dir=project_root / "data" / "hypermute3_pairwise_results",
        python_executable=sys.executable,
        mafft_binary="mafft",
        mafft_threads=-1,
        mutation_from="G",
        mutation_to="A",
        upstream_context="",
        downstream_context="RD",
        enforce="D",
        begin=0,
        finish=None,
        number_of_sequence=args.number_of_sequence,
        run_all_mutation_directions=args.full_run,
        reuse_existing_pairwise_alignments=args.reuse_existing_pairwise_alignments,
        parallel_workers=args.parallel_workers,
        hxb2_reference_fasta=hxb2_fasta,
        log_level="INFO",
        full_run=args.full_run,
        problematic_only=args.problematic_only,
    )

    configure_logging(cfg.log_level)
    results = run_pipeline(cfg)

    print("Pairwise Hypermut3 completed. Key outputs:")
    for key, value in results.items():
        print(f"- {key}: {value}")


if __name__ == "__main__":
    main()

# ──────────────────────────────────────────────────────────────
# Usage:
#
#   Default (G→A, strict, keepgaps — 1 check per sequence):
#     python 02_hypermut3.py
#
#   Full run (all 12 mutation directions × all 4 modes):
#     python 02_hypermut3.py --full-run true
#
#   Problematic sequences only (skip non-problematic):
#     python 02_hypermut3.py --problematic-only true
#
#   Other flags:
#     --number-of-sequence N                          Run first N sequences only
#     --parallel-workers N                            Number of parallel workers
#     --reuse-existing-pairwise-alignments true       Reuse existing pairwise FASTA
#     --hxb2-reference path/to/hxb2.fasta             Custom HXB2 reference
# ──────────────────────────────────────────────────────────────

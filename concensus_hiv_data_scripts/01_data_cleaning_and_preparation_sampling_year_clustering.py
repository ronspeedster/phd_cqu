from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import logging
import shutil
import subprocess
import tempfile
from collections import Counter
from typing import Iterable


LOGGER = logging.getLogger("grouped_preparation")
DNA_STRICT_ALLOWED = set("ACGT-")

HEADER_FIELDS = [
    "Se ID",
    "Patient Code",
    "PAT id(SSAM)",
    "Accession",
    "Name",
    "Subtype",
    "Country",
    "Sampling Year",
    "Problematic Sequence",
    "HXB2/MAC239 start",
    "HXB2/MAC239 stop",
    "Sequence Length",
    "Organism",
]

YEAR_FIELD = "Sampling Year"
START_FIELD = "HXB2/MAC239 start"


@dataclass
class SequenceRecord:
    header: str
    sequence: str


@dataclass
class GroupedConfig:
    input_fasta: Path
    output_dir: Path
    min_bin_size: int = 20
    force_mafft_alignment: bool = True
    mafft_binary: str = "mafft"
    mafft_threads: int = -1
    vsearch_binary: str = "vsearch"
    vsearch_identity: float = 0.95
    max_bin_size: int = 30
    unknown_year_handling: str = "separate"  # separate | discard
    consensus_tie_priority: str = "ACGT"
    log_level: str = "INFO"


def configure_logging(level_name: str) -> None:
    level = getattr(logging, (level_name or "INFO").upper(), logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def read_fasta(file_path: Path) -> list[SequenceRecord]:
    LOGGER.info("Reading FASTA: %s", file_path)
    records: list[SequenceRecord] = []
    header: str | None = None
    seq_parts: list[str] = []

    with file_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(
                        SequenceRecord(header=header, sequence="".join(seq_parts))
                    )
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        records.append(SequenceRecord(header=header, sequence="".join(seq_parts)))

    LOGGER.info("Loaded %d sequences from FASTA", len(records))
    return records


def write_fasta(
    records: Iterable[SequenceRecord],
    file_path: Path,
    line_width: int = 80,
) -> None:
    with file_path.open("w", encoding="utf-8") as handle:
        for rec in records:
            handle.write(f">{rec.header}\n")
            seq = rec.sequence
            for idx in range(0, len(seq), line_width):
                handle.write(seq[idx : idx + line_width] + "\n")


def normalize_sequence(seq: str) -> str:
    seq = seq.upper().replace("U", "T")
    return "".join(base if base in DNA_STRICT_ALLOWED else "-" for base in seq)


def is_aligned(records: list[SequenceRecord]) -> bool:
    if not records:
        return False
    lengths = {len(rec.sequence) for rec in records}
    return len(lengths) == 1


def parse_header_fields(header: str) -> dict[str, str | None]:
    values = [part.strip() for part in header.split(",")]
    parsed: dict[str, str | None] = {}

    for idx, field in enumerate(HEADER_FIELDS):
        parsed[field] = values[idx] if idx < len(values) and values[idx] else None

    return parsed


def safe_int(text: str | None) -> int | None:
    if text is None:
        return None
    try:
        return int(float(text))
    except (TypeError, ValueError):
        return None


def initial_group_key(header: str) -> tuple[int | None, int | None]:
    fields = parse_header_fields(header)
    year = safe_int(fields.get(YEAR_FIELD))
    start = safe_int(fields.get(START_FIELD))
    if year is None or start is None:
        return (None, None)
    return (year, start)


def key_sort_value(key: tuple[int | None, int | None]) -> tuple[int, int]:
    year, start = key
    if year is None or start is None:
        return (10**9, 10**9)
    return (year, start)


def key_to_text(key: tuple[int | None, int | None]) -> str:
    year, start = key
    year_text = str(year) if year is not None else "UNKNOWN"
    start_text = str(start) if start is not None else "UNKNOWN"
    return f"SamplingYear={year_text}|HXB2Start={start_text}"


def build_initial_bins(
    records: list[SequenceRecord],
) -> tuple[dict[tuple[int | None, int | None], list[SequenceRecord]], list[dict[str, str]]]:
    LOGGER.info(
        "Building initial bins using %s + %s",
        YEAR_FIELD,
        START_FIELD,
    )
    bins: dict[tuple[int | None, int | None], list[SequenceRecord]] = {}
    seq_map_rows: list[dict[str, str]] = []

    for rec in records:
        key = initial_group_key(rec.header)
        bins.setdefault(key, []).append(rec)
        seq_map_rows.append(
            {
                "header": rec.header,
                "initial_group_key": key_to_text(key),
            }
        )

    unknown_count = len(bins.get((None, None), []))
    LOGGER.info(
        "Initial binning complete: %d bins (%d UNKNOWN sequences)",
        len(bins),
        unknown_count,
    )
    return bins, seq_map_rows


def build_year_bins(
    records: list[SequenceRecord],
) -> tuple[dict[int | None, list[SequenceRecord]], list[dict[str, str]]]:
    LOGGER.info("Building initial bins using %s only", YEAR_FIELD)
    bins: dict[int | None, list[SequenceRecord]] = {}
    seq_map_rows: list[dict[str, str]] = []

    for rec in records:
        fields = parse_header_fields(rec.header)
        year = safe_int(fields.get(YEAR_FIELD))
        bins.setdefault(year, []).append(rec)
        year_text = str(year) if year is not None else "UNKNOWN"
        seq_map_rows.append(
            {
                "header": rec.header,
                "initial_year_key": f"SamplingYear={year_text}",
            }
        )

    unknown_count = len(bins.get(None, []))
    LOGGER.info(
        "Year binning complete: %d bins (%d UNKNOWN sequences)",
        len(bins),
        unknown_count,
    )
    return bins, seq_map_rows


def cluster_sequences_vsearch(
    records: list[SequenceRecord],
    cfg: GroupedConfig,
    year_label: str,
) -> list[tuple[str, list[SequenceRecord]]]:
    if not records:
        return []

    vsearch_exe = shutil.which(cfg.vsearch_binary)
    if not vsearch_exe:
        raise RuntimeError(
            "vsearch was not found on PATH. Install vsearch to enable year-wise clustering."
        )

    with tempfile.TemporaryDirectory() as tmpdir:
        input_fasta = Path(tmpdir) / "input.fasta"
        output_uc = Path(tmpdir) / "clusters.uc"

        write_fasta(records, input_fasta)
        cmd = [
            vsearch_exe,
            "--cluster_fast",
            str(input_fasta),
            "--id",
            str(cfg.vsearch_identity),
            "--uc",
            str(output_uc),
        ]

        LOGGER.info(
            "[%s] Running vsearch clustering (id=%s)",
            year_label,
            cfg.vsearch_identity,
        )
        proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
        if proc.returncode != 0:
            raise RuntimeError(
                f"vsearch clustering failed ({proc.returncode}): {proc.stderr}"
            )

        header_map = {rec.header: rec for rec in records}
        clusters: dict[str, list[SequenceRecord]] = {}
        if output_uc.exists():
            with output_uc.open("r", encoding="utf-8") as handle:
                for line in handle:
                    if not line:
                        continue
                    rec_type = line[0]
                    if rec_type not in {"S", "H"}:
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 9:
                        continue
                    cluster_id = parts[1]
                    seq_header = parts[8]
                    seq_rec = header_map.get(seq_header)
                    if seq_rec is not None:
                        clusters.setdefault(cluster_id, []).append(seq_rec)

        # Safety net: if vsearch output parsing misses any headers, keep all records.
        if not clusters:
            return [("0", records)]

        ordered = sorted(clusters.items(), key=lambda item: int(item[0]))
        return ordered


def split_cluster_to_bins(
    cluster_records: list[SequenceRecord],
    max_size: int,
) -> list[list[SequenceRecord]]:
    bins: list[list[SequenceRecord]] = []
    for idx in range(0, len(cluster_records), max_size):
        bins.append(cluster_records[idx : idx + max_size])

    # Avoid trailing singleton bins (e.g., 31 -> 30+1 becomes 31).
    if len(bins) >= 2 and len(bins[-1]) == 1:
        bins[-2].extend(bins[-1])
        bins.pop()

    return bins


def rebalance_bins_by_similarity(
    bins: list[list[SequenceRecord]],
    min_size: int,
    max_size: int,
) -> list[list[SequenceRecord]]:
    if not bins:
        return []

    def centroid(records: list[SequenceRecord]) -> str:
        return records[0].sequence if records else ""

    def similarity(seq1: str, seq2: str) -> float:
        length = min(len(seq1), len(seq2))
        matches = sum(1 for idx in range(length) if seq1[idx] == seq2[idx])
        return matches / length if length > 0 else 0.0

    working = [list(group) for group in bins if group]
    stable = [group for group in working if len(group) >= min_size]
    unassigned: list[SequenceRecord] = []
    for group in working:
        if len(group) < min_size:
            unassigned.extend(group)

    def best_target_for_sequence(seq: SequenceRecord) -> int | None:
        best_idx = None
        best_score = -1.0
        for idx, group in enumerate(stable):
            if len(group) >= max_size:
                continue
            score = similarity(seq.sequence, centroid(group))
            if score > best_score:
                best_score = score
                best_idx = idx
        return best_idx

    carry: list[SequenceRecord] = []
    for rec in unassigned:
        target_idx = best_target_for_sequence(rec)
        if target_idx is None:
            carry.append(rec)
            continue
        stable[target_idx].append(rec)

    if carry:
        for idx in range(0, len(carry), max_size):
            stable.append(carry[idx : idx + max_size])

    # Repair undersized bins by redistributing their sequences to nearest bins.
    changed = True
    while changed and len(stable) > 1:
        changed = False
        undersized_indices = [idx for idx, group in enumerate(stable) if len(group) < min_size]
        if not undersized_indices:
            break

        source_idx = undersized_indices[0]
        source_group = stable[source_idx]
        moved_any = False
        leftovers: list[SequenceRecord] = []

        for rec in source_group:
            best_idx = None
            best_score = -1.0
            for idx, group in enumerate(stable):
                if idx == source_idx or len(group) >= max_size:
                    continue
                score = similarity(rec.sequence, centroid(group))
                if score > best_score:
                    best_score = score
                    best_idx = idx
            if best_idx is None:
                leftovers.append(rec)
            else:
                stable[best_idx].append(rec)
                moved_any = True

        if moved_any:
            stable.pop(source_idx)
            if leftovers:
                stable.append(leftovers)
            changed = True

    # Final singleton prevention from new bins if possible.
    if len(stable) >= 2 and len(stable[-1]) == 1 and len(stable[-2]) < max_size:
        stable[-2].extend(stable[-1])
        stable.pop()

    return stable


def merge_small_clusters(
    clusters: list[tuple[str, list[SequenceRecord]]],
    min_size: int,
) -> list[tuple[str, list[SequenceRecord]]]:
    if not clusters:
        return []

    large: list[tuple[str, list[SequenceRecord]]] = []
    small: list[tuple[str, list[SequenceRecord]]] = []

    for cid, recs in clusters:
        if len(recs) >= min_size:
            large.append((cid, recs))
        else:
            small.append((cid, recs))

    def get_centroid(rec_list: list[SequenceRecord]) -> str:
        return rec_list[0].sequence

    def similarity(seq1: str, seq2: str) -> float:
        length = min(len(seq1), len(seq2))
        matches = sum(1 for idx in range(length) if seq1[idx] == seq2[idx])
        return matches / length if length > 0 else 0.0

    # If no large cluster exists, iteratively merge the smallest cluster into
    # its nearest neighbor until all clusters reach min_size or only one cluster remains.
    # This handles years where all initial vsearch clusters are small.
    if not large:
        working = [[cid, recs] for cid, recs in small]

        if len(working) == 1:
            return [(str(working[0][0]), working[0][1])]

        while True:
            small_indices = [idx for idx, (_cid, recs) in enumerate(working) if len(recs) < min_size]
            if not small_indices or len(working) == 1:
                break

            source_idx = min(small_indices, key=lambda idx: len(working[idx][1]))
            source_cid, source_recs = working[source_idx]
            source_centroid = get_centroid(source_recs)

            best_idx = None
            best_score = -1.0
            for idx, (_target_cid, target_recs) in enumerate(working):
                if idx == source_idx:
                    continue
                score = similarity(source_centroid, get_centroid(target_recs))
                if score > best_score:
                    best_score = score
                    best_idx = idx

            if best_idx is None:
                break

            working[best_idx][1].extend(source_recs)
            working.pop(source_idx)

        return [(str(cid), recs) for cid, recs in working]

    large_centroids: list[list[object]] = [
        [cid, get_centroid(recs), recs] for cid, recs in large
    ]

    for _small_id, small_recs in small:
        small_centroid = get_centroid(small_recs)
        best_idx = 0
        best_score = -1.0

        for idx, (_cid, large_centroid, _recs) in enumerate(large_centroids):
            score = similarity(small_centroid, str(large_centroid))
            if score > best_score:
                best_score = score
                best_idx = idx

        large_centroids[best_idx][2].extend(small_recs)

    return [(str(cid), recs) for cid, _centroid, recs in large_centroids]


def merge_small_known_bins(
    bins: dict[tuple[int | None, int | None], list[SequenceRecord]],
    min_size: int,
) -> tuple[list[tuple[int | None, int | None]], dict[tuple[int | None, int | None], list[str]]]:
    LOGGER.info("Merging known bins with minimum bin size = %d", min_size)
    known_keys = [k for k in bins if k != (None, None)]
    known_keys.sort(key=key_sort_value)
    merge_trace: dict[tuple[int | None, int | None], list[str]] = {
        key: [key_to_text(key)] for key in known_keys
    }

    idx = 0
    while idx < len(known_keys) - 1:
        key = known_keys[idx]
        if len(bins[key]) < min_size:
            next_key = known_keys[idx + 1]
            LOGGER.info(
                "Merging small bin %s (n=%d) into next bin %s (n=%d)",
                key_to_text(key),
                len(bins[key]),
                key_to_text(next_key),
                len(bins[next_key]),
            )
            bins[next_key] = bins[key] + bins[next_key]
            merge_trace.setdefault(next_key, [key_to_text(next_key)]).extend(
                merge_trace.get(key, [key_to_text(key)])
            )
            bins.pop(key, None)
            merge_trace.pop(key, None)
            known_keys.pop(idx)
            continue
        idx += 1

    if len(known_keys) >= 2:
        last_key = known_keys[-1]
        if len(bins[last_key]) < min_size:
            prev_key = known_keys[-2]
            LOGGER.info(
                "Merging trailing small bin %s (n=%d) into previous bin %s (n=%d)",
                key_to_text(last_key),
                len(bins[last_key]),
                key_to_text(prev_key),
                len(bins[prev_key]),
            )
            bins[prev_key] = bins[prev_key] + bins[last_key]
            merge_trace.setdefault(prev_key, [key_to_text(prev_key)]).extend(
                merge_trace.get(last_key, [key_to_text(last_key)])
            )
            bins.pop(last_key, None)
            merge_trace.pop(last_key, None)
            known_keys.pop()

    final_keys = sorted(known_keys, key=key_sort_value)
    if (None, None) in bins:
        final_keys.append((None, None))
        merge_trace[(None, None)] = [key_to_text((None, None))]
        LOGGER.info(
            "Keeping UNKNOWN bin as separate final group (n=%d)",
            len(bins[(None, None)]),
        )

    LOGGER.info("Final groups after merging: %d", len(final_keys))

    return final_keys, merge_trace


def align_group_records(
    group_records: list[SequenceRecord],
    cfg: GroupedConfig,
    group_tmp_input: Path,
    group_aligned_output: Path,
    group_label: str,
) -> list[SequenceRecord]:
    LOGGER.info(
        "[%s] Preparing group alignment for %d sequences",
        group_label,
        len(group_records),
    )
    normalized = [
        SequenceRecord(header=r.header, sequence=normalize_sequence(r.sequence))
        for r in group_records
    ]

    # MAFFT requires at least 2 sequences. Keep singleton bins as-is.
    if len(normalized) == 1:
        LOGGER.warning(
            "[%s] Singleton bin detected; MAFFT requires >=2 sequences, using input as aligned output",
            group_label,
        )
        write_fasta(normalized, group_aligned_output)
        return normalized

    if is_aligned(normalized) and not cfg.force_mafft_alignment:
        LOGGER.info("[%s] Group already aligned; MAFFT skipped", group_label)
        write_fasta(normalized, group_aligned_output)
        return normalized

    mafft_exe = shutil.which(cfg.mafft_binary)
    if not mafft_exe:
        raise RuntimeError(
            "MAFFT was not found on PATH. Install it or provide aligned input."
        )

    LOGGER.info("[%s] Running MAFFT: %s", group_label, mafft_exe)
    write_fasta(normalized, group_tmp_input)
    cmd = [
        mafft_exe,
        "--thread",
        str(cfg.mafft_threads),
        "--auto",
        str(group_tmp_input),
    ]
    proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"MAFFT failed ({proc.returncode}): {proc.stderr}")

    group_aligned_output.write_text(proc.stdout, encoding="utf-8")
    aligned = read_fasta(group_aligned_output)
    LOGGER.info(
        "[%s] MAFFT complete; aligned output: %s",
        group_label,
        group_aligned_output,
    )
    return [
        SequenceRecord(header=r.header, sequence=normalize_sequence(r.sequence))
        for r in aligned
    ]


def consensus_from_alignment(
    aligned_records: list[SequenceRecord],
    priority: str,
) -> str:
    if not aligned_records:
        raise ValueError("Cannot build consensus from empty group")
    if not is_aligned(aligned_records):
        raise ValueError("Consensus requires aligned records")

    seq_len = len(aligned_records[0].sequence)
    consensus: list[str] = []

    for idx in range(seq_len):
        col = [rec.sequence[idx] for rec in aligned_records]
        bases = [b for b in col if b in {"A", "C", "G", "T"}]
        if not bases:
            consensus.append("-")
            continue

        counts = Counter(bases)
        top = max(counts.values())
        top_bases = {base for base, count in counts.items() if count == top}

        chosen = None
        for candidate in priority:
            if candidate in top_bases:
                chosen = candidate
                break
        consensus.append(chosen or "N")

    LOGGER.debug(
        "Consensus built from %d aligned sequences with length %d",
        len(aligned_records),
        seq_len,
    )
    return "".join(consensus)


def write_csv(rows: list[dict[str, str]], out_path: Path) -> None:
    if not rows:
        LOGGER.warning("Skipping CSV write for %s (no rows)", out_path)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    LOGGER.info("Wrote CSV: %s (%d rows)", out_path, len(rows))


def run_grouped_pipeline(cfg: GroupedConfig) -> dict[str, str | int]:
    LOGGER.info("Grouped pipeline started")
    LOGGER.info("Input FASTA: %s", cfg.input_fasta)
    LOGGER.info("Output directory: %s", cfg.output_dir)
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    bins_raw_dir = cfg.output_dir / "bins_raw"
    bins_aligned_dir = cfg.output_dir / "bins_aligned"
    bins_consensus_dir = cfg.output_dir / "consensus_per_group"
    for path in [bins_raw_dir, bins_aligned_dir, bins_consensus_dir]:
        path.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Output subdirectories ready")

    records = read_fasta(cfg.input_fasta)
    if not records:
        raise ValueError(f"No FASTA records found in {cfg.input_fasta}")

    if cfg.max_bin_size <= 0:
        raise ValueError("max_bin_size must be > 0")
    if cfg.unknown_year_handling not in {"separate", "discard"}:
        raise ValueError("unknown_year_handling must be 'separate' or 'discard'")

    year_bins, seq_map_rows = build_year_bins(records)

    final_bins: list[tuple[int | None, str, int, list[SequenceRecord], str]] = []
    seq_to_final_key: dict[str, str] = {}

    known_years = sorted([year for year in year_bins if year is not None])
    for year in known_years:
        year_records = year_bins.get(year, [])
        year_text = str(year)
        year_label = f"year_{year_text}"
        clusters = cluster_sequences_vsearch(year_records, cfg, year_label)
        clusters = merge_small_clusters(clusters, cfg.min_bin_size)
        year_candidate_bins: list[list[SequenceRecord]] = []
        for _cluster_raw_id, cluster_records in clusters:
            year_candidate_bins.extend(
                split_cluster_to_bins(cluster_records, cfg.max_bin_size)
            )
        year_final_bins = rebalance_bins_by_similarity(
            year_candidate_bins,
            min_size=cfg.min_bin_size,
            max_size=cfg.max_bin_size,
        )

        small_count = sum(1 for group in year_final_bins if len(group) < cfg.min_bin_size)
        if small_count:
            LOGGER.warning(
                "[%s] %d bins remain below min_bin_size=%d (likely due insufficient sequences)",
                year_label,
                small_count,
                cfg.min_bin_size,
            )

        LOGGER.info(
            "[%s] Found %d clusters and %d final bins for %d sequences",
            year_label,
            len(clusters),
            len(year_final_bins),
            len(year_records),
        )

        for bin_idx, bin_records in enumerate(year_final_bins, start=1):
            cluster_id = str(bin_idx)
            final_key = (
                f"SamplingYear={year_text}|ClusterId={cluster_id}|BinId=1"
            )
            for rec in bin_records:
                seq_to_final_key[rec.header] = final_key
            final_bins.append(
                (
                    year,
                    cluster_id,
                    1,
                    bin_records,
                    f"SamplingYear={year_text}",
                )
            )

    unknown_records = year_bins.get(None, [])
    if unknown_records:
        LOGGER.warning("Found %d sequences with UNKNOWN sampling year", len(unknown_records))
        if cfg.unknown_year_handling == "discard":
            LOGGER.warning("Discarding UNKNOWN year sequences by configuration")
        else:
            unknown_bins = split_cluster_to_bins(unknown_records, cfg.max_bin_size)
            unknown_bins = rebalance_bins_by_similarity(
                unknown_bins,
                min_size=cfg.min_bin_size,
                max_size=cfg.max_bin_size,
            )
            for bin_idx, bin_records in enumerate(unknown_bins, start=1):
                final_key = f"SamplingYear=UNKNOWN|ClusterId=UNKNOWN|BinId={bin_idx}"
                for rec in bin_records:
                    seq_to_final_key[rec.header] = final_key
                final_bins.append(
                    (
                        None,
                        "UNKNOWN",
                        bin_idx,
                        bin_records,
                        "SamplingYear=UNKNOWN",
                    )
                )

    LOGGER.info(
        "Proceeding with %d final bins from %d sequences",
        len(final_bins),
        len(records),
    )

    all_consensus_records: list[SequenceRecord] = []
    group_manifest_rows: list[dict[str, str]] = []

    for row in seq_map_rows:
        row["final_group_key"] = seq_to_final_key.get(row["header"], "")

    for group_index, (year, cluster_id, bin_id, group_records, source_key) in enumerate(final_bins, start=1):
        group_size = len(group_records)
        year_text = str(year) if year is not None else "UNKNOWN"
        source_keys = source_key

        group_tag = f"group_{group_index:03d}"
        final_group_key = f"SamplingYear={year_text}|ClusterId={cluster_id}|BinId={bin_id}"
        full_key = f"{final_group_key}|n={group_size}"
        consensus_header = f"Consensus|{group_tag}|{full_key}"

        raw_fasta_path = bins_raw_dir / f"{group_tag}.fasta"
        aligned_fasta_path = bins_aligned_dir / f"{group_tag}_aligned.fasta"
        consensus_fasta_path = bins_consensus_dir / f"{group_tag}_consensus.fasta"
        mafft_tmp_input = bins_aligned_dir / f"_{group_tag}_mafft_input.fasta"

        LOGGER.info(
            "[%s] Processing group key=%s with %d sequences",
            group_tag,
            final_group_key,
            group_size,
        )

        write_fasta(group_records, raw_fasta_path)
        LOGGER.info("[%s] Wrote raw group FASTA: %s", group_tag, raw_fasta_path)

        aligned_records = align_group_records(
            group_records=group_records,
            cfg=cfg,
            group_tmp_input=mafft_tmp_input,
            group_aligned_output=aligned_fasta_path,
            group_label=group_tag,
        )
        if mafft_tmp_input.exists():
            mafft_tmp_input.unlink()
            LOGGER.debug("[%s] Removed temporary MAFFT input", group_tag)

        consensus_seq = consensus_from_alignment(
            aligned_records,
            priority=cfg.consensus_tie_priority,
        )
        consensus_record = SequenceRecord(
            header=consensus_header,
            sequence=consensus_seq,
        )
        write_fasta([consensus_record], consensus_fasta_path)
        LOGGER.info(
            "[%s] Wrote consensus FASTA: %s",
            group_tag,
            consensus_fasta_path,
        )
        all_consensus_records.append(consensus_record)

        group_manifest_rows.append(
            {
                "group_id": group_tag,
                "full_key": full_key,
                "sampling_year": year_text,
                "cluster_id": str(cluster_id),
                "bin_id": str(bin_id),
                "hxb2_start": "VSEARCH",
                "n_sequences": str(group_size),
                "source_initial_keys": source_keys,
                "raw_fasta": str(raw_fasta_path),
                "aligned_fasta": str(aligned_fasta_path),
                "consensus_fasta": str(consensus_fasta_path),
            }
        )

    all_consensus_path = cfg.output_dir / "all_group_consensus.fasta"
    write_fasta(all_consensus_records, all_consensus_path)
    LOGGER.info(
        "Wrote combined consensus FASTA: %s (%d records)",
        all_consensus_path,
        len(all_consensus_records),
    )

    group_manifest_path = cfg.output_dir / "group_manifest.csv"
    write_csv(group_manifest_rows, group_manifest_path)

    seq_to_group_path = cfg.output_dir / "sequence_to_group_mapping.csv"
    write_csv(seq_map_rows, seq_to_group_path)

    run_manifest = {
        "input_fasta": str(cfg.input_fasta),
        "output_dir": str(cfg.output_dir),
        "min_bin_size": cfg.min_bin_size,
        "max_bin_size": cfg.max_bin_size,
        "vsearch_binary": cfg.vsearch_binary,
        "vsearch_identity": cfg.vsearch_identity,
        "unknown_year_handling": cfg.unknown_year_handling,
        "total_input_sequences": len(records),
        "n_final_groups": len(final_bins),
        "group_manifest_csv": str(group_manifest_path),
        "sequence_to_group_mapping_csv": str(seq_to_group_path),
        "all_group_consensus_fasta": str(all_consensus_path),
        "bins_raw_dir": str(bins_raw_dir),
        "bins_aligned_dir": str(bins_aligned_dir),
        "consensus_per_group_dir": str(bins_consensus_dir),
    }

    run_manifest_path = cfg.output_dir / "run_manifest.csv"
    with run_manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["key", "value"])
        for key, value in run_manifest.items():
            writer.writerow([key, value])
    LOGGER.info("Wrote run manifest: %s", run_manifest_path)

    LOGGER.info("Grouped pipeline complete")
    return run_manifest


def main() -> None:
    project_root = Path(__file__).resolve().parent
    cfg = GroupedConfig(
        input_fasta=(
            project_root
            / "data"
            / "hiv-db-any-unaligned.fasta"
        ),
        output_dir=project_root / "data" / "processed_grouped",
        min_bin_size=20,
        force_mafft_alignment=True,
        mafft_binary="mafft",
        mafft_threads=-1,
        vsearch_binary="vsearch",
        vsearch_identity=0.95,
        max_bin_size=30,
        unknown_year_handling="separate",
        consensus_tie_priority="ACGT",
        log_level="INFO",
    )

    configure_logging(cfg.log_level)
    results = run_grouped_pipeline(cfg)

    print("Grouped pipeline completed. Outputs:")
    for key, value in results.items():
        print(f"- {key}: {value}")


if __name__ == "__main__":
    main()

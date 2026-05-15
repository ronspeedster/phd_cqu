from __future__ import annotations

import argparse
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

DONOR_FIELD = "PAT id(SSAM)"
START_FIELD = "HXB2/MAC239 start"
STOP_FIELD = "HXB2/MAC239 stop"


@dataclass
class SequenceRecord:
    header: str
    sequence: str


@dataclass
class GroupedConfig:
    input_fasta: Path
    output_dir: Path
    donor_field: str = DONOR_FIELD
    coordinate_buffer: int = 50
    min_bin_size: int = 5  # Minimum sequences per group for meaningful consensus
    max_bin_size: int = 10 # Adviser requested max of 10
    force_mafft_alignment: bool = True
    mafft_binary: str = "mafft"
    mafft_threads: int = -1
    unknown_donor_handling: str = "separate"  # separate | discard
    consensus_tie_priority: str = "ACGT"
    max_sequences: int | None = None
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


def build_donor_bins(
    records: list[SequenceRecord],
    donor_field: str
) -> tuple[dict[str | None, list[SequenceRecord]], list[dict[str, str]]]:
    LOGGER.info("Building Level 1 bins using Donor ID: %s", donor_field)
    bins: dict[str | None, list[SequenceRecord]] = {}
    seq_map_rows: list[dict[str, str]] = []

    for rec in records:
        fields = parse_header_fields(rec.header)
        donor_id = fields.get(donor_field)
        
        # Treat empty strings as None
        if donor_id is not None and not donor_id.strip():
            donor_id = None
            
        bins.setdefault(donor_id, []).append(rec)
        
        donor_text = donor_id if donor_id is not None else "UNKNOWN"
        seq_map_rows.append(
            {
                "header": rec.header,
                "initial_donor_key": f"DonorId={donor_text}",
            }
        )

    unknown_count = len(bins.get(None, []))
    LOGGER.info(
        "Donor binning complete: %d unique donors (%d UNKNOWN sequences)",
        len(bins),
        unknown_count,
    )
    return bins, seq_map_rows


def cluster_by_region(
    records: list[SequenceRecord], 
    buffer: int
) -> tuple[list[tuple[list[SequenceRecord], int, int]], list[tuple[list[SequenceRecord], int, int]]]:
    """
    Clusters sequences based on HXB2/MAC239 start and stop coordinates.
    Two sequences belong to the same cluster if their buffered coordinate ranges overlap.
    Returns: (valid_clusters, missing_coordinate_clusters)
    Each cluster entry is (records, cluster_start, cluster_stop).
    """
    valid_recs = []
    missing_recs = []
    
    for r in records:
        fields = parse_header_fields(r.header)
        start = safe_int(fields.get(START_FIELD))
        stop = safe_int(fields.get(STOP_FIELD))
        
        if start is not None and stop is not None:
            actual_start = min(start, stop)
            actual_stop = max(start, stop)
            valid_recs.append((actual_start, actual_stop, r))
        else:
            missing_recs.append(r)
            
    valid_recs.sort(key=lambda x: x[0])
    
    clusters: list[tuple[list[SequenceRecord], int, int]] = []
    
    if valid_recs:
        current_cluster = [valid_recs[0][2]]
        current_start = valid_recs[0][0] - buffer
        current_stop = valid_recs[0][1] + buffer
        cluster_real_start = valid_recs[0][0]
        cluster_real_stop = valid_recs[0][1]
        
        for start, stop, rec in valid_recs[1:]:
            rec_start = start - buffer
            rec_stop = stop + buffer
            
            if rec_start <= current_stop:
                current_cluster.append(rec)
                current_stop = max(current_stop, rec_stop)
                cluster_real_stop = max(cluster_real_stop, stop)
            else:
                clusters.append((current_cluster, cluster_real_start, cluster_real_stop))
                current_cluster = [rec]
                current_start = rec_start
                current_stop = rec_stop
                cluster_real_start = start
                cluster_real_stop = stop
                
        clusters.append((current_cluster, cluster_real_start, cluster_real_stop))

    missing_clusters: list[tuple[list[SequenceRecord], int, int]] = []
    if missing_recs:
        missing_clusters.append((missing_recs, 0, 0))

    return clusters, missing_clusters


def merge_small_region_clusters(
    clusters: list[tuple[list[SequenceRecord], int, int]],
    min_size: int,
) -> list[tuple[list[SequenceRecord], int, int]]:
    """
    Merges region clusters with fewer than min_size sequences into the nearest
    cluster (by coordinate center distance). If all clusters are small, merges
    them all into one.
    """
    if not clusters:
        return clusters

    if len(clusters) == 1:
        return clusters

    records_list, starts, stops = [], [], []
    for recs, s, e in clusters:
        records_list.append(recs)
        starts.append(s)
        stops.append(e)

    centers = [(s + e) / 2.0 for s, e in zip(starts, stops)]
    sizes = [len(recs) for recs in records_list]
    merged_into: dict[int, int] = {}

    small_indices = [i for i, sz in enumerate(sizes) if sz < min_size]
    large_indices = [i for i, sz in enumerate(sizes) if sz >= min_size]

    if not large_indices:
        all_records: list[SequenceRecord] = []
        all_start = min(starts)
        all_stop = max(stops)
        for recs in records_list:
            all_records.extend(recs)
        return [(all_records, all_start, all_stop)]

    for si in small_indices:
        best_target = min(large_indices, key=lambda li: abs(centers[si] - centers[li]))
        merged_into[si] = best_target

    result: list[tuple[list[SequenceRecord], int, int]] = []
    bucket_records: dict[int, list[SequenceRecord]] = {}
    bucket_start: dict[int, int] = {}
    bucket_stop: dict[int, int] = {}

    for i in range(len(clusters)):
        if i in merged_into:
            target = merged_into[i]
            bucket_records.setdefault(target, []).extend(records_list[i])
            bucket_start[target] = min(bucket_start.get(target, starts[i]), starts[i])
            bucket_stop[target] = max(bucket_stop.get(target, stops[i]), stops[i])
        else:
            bucket_records.setdefault(i, []).extend(records_list[i])
            bucket_start[i] = min(bucket_start.get(i, starts[i]), starts[i])
            bucket_stop[i] = max(bucket_stop.get(i, stops[i]), stops[i])

    for i in sorted(bucket_records.keys()):
        result.append((bucket_records[i], bucket_start[i], bucket_stop[i]))

    LOGGER.info(
        "Merged %d small clusters (min_size=%d) into %d remaining clusters",
        len(small_indices),
        min_size,
        len(result),
    )
    return result


def split_cluster_to_bins(
    cluster_records: list[SequenceRecord],
    max_size: int,
    min_size: int = 3,
) -> list[list[SequenceRecord]]:
    if len(cluster_records) <= max_size:
        return [cluster_records]

    bins: list[list[SequenceRecord]] = []
    for idx in range(0, len(cluster_records), max_size):
        bins.append(cluster_records[idx : idx + max_size])

    while len(bins) >= 2 and len(bins[-1]) < min_size:
        shortfall = min_size - len(bins[-1])
        available = len(bins[-2])
        if available - shortfall >= min_size:
            transfer = bins[-2][-shortfall:]
            bins[-2] = bins[-2][:-shortfall]
            bins[-1] = transfer + bins[-1]
            break
        else:
            bins[-2].extend(bins[-1])
            bins.pop()

    if len(bins) >= 2 and len(bins[-1]) < min_size:
        bins[-2].extend(bins[-1])
        bins.pop()

    return bins


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

    if cfg.max_sequences is not None and cfg.max_sequences > 0:
        records = records[: cfg.max_sequences]
        LOGGER.info("Limited to first %d sequences", len(records))

    if cfg.max_bin_size <= 0:
        raise ValueError("max_bin_size must be > 0")

    # LEVEL 1: Group by Donor ID
    donor_bins, seq_map_rows = build_donor_bins(records, cfg.donor_field)

    final_bins: list[tuple[str | None, int, int, list[SequenceRecord], str]] = []
    seq_to_final_key: dict[str, str] = {}

    known_donors = sorted([d for d in donor_bins if d is not None])
    
    for donor in known_donors:
        donor_records = donor_bins.get(donor, [])
        donor_label = str(donor)
        
        # LEVEL 2: Cluster by Genomic Region Overlap (+/- 50 buffer)
        valid_region_clusters, missing_coords_clusters = cluster_by_region(
            donor_records, 
            buffer=cfg.coordinate_buffer
        )
        
        all_donor_clusters_raw = valid_region_clusters + missing_coords_clusters
        
        LOGGER.info(
            "[Donor %s] Found %d regional clusters for %d sequences",
            donor_label,
            len(all_donor_clusters_raw),
            len(donor_records),
        )

        all_donor_clusters = merge_small_region_clusters(
            all_donor_clusters_raw, cfg.min_bin_size
        )

        # LEVEL 3: Enforce max_bin_size chunks (10 sequences)
        for region_idx, (cluster_records, _cstart, _cstop) in enumerate(all_donor_clusters, start=1):
            chunked_bins = split_cluster_to_bins(cluster_records, cfg.max_bin_size, cfg.min_bin_size)
            
            for bin_idx, bin_records in enumerate(chunked_bins, start=1):
                final_key = f"Donor={donor_label}|RegionId={region_idx}|BinId={bin_idx}"
                
                for rec in bin_records:
                    seq_to_final_key[rec.header] = final_key
                    
                final_bins.append(
                    (
                        donor,
                        region_idx,
                        bin_idx,
                        bin_records,
                        f"Donor={donor_label}",
                    )
                )

    unknown_records = donor_bins.get(None, [])
    if unknown_records:
        LOGGER.warning("Found %d sequences with UNKNOWN Donor", len(unknown_records))
        if cfg.unknown_donor_handling == "discard":
            LOGGER.warning("Discarding UNKNOWN donor sequences by configuration")
        else:
            valid_region_clusters, missing_coords_clusters = cluster_by_region(
                unknown_records, 
                buffer=cfg.coordinate_buffer
            )
            all_unknown_clusters_raw = valid_region_clusters + missing_coords_clusters

            all_unknown_clusters = merge_small_region_clusters(
                all_unknown_clusters_raw, cfg.min_bin_size
            )
            
            for region_idx, (cluster_records, _cstart, _cstop) in enumerate(all_unknown_clusters, start=1):
                chunked_bins = split_cluster_to_bins(cluster_records, cfg.max_bin_size, cfg.min_bin_size)
                
                for bin_idx, bin_records in enumerate(chunked_bins, start=1):
                    final_key = f"Donor=UNKNOWN|RegionId={region_idx}|BinId={bin_idx}"
                    for rec in bin_records:
                        seq_to_final_key[rec.header] = final_key
                        
                    final_bins.append(
                        (
                            None,
                            region_idx,
                            bin_idx,
                            bin_records,
                            "Donor=UNKNOWN",
                        )
                    )

    # CROSS-DONOR MERGE: collect undersized groups and re-cluster by region
    valid_bins = []
    undersized_records: list[SequenceRecord] = []
    for entry in final_bins:
        _donor, _rid, _bid, group_recs, _src = entry
        if len(group_recs) < cfg.min_bin_size:
            undersized_records.extend(group_recs)
        else:
            valid_bins.append(entry)

    if undersized_records:
        LOGGER.info(
            "Cross-donor merge: %d undersized sequences from %d groups, re-clustering by region",
            len(undersized_records),
            len(final_bins) - len(valid_bins),
        )
        cross_valid, cross_missing = cluster_by_region(
            undersized_records, buffer=cfg.coordinate_buffer
        )
        cross_all = cross_valid + cross_missing
        cross_merged = merge_small_region_clusters(cross_all, cfg.min_bin_size)

        for region_idx, (cluster_records, _cs, _ce) in enumerate(cross_merged, start=1):
            chunked_bins = split_cluster_to_bins(cluster_records, cfg.max_bin_size, cfg.min_bin_size)
            for bin_idx, bin_records in enumerate(chunked_bins, start=1):
                for rec in bin_records:
                    seq_to_final_key[rec.header] = f"Donor=CROSSDONOR|RegionId={region_idx}|BinId={bin_idx}"
                valid_bins.append(
                    (None, region_idx, bin_idx, bin_records, "Donor=CROSSDONOR")
                )

        LOGGER.info(
            "Cross-donor merge complete: %d final groups (%d undersized sequences absorbed)",
            len(valid_bins),
            len(undersized_records),
        )
    final_bins = valid_bins

    LOGGER.info(
        "Proceeding with %d final bins from %d sequences",
        len(final_bins),
        len(records),
    )

    all_consensus_records: list[SequenceRecord] = []
    group_manifest_rows: list[dict[str, str]] = []

    for row in seq_map_rows:
        row["final_group_key"] = seq_to_final_key.get(row["header"], "")

    for group_index, (donor, region_id, bin_id, group_records, source_key) in enumerate(final_bins, start=1):
        group_size = len(group_records)
        donor_text = "CROSSDONOR" if source_key == "Donor=CROSSDONOR" else (str(donor) if donor is not None else "UNKNOWN")

        group_tag = f"group_{group_index:03d}"
        final_group_key = f"Donor={donor_text}|RegionId={region_id}|BinId={bin_id}"
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

        aligned_records = align_group_records(
            group_records=group_records,
            cfg=cfg,
            group_tmp_input=mafft_tmp_input,
            group_aligned_output=aligned_fasta_path,
            group_label=group_tag,
        )
        if mafft_tmp_input.exists():
            mafft_tmp_input.unlink()

        consensus_seq = consensus_from_alignment(
            aligned_records,
            priority=cfg.consensus_tie_priority,
        )
        consensus_record = SequenceRecord(
            header=consensus_header,
            sequence=consensus_seq,
        )
        write_fasta([consensus_record], consensus_fasta_path)
        all_consensus_records.append(consensus_record)

        group_manifest_rows.append(
            {
                "group_id": group_tag,
                "full_key": full_key,
                "donor_id": donor_text,
                "region_cluster_id": str(region_id),
                "bin_id": str(bin_id),
                "n_sequences": str(group_size),
                "raw_fasta": str(raw_fasta_path),
                "aligned_fasta": str(aligned_fasta_path),
                "consensus_fasta": str(consensus_fasta_path),
            }
        )

    all_consensus_path = cfg.output_dir / "all_group_consensus.fasta"
    write_fasta(all_consensus_records, all_consensus_path)

    group_manifest_path = cfg.output_dir / "group_manifest.csv"
    write_csv(group_manifest_rows, group_manifest_path)

    seq_to_group_path = cfg.output_dir / "sequence_to_group_mapping.csv"
    write_csv(seq_map_rows, seq_to_group_path)

    run_manifest = {
        "input_fasta": str(cfg.input_fasta),
        "output_dir": str(cfg.output_dir),
        "donor_field": cfg.donor_field,
        "coordinate_buffer": cfg.coordinate_buffer,
        "min_bin_size": cfg.min_bin_size,
        "max_bin_size": cfg.max_bin_size,
        "unknown_donor_handling": cfg.unknown_donor_handling,
        "total_input_sequences": len(records),
        "n_final_groups": len(final_bins),
        "group_manifest_csv": str(group_manifest_path),
        "sequence_to_group_mapping_csv": str(seq_to_group_path),
        "all_group_consensus_fasta": str(all_consensus_path),
    }

    run_manifest_path = cfg.output_dir / "run_manifest.csv"
    with run_manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["key", "value"])
        for key, value in run_manifest.items():
            writer.writerow([key, value])

    LOGGER.info("Grouped pipeline complete")
    return run_manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Donor + region clustering pipeline for HIV FASTA data."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Path to input FASTA file."
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path to output directory."
    )
    parser.add_argument(
        "-n",
        "--limit",
        type=int,
        default=None,
        help="Process only the first N sequences from the input FASTA.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    cfg = GroupedConfig(
        input_fasta=Path(args.input),
        output_dir=Path(args.output),
        donor_field=DONOR_FIELD,
        coordinate_buffer=50,
        min_bin_size=5,
        max_bin_size=10,
        force_mafft_alignment=True,
        mafft_binary="mafft",
        mafft_threads=-1,
        unknown_donor_handling="separate",
        consensus_tie_priority="ACGT",
        max_sequences=args.limit,
        log_level="INFO",
    )

    configure_logging(cfg.log_level)
    results = run_grouped_pipeline(cfg)

    print("Grouped pipeline completed. Outputs:")
    for key, value in results.items():
        print(f"- {key}: {value}")


# Usage:
# python 01_data_preparation_donor_region_clustering.py -i data/hiv-db-any-unaligned.fasta -o data/processed_grouped
# python 01_data_preparation_donor_region_clustering.py -i data/hiv-db-any-unaligned.fasta -o data/processed_grouped -n 10


if __name__ == "__main__":
    main()
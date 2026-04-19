from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import logging
import shutil
import subprocess
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
    min_bin_size: int = 5
    force_mafft_alignment: bool = False
    mafft_binary: str = "mafft"
    mafft_threads: int = -1
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

    initial_bins, seq_map_rows = build_initial_bins(records)
    final_keys, merge_trace = merge_small_known_bins(initial_bins, cfg.min_bin_size)
    LOGGER.info(
        "Proceeding with %d final groups from %d sequences",
        len(final_keys),
        len(records),
    )

    all_consensus_records: list[SequenceRecord] = []
    group_manifest_rows: list[dict[str, str]] = []

    initial_to_final: dict[str, str] = {}
    for final_key in final_keys:
        final_label = key_to_text(final_key)
        for source_label in merge_trace.get(final_key, [final_label]):
            initial_to_final[source_label] = final_label

    for row in seq_map_rows:
        src = row["initial_group_key"]
        row["final_group_key"] = initial_to_final.get(src, src)

    for group_index, key in enumerate(final_keys, start=1):
        group_records = initial_bins[key]
        group_size = len(group_records)
        year, start = key
        year_text = str(year) if year is not None else "UNKNOWN"
        start_text = str(start) if start is not None else "UNKNOWN"
        source_keys = ";".join(merge_trace.get(key, [key_to_text(key)]))

        group_tag = f"group_{group_index:03d}"
        full_key = (
            f"SamplingYear={year_text}|HXB2Start={start_text}|n={group_size}"
        )
        consensus_header = f"Consensus|{group_tag}|{full_key}"

        raw_fasta_path = bins_raw_dir / f"{group_tag}.fasta"
        aligned_fasta_path = bins_aligned_dir / f"{group_tag}_aligned.fasta"
        consensus_fasta_path = bins_consensus_dir / f"{group_tag}_consensus.fasta"
        mafft_tmp_input = bins_aligned_dir / f"_{group_tag}_mafft_input.fasta"

        LOGGER.info(
            "[%s] Processing group key=%s with %d sequences",
            group_tag,
            key_to_text(key),
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
                "hxb2_start": start_text,
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
        "total_input_sequences": len(records),
        "n_final_groups": len(final_keys),
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
        min_bin_size=5,
        force_mafft_alignment=False,
        mafft_binary="mafft",
        mafft_threads=-1,
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

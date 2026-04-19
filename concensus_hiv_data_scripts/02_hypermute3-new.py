from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
from decimal import Decimal, InvalidOperation
import logging
import re
import shutil
import subprocess
import sys


LOGGER = logging.getLogger("hypermute3_pairwise")
DNA_STRICT_ALLOWED = set("ACGT-")


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
    log_level: str = "INFO"


def configure_logging(log_level: str) -> None:
    level_name = (log_level or "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


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


def sanitize_reference_for_partial(
    input_fasta: Path,
    output_fasta: Path,
) -> None:
    records = read_fasta_records(input_fasta)
    if not records:
        raise ValueError(f"No records found in {input_fasta}")

    ref_header, ref_seq = records[0]
    records[0] = (ref_header, normalize_sequence(ref_seq))
    write_fasta_records(records, output_fasta)


def run_hypermut_once(
    cfg: PipelineConfig,
    run_cfg: HypermutRunConfig,
    aligned_pair_fasta: Path,
    seq_output_dir: Path,
) -> dict[str, str]:
    prefix = str(seq_output_dir / f"{run_cfg.name}-")
    run_input_fasta = aligned_pair_fasta

    if run_cfg.match == "partial":
        run_input_fasta = seq_output_dir / f"{run_cfg.name}-input-sanitized.fasta"
        sanitize_reference_for_partial(aligned_pair_fasta, run_input_fasta)

    cmd = [
        cfg.python_executable,
        str(cfg.hypermut_script),
        str(run_input_fasta),
        cfg.mutation_from,
        cfg.mutation_to,
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

    summary_csv = seq_output_dir / f"{run_cfg.name}-summary.csv"
    positions_csv = seq_output_dir / f"{run_cfg.name}-positions.csv"
    args_csv = seq_output_dir / f"{run_cfg.name}-args.csv"

    return {
        "run_name": run_cfg.name,
        "match": run_cfg.match,
        "keepgaps": str(run_cfg.keepgaps),
        "summary_csv": str(summary_csv),
        "positions_csv": str(positions_csv),
        "args_csv": str(args_csv),
    }


def read_single_hypermut_summary(summary_csv: Path) -> dict[str, str]:
    if not summary_csv.exists():
        return {
            "seq_name": "",
            "primary_matches": "",
            "potential_primaries": "",
            "control_matches": "",
            "potential_controls": "",
            "rate_ratio": "",
            "fisher_p": "",
        }

    with summary_csv.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        first = next(reader, None)
        if first is None:
            return {
                "seq_name": "",
                "primary_matches": "",
                "potential_primaries": "",
                "control_matches": "",
                "potential_controls": "",
                "rate_ratio": "",
                "fisher_p": "",
            }
        return {
            "seq_name": first.get("seq_name", ""),
            "primary_matches": to_non_scientific(first.get("primary_matches", "")),
            "potential_primaries": to_non_scientific(
                first.get("potential_primaries", "")
            ),
            "control_matches": to_non_scientific(first.get("control_matches", "")),
            "potential_controls": to_non_scientific(
                first.get("potential_controls", "")
            ),
            "rate_ratio": to_non_scientific(first.get("rate_ratio", "")),
            "fisher_p": to_non_scientific(first.get("fisher_p", "")),
        }


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
            reduced_key = f"SamplingYear={sampling_year}|HXB2Start={hxb2_start}"
            consensus_fasta = Path(row.get("consensus_fasta", "").strip())
            if not consensus_fasta.exists():
                continue
            records = read_fasta_records(consensus_fasta)
            if not records:
                continue

            # Key used by sequence_to_group_mapping.csv
            out[reduced_key] = (records[0][0], records[0][1])
            # Keep compatibility with full-key lookups when available.
            if full_key:
                out[full_key] = (records[0][0], records[0][1])

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
    per_seq_dir = cfg.output_dir / "per_sequence_hypermut"
    for path in [pair_input_dir, pair_aligned_dir, per_seq_dir]:
        path.mkdir(parents=True, exist_ok=True)

    group_consensus = load_group_consensus_map(cfg.group_manifest_csv)
    seq_to_group = load_sequence_to_group_map(cfg.sequence_to_group_csv)
    LOGGER.info(
        "Loaded %d group consensus entries and %d sequence mappings",
        len(group_consensus),
        len(seq_to_group),
    )

    run_modes = [
        HypermutRunConfig(name="strict-keepgaps", match="strict", keepgaps=True),
        HypermutRunConfig(name="strict-skipgaps", match="strict", keepgaps=False),
        HypermutRunConfig(name="partial-keepgaps", match="partial", keepgaps=True),
        HypermutRunConfig(name="partial-skipgaps", match="partial", keepgaps=False),
    ]

    sequence_records: list[tuple[str, str, str]] = []
    for seq_header, final_key in seq_to_group.items():
        sequence_records.append((seq_header, "", final_key))

    # We reconstruct sequences from group raw FASTA files listed in group manifest.
    with cfg.group_manifest_csv.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        raw_map: dict[str, str] = {}
        for row in reader:
            full_key = row.get("full_key", "").strip()
            raw_fasta = Path(row.get("raw_fasta", "").strip())
            if not full_key or not raw_fasta.exists():
                continue
            for h, s in read_fasta_records(raw_fasta):
                raw_map[h] = s

    rows_merged: list[dict[str, str]] = []
    rows_failures: list[dict[str, str]] = []

    total = len(sequence_records)
    LOGGER.info("Starting per-sequence pairwise alignment + Hypermut3 for %d sequences", total)

    for idx, (seq_header, _, final_key) in enumerate(sequence_records, start=1):
        seq_body = raw_map.get(seq_header)
        if seq_body is None:
            rows_failures.append(
                {
                    "seq_header": seq_header,
                    "final_group_key": final_key,
                    "error": "Sequence not found in raw group FASTA map",
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
        seq_slug = f"{idx:04d}_{sanitize_file_component(seq_header)}"
        pair_input = pair_input_dir / f"{seq_slug}.fasta"
        pair_aligned = pair_aligned_dir / f"{seq_slug}_aligned.fasta"
        seq_run_dir = per_seq_dir / seq_slug
        seq_run_dir.mkdir(parents=True, exist_ok=True)

        if idx == 1 or idx % 25 == 0 or idx == total:
            LOGGER.info("Processing sequence %d/%d: %s", idx, total, seq_header)

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

            for mode in run_modes:
                run_info = run_hypermut_once(
                    cfg=cfg,
                    run_cfg=mode,
                    aligned_pair_fasta=pair_aligned,
                    seq_output_dir=seq_run_dir,
                )
                sum_row = read_single_hypermut_summary(Path(run_info["summary_csv"]))
                rows_merged.append(
                    {
                        "seq_header": seq_header,
                        "final_group_key": final_key,
                        "consensus_header": consensus_header,
                        "pair_aligned_fasta": str(pair_aligned),
                        "run_name": run_info["run_name"],
                        "match": run_info["match"],
                        "keepgaps": run_info["keepgaps"],
                        "seq_name": sum_row["seq_name"],
                        "primary_matches": sum_row["primary_matches"],
                        "potential_primaries": sum_row["potential_primaries"],
                        "control_matches": sum_row["control_matches"],
                        "potential_controls": sum_row["potential_controls"],
                        "rate_ratio": sum_row["rate_ratio"],
                        "fisher_p": sum_row["fisher_p"],
                        "summary_csv": run_info["summary_csv"],
                        "positions_csv": run_info["positions_csv"],
                    }
                )

        except Exception as exc:
            LOGGER.exception("Failed sequence: %s", seq_header)
            rows_failures.append(
                {
                    "seq_header": seq_header,
                    "final_group_key": final_key,
                    "error": str(exc),
                }
            )

    merged_csv = cfg.output_dir / "per_sequence_hypermut_merged.csv"
    failed_csv = cfg.output_dir / "per_sequence_hypermut_failures.csv"
    write_csv(rows_merged, merged_csv)
    write_csv(rows_failures, failed_csv)

    LOGGER.info("Pairwise Hypermut3 pipeline finished")
    LOGGER.info("Success rows: %d", len(rows_merged))
    LOGGER.info("Failure rows: %d", len(rows_failures))

    return {
        "group_manifest_csv": str(cfg.group_manifest_csv),
        "sequence_to_group_csv": str(cfg.sequence_to_group_csv),
        "output_dir": str(cfg.output_dir),
        "merged_csv": str(merged_csv),
        "failures_csv": str(failed_csv),
        "n_sequences": total,
        "n_success_rows": len(rows_merged),
        "n_failure_rows": len(rows_failures),
    }


def main() -> None:
    project_root = Path(__file__).resolve().parent

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
        log_level="INFO",
    )

    configure_logging(cfg.log_level)
    results = run_pipeline(cfg)

    print("Pairwise Hypermut3 completed. Key outputs:")
    for key, value in results.items():
        print(f"- {key}: {value}")


if __name__ == "__main__":
    main()

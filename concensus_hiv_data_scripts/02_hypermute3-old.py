# section 1

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
from decimal import Decimal, InvalidOperation
import logging
import subprocess
import sys
from typing import List


# section 2

LOGGER = logging.getLogger("hypermut3_runner")
DNA_STRICT_ALLOWED = set("ACGT-")


@dataclass
class HypermutRunConfig:
    name: str
    match: str  # strict | partial
    keepgaps: bool


@dataclass
class PipelineConfig:
    input_fasta: Path
    hypermut_script: Path
    output_dir: Path
    python_executable: str = sys.executable
    mutation_from: str = "G"
    mutation_to: str = "A"
    upstream_context: str = ""
    downstream_context: str = "RD"
    enforce: str = "D"
    begin: int = 0
    finish: int | None = None
    log_level: str = "INFO"


# section 3

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
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def sanitize_reference_for_partial(
    input_fasta: Path,
    output_fasta: Path,
) -> Path:
    records = read_fasta_records(input_fasta)
    if not records:
        raise ValueError(f"No records found in {input_fasta}")

    ref_header, ref_seq = records[0]
    sanitized_ref = "".join(
        base if base in DNA_STRICT_ALLOWED else "-"
        for base in ref_seq.upper().replace("U", "T")
    )
    records[0] = (ref_header, sanitized_ref)
    write_fasta_records(records, output_fasta)
    return output_fasta


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


def run_hypermut_once(cfg: PipelineConfig, run_cfg: HypermutRunConfig) -> dict:
    prefix = str(cfg.output_dir / f"{run_cfg.name}-")
    run_input_fasta = cfg.input_fasta

    if run_cfg.match == "partial":
        run_input_fasta = (
            cfg.output_dir / f"{run_cfg.name}-input-sanitized.fasta"
        )
        sanitize_reference_for_partial(cfg.input_fasta, run_input_fasta)
        LOGGER.info(
            "Partial mode input sanitized for A/C/G/T/- reference: %s",
            run_input_fasta,
        )

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

    LOGGER.info("Running Hypermut3 mode: %s", run_cfg.name)
    proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            "Hypermut3 failed in mode "
            f"'{run_cfg.name}' (code {proc.returncode})\n"
            f"STDERR:\n{proc.stderr}\nSTDOUT:\n{proc.stdout}"
        )

    summary_csv = cfg.output_dir / f"{run_cfg.name}-summary.csv"
    positions_csv = cfg.output_dir / f"{run_cfg.name}-positions.csv"
    args_csv = cfg.output_dir / f"{run_cfg.name}-args.csv"

    LOGGER.info("Completed mode: %s", run_cfg.name)
    return {
        "run_name": run_cfg.name,
        "match": run_cfg.match,
        "keepgaps": run_cfg.keepgaps,
        "summary_csv": str(summary_csv),
        "positions_csv": str(positions_csv),
        "args_csv": str(args_csv),
    }


def aggregate_summaries(run_outputs: List[dict], out_file: Path) -> None:
    merged_rows = []
    for run in run_outputs:
        summary_path = Path(run["summary_csv"])
        if not summary_path.exists():
            LOGGER.warning(
                "Missing summary file for %s: %s",
                run["run_name"],
                summary_path,
            )
            continue

        with summary_path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                merged_rows.append(
                    {
                        "run_name": run["run_name"],
                        "match": run["match"],
                        "keepgaps": run["keepgaps"],
                        "seq_name": row.get("seq_name", ""),
                        "primary_matches": to_non_scientific(
                            row.get("primary_matches", "")
                        ),
                        "potential_primaries": to_non_scientific(
                            row.get("potential_primaries", "")
                        ),
                        "control_matches": to_non_scientific(
                            row.get("control_matches", "")
                        ),
                        "potential_controls": to_non_scientific(
                            row.get("potential_controls", "")
                        ),
                        "rate_ratio": to_non_scientific(
                            row.get("rate_ratio", "")
                        ),
                        "fisher_p": to_non_scientific(
                            row.get("fisher_p", "")
                        ),
                    }
                )

    if not merged_rows:
        LOGGER.warning("No summary rows found to aggregate")
        return

    merged_rows.sort(
        key=lambda row: (
            row.get("seq_name", ""),
            row.get("run_name", ""),
        )
    )

    fields = list(merged_rows[0].keys())
    with out_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(merged_rows)

    LOGGER.info("Aggregated summary written to %s", out_file)


# section 4

def run_pipeline(cfg: PipelineConfig) -> dict:
    LOGGER.info("Hypermut3 pipeline started")

    if not cfg.input_fasta.exists():
        raise FileNotFoundError(f"Input FASTA not found: {cfg.input_fasta}")
    if not cfg.hypermut_script.exists():
        raise FileNotFoundError(
            f"Hypermut3 script not found: {cfg.hypermut_script}"
        )

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Output directory: %s", cfg.output_dir)

    run_modes = [
        HypermutRunConfig(
            name="strict-keepgaps",
            match="strict",
            keepgaps=True,
        ),
        HypermutRunConfig(
            name="strict-skipgaps",
            match="strict",
            keepgaps=False,
        ),
        HypermutRunConfig(
            name="partial-keepgaps",
            match="partial",
            keepgaps=True,
        ),
        HypermutRunConfig(
            name="partial-skipgaps",
            match="partial",
            keepgaps=False,
        ),
    ]

    run_outputs = []
    for mode in run_modes:
        run_outputs.append(run_hypermut_once(cfg, mode))

    aggregate_path = cfg.output_dir / "all_runs_summary_merged.csv"
    aggregate_summaries(run_outputs, aggregate_path)

    manifest_path = cfg.output_dir / "hypermut3_runs_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(run_outputs[0].keys()))
        writer.writeheader()
        writer.writerows(run_outputs)

    LOGGER.info("Runs manifest written to %s", manifest_path)
    LOGGER.info("Hypermut3 pipeline finished successfully")

    return {
        "input_fasta": str(cfg.input_fasta),
        "manifest_csv": str(manifest_path),
        "merged_summary_csv": str(aggregate_path),
        "n_modes_run": len(run_outputs),
    }


def main() -> None:
    project_root = Path(__file__).resolve().parent

    cfg = PipelineConfig(
        input_fasta=project_root
        / "data"
        / "processed"
        / "hypermut3_input_with_consensus_first.fasta",
        hypermut_script=project_root / "Hypermut3" / "hypermut.py",
        output_dir=project_root / "data" / "hypermut3_results",
        python_executable=sys.executable,
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

    print("Hypermut3 run completed. Key outputs:")
    for key, value in results.items():
        print(f"- {key}: {value}")


if __name__ == "__main__":
    main()

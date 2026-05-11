from __future__ import annotations

import argparse
import logging
import random
import string
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Literal, Optional, Tuple

LOGGER = logging.getLogger("fasta_file_generator")

DEFAULT_ALPHABET = "ACGT"
DEFAULT_LINE_WIDTH = 80

HEADER_FIELD_NAMES: Tuple[str, ...] = (
    "se_id",
    "patient_code",
    "pat_id_ssam",
    "accession",
    "name",
    "subtype",
    "country",
    "sampling_year",
    "problematic_sequence",
    "hxb2_start",
    "hxb2_stop",
    "sequence_length",
    "organism",
)

COUNTRY_CODES: Tuple[str, ...] = (
    "CI", "US", "KE", "UG", "TZ", "ZA", "NG", "CM", "CD", "BR",
    "IN", "TH", "CN", "GB", "DE", "FR", "SE", "AU", "JP", "MX",
)

SUBTYPES: Tuple[str, ...] = ("A", "B", "C", "D", "F", "G", "H", "J", "K", "CRF01_AE", "CRF02_AG")

HXB2_REGIONS: Tuple[Tuple[int, int], ...] = (
    (790, 2292),
    (2253, 2549),
    (4483, 5568),
    (6062, 7050),
    (7050, 7400),
    (5831, 8469),
    (4550, 5850),
)


@dataclass
class LANLHeaderTemplate:
    se_id_start: int = 150000
    patient_code: str = "-"
    pat_id_ssam: str = "-"
    accession_prefix: str = "AF"
    accession_start: int = 0
    name_prefix: str = "IC"
    name_start: int = 1000
    subtype: str = "A"
    country: str = "CI"
    sampling_year: str = "1995"
    problematic_sequence: int = 0
    hxb2_start: int = 7050
    hxb2_stop: int = 7400
    organism: str = "HIV-1"
    field_separator: str = ","
    missing_char: str = "-"


@dataclass
class SequenceRecord:
    header: str
    sequence: str


# section 2 — header construction


def _build_header(
    template: LANLHeaderTemplate,
    index: int,
    seq_length: int,
) -> str:
    se_id = template.se_id_start + index
    accession_num = template.accession_start + index
    accession = f"{template.accession_prefix}{accession_num:06d}"
    name = f"{template.name_prefix}{template.name_start + index}"

    fields = [
        str(se_id),
        template.patient_code,
        template.pat_id_ssam,
        accession,
        name,
        template.subtype,
        template.country,
        template.sampling_year,
        str(template.problematic_sequence),
        str(template.hxb2_start),
        str(template.hxb2_stop),
        str(seq_length),
        template.organism,
    ]

    assert len(fields) == len(HEADER_FIELD_NAMES), (
        f"Field count mismatch: {len(fields)} vs {len(HEADER_FIELD_NAMES)}"
    )

    return template.field_separator.join(fields)


# section 3 — sequence generation


def _generate_sequence(
    length: int,
    alphabet: str,
    rng: random.Random,
) -> str:
    return "".join(rng.choice(alphabet) for _ in range(length))


def _mutate_sequence(
    parent_seq: str,
    mutation_rate: float,
    alphabet: str,
    rng: random.Random,
) -> str:
    seq = list(parent_seq)
    num_mutations = max(1, int(len(seq) * mutation_rate))
    positions = rng.sample(range(len(seq)), num_mutations)
    for pos in positions:
        original = seq[pos]
        candidates = [b for b in alphabet if b != original]
        seq[pos] = rng.choice(candidates)
    return "".join(seq)


def _random_subsequence(parent_seq: str, length: int, rng: random.Random) -> str:
    max_start = len(parent_seq) - length
    start = rng.randint(0, max_start)
    return parent_seq[start : start + length]


# section 4 — FASTA I/O


def write_fasta(
    records: List[SequenceRecord],
    file_path: Path,
    line_width: int = DEFAULT_LINE_WIDTH,
) -> None:
    with file_path.open("w", encoding="utf-8") as handle:
        for rec in records:
            handle.write(f">{rec.header}\n")
            seq = rec.sequence
            for i in range(0, len(seq), line_width):
                handle.write(seq[i : i + line_width] + "\n")


def read_fasta(file_path: Path) -> List[SequenceRecord]:
    records: List[SequenceRecord] = []
    header = None
    seq_parts: List[str] = []

    with file_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(SequenceRecord(header=header, sequence="".join(seq_parts)))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        records.append(SequenceRecord(header=header, sequence="".join(seq_parts)))

    return records


# section 5 — template parsing


def parse_template_fasta(fasta_path: Path) -> LANLHeaderTemplate:
    records = read_fasta(fasta_path)
    if not records:
        raise ValueError(f"No records found in {fasta_path}")

    first_header = records[0].header

    sep = ","
    if "\t" in first_header:
        sep = "\t"

    fields = first_header.split(sep)
    if len(fields) < len(HEADER_FIELD_NAMES):
        raise ValueError(
            f"Header has {len(fields)} fields, expected {len(HEADER_FIELD_NAMES)}"
        )

    missing_char = "-"
    for f in fields:
        stripped = f.strip()
        if stripped and all(c == "-" for c in stripped):
            missing_char = stripped
            break

    def _int_or_default(val: str, default: int) -> int:
        try:
            return int(val.strip())
        except (ValueError, TypeError):
            return default

    template = LANLHeaderTemplate(
        se_id_start=_int_or_default(fields[0], 150000),
        patient_code=fields[1].strip(),
        pat_id_ssam=fields[2].strip(),
        accession_prefix="".join(c for c in fields[3].strip() if c.isalpha()),
        accession_start=_int_or_default("".join(c for c in fields[3] if c.isdigit()), 0),
        name_prefix="".join(c for c in fields[4].strip() if c.isalpha()),
        name_start=_int_or_default("".join(c for c in fields[4] if c.isdigit()), 1000),
        subtype=fields[5].strip(),
        country=fields[6].strip(),
        sampling_year=fields[7].strip(),
        problematic_sequence=_int_or_default(fields[8], 0),
        hxb2_start=_int_or_default(fields[9], 7050),
        hxb2_stop=_int_or_default(fields[10], 7400),
        organism=fields[12].strip() if len(fields) > 12 else "HIV-1",
        field_separator=sep,
        missing_char=missing_char,
    )

    LOGGER.info("Parsed template from %s: se_id_start=%d, sep=%r", fasta_path, template.se_id_start, sep)
    return template


# section 6 — main generation function


def generate_synthetic_fasta(
    output_path: Path,
    num_sequences: int,
    min_length: int = 250,
    max_length: int = 500,
    alphabet: str = DEFAULT_ALPHABET,
    template: Optional[LANLHeaderTemplate] = None,
    line_width: int = DEFAULT_LINE_WIDTH,
    seed: Optional[int] = None,
    mutation_rate: Optional[float] = None,
    mutation_strategy: Literal["star", "drift"] = "star",
) -> List[SequenceRecord]:
    if num_sequences <= 0:
        raise ValueError(f"num_sequences must be positive, got {num_sequences}")
    if min_length <= 0:
        raise ValueError(f"min_length must be positive, got {min_length}")
    if max_length < min_length:
        raise ValueError(f"max_length ({max_length}) must be >= min_length ({min_length})")
    if not alphabet:
        raise ValueError("alphabet must not be empty")
    if mutation_rate is not None and not (0.0 < mutation_rate <= 1.0):
        raise ValueError(f"mutation_rate must be in (0.0, 1.0], got {mutation_rate}")
    if mutation_strategy not in ("star", "drift"):
        raise ValueError(f"mutation_strategy must be 'star' or 'drift', got {mutation_strategy!r}")

    if template is None:
        template = LANLHeaderTemplate()

    rng = random.Random(seed)
    alphabet_upper = alphabet.upper()

    use_mutation = mutation_rate is not None

    records: List[SequenceRecord] = []

    if use_mutation:
        parent_seq = _generate_sequence(max_length, alphabet_upper, rng)
        prev_seq = parent_seq

        for i in range(num_sequences):
            seq_len = rng.randint(min_length, max_length)

            if i == 0 or mutation_strategy == "star":
                base = _random_subsequence(parent_seq, seq_len, rng)
            else:
                if len(prev_seq) >= seq_len:
                    base = _random_subsequence(prev_seq, seq_len, rng)
                else:
                    base = _generate_sequence(seq_len, alphabet_upper, rng)

            seq = _mutate_sequence(base, mutation_rate, alphabet_upper, rng)
            prev_seq = seq

            header = _build_header(template, i, seq_len)
            records.append(SequenceRecord(header=header, sequence=seq))
    else:
        for i in range(num_sequences):
            seq_len = rng.randint(min_length, max_length)
            seq = _generate_sequence(seq_len, alphabet_upper, rng)
            header = _build_header(template, i, seq_len)
            records.append(SequenceRecord(header=header, sequence=seq))

    write_fasta(records, output_path, line_width=line_width)
    LOGGER.info(
        "Generated %d sequences [%d-%d bp] from alphabet %r -> %s",
        num_sequences, min_length, max_length, alphabet, output_path,
    )

    return records


# section 7 — CLI entry point


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate synthetic FASTA files mimicking LANL HIV database structure.",
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default=Path("synthetic_output.fasta"),
        help="Output FASTA file path.",
    )
    parser.add_argument(
        "--num", "-n",
        type=int,
        default=10,
        help="Number of sequences to generate.",
    )
    parser.add_argument(
        "--min-len",
        type=int,
        default=250,
        help="Minimum sequence length.",
    )
    parser.add_argument(
        "--max-len",
        type=int,
        default=500,
        help="Maximum sequence length.",
    )
    parser.add_argument(
        "--alphabet",
        type=str,
        default=DEFAULT_ALPHABET,
        help="Nucleotide alphabet for sequence generation.",
    )
    parser.add_argument(
        "--template",
        type=Path,
        default=None,
        help="Existing FASTA file to use as header template.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility.",
    )
    parser.add_argument(
        "--line-width",
        type=int,
        default=DEFAULT_LINE_WIDTH,
        help="Line width for sequence wrapping.",
    )
    parser.add_argument(
        "--mutation-rate",
        type=float,
        default=None,
        help="Fraction of positions to mutate per sequence (e.g., 0.02 = 2%%). Enables mutation mode.",
    )
    parser.add_argument(
        "--mutation-strategy",
        choices=["star", "drift"],
        default="star",
        help="Mutation strategy: 'star' = all children from one parent, 'drift' = sequential. Default: star.",
    )
    return parser


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    parser = _build_parser()
    args = parser.parse_args()

    template = None
    if args.template is not None:
        template = parse_template_fasta(args.template)

    records = generate_synthetic_fasta(
        output_path=args.output,
        num_sequences=args.num,
        min_length=args.min_len,
        max_length=args.max_len,
        alphabet=args.alphabet,
        template=template,
        line_width=args.line_width,
        seed=args.seed,
        mutation_rate=args.mutation_rate,
        mutation_strategy=args.mutation_strategy,
    )

    print(f"Generated {len(records)} sequences -> {args.output}")


if __name__ == "__main__":
    main()


# ============================================================================
# Usage
# ============================================================================
#
# CLI:
#   python utils/fasta_file_generator.py -n 10 -o output.fasta
#   python utils/fasta_file_generator.py -n 50 --min-len 200 --max-len 800 --seed 42 -o test.fasta
#   python utils/fasta_file_generator.py -n 20 --template data/synthetic_data/sample.fasta -o from_template.fasta
#   python utils/fasta_file_generator.py -n 5 --alphabet ACGT --min-len 100 --max-len 300 -o custom.fasta
#
# Mutation-based (star phylogeny — default):
#   python utils/fasta_file_generator.py -n 20 --mutation-rate 0.02 --seed 42 -o star.fasta
#
# Mutation-based (sequential drift):
#   python utils/fasta_file_generator.py -n 20 --mutation-rate 0.03 --mutation-strategy drift --seed 42 -o drift.fasta
#
# Programmatic:
#   from pathlib import Path
#   from utils.fasta_file_generator import generate_synthetic_fasta, LANLHeaderTemplate
#
#   records = generate_synthetic_fasta(
#       output_path=Path("synthetic.fasta"),
#       num_sequences=10,
#       min_length=250,
#       max_length=500,
#       alphabet="ACGT",
#       seed=42,
#   )
#
#   template = LANLHeaderTemplate(subtype="B", country="US", sampling_year="2005")
#   records = generate_synthetic_fasta(
#       output_path=Path("subtype_b.fasta"),
#       num_sequences=20,
#       template=template,
#   )
#
#   # Star phylogeny (default): all children mutate independently from one parent
#   records = generate_synthetic_fasta(
#       output_path=Path("star.fasta"),
#       num_sequences=20,
#       mutation_rate=0.02,
#       seed=42,
#   )
#
#   # Sequential drift: each child mutates from the previous one
#   records = generate_synthetic_fasta(
#       output_path=Path("drift.fasta"),
#       num_sequences=20,
#       mutation_rate=0.03,
#       mutation_strategy="drift",
#       seed=42,
#   )

# %% Cell 1 - Imports and configuration
from itertools import combinations
from pathlib import Path

import pandas as pd


# Position-based composition labels from FASTA header
HEADER_FIELDS = [
    "Se ID",  # 1
    "Patient Code",  # 2
    "PAT id(SSAM)",  # 3
    "Accession",  # 4
    "Name",  # 5
    "Subtype",  # 6
    "Country",  # 7
    "Sampling Year",  # 8
    "Problematic Sequence",  # 9
    "HXB2/MAC239 start",  # 10
    "HXB2/MAC239 stop",  # 11
    "Sequence Length",  # 12
    "Organism",  # 13
]

EXPECTED_HEADER_PARTS = len(HEADER_FIELDS)
FIELD_SEPARATOR = ","
MISSING_VALUES = {"", "-", "NA", "N/A", "na", "n/a", "null", "None"}

# Requested input
fasta_path = Path("data/hiv-db-complete-genome-problematic-nucleutides.fasta")
print(f"Configured input file: {fasta_path}")


# %% Cell 2 - FASTA header parser helpers
def normalize_value(value: str):
    if value is None:
        return None
    value = value.strip()
    if value in MISSING_VALUES:
        return None
    return value


def parse_header(header_line: str) -> dict:
    """Parse one FASTA header line into the 13 known composition fields."""
    header = header_line.strip()
    if header.startswith(">"):
        header = header[1:]

    parts = [normalize_value(x) for x in header.split(FIELD_SEPARATOR)]

    row = {}
    for i, field in enumerate(HEADER_FIELDS):
        row[field] = parts[i] if i < len(parts) else None

    # Keep extra columns (if any) so metadata is not silently discarded.
    row["extra_fields_count"] = max(0, len(parts) - EXPECTED_HEADER_PARTS)
    return row


def iter_fasta_headers(path: Path):
    with path.open("r", encoding="utf-8") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                yield line.rstrip("\n")


def load_headers_dataframe(path: Path) -> pd.DataFrame:
    rows = [parse_header(h) for h in iter_fasta_headers(path)]
    if not rows:
        raise ValueError(f"No FASTA headers found in: {path}")
    df = pd.DataFrame(rows)

    # Convert numeric-like columns when possible.
    numeric_cols = [
        "Sampling Year",
        "Problematic Sequence",
        "HXB2/MAC239 start",
        "HXB2/MAC239 stop",
        "Sequence Length",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


# %% Cell 3 - Load and preview parsed label table
if not fasta_path.exists():
    raise FileNotFoundError(f"Input FASTA not found: {fasta_path}")

labels_df = load_headers_dataframe(fasta_path)

print(f"Sequences parsed: {len(labels_df):,}")
print("\nPreview of parsed composition labels:")
print(labels_df.head(10))


# %% Cell 4 - Field quality profile (missingness + uniqueness)
def profile_fields(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    total = len(df)

    for col in HEADER_FIELDS:
        missing_count = int(df[col].isna().sum())
        non_missing = total - missing_count
        unique_non_missing = int(df[col].dropna().nunique())
        fill_rate = (non_missing / total) if total else 0.0
        uniqueness_rate = (
            (unique_non_missing / non_missing)
            if non_missing
            else 0.0
        )

        rows.append(
            {
                "field": col,
                "position": HEADER_FIELDS.index(col) + 1,
                "missing_count": missing_count,
                "missing_pct": (
                    round((missing_count / total) * 100, 2)
                    if total
                    else 0.0
                ),
                "fill_rate": round(fill_rate, 4),
                "unique_non_missing": unique_non_missing,
                "uniqueness_rate": round(uniqueness_rate, 4),
            }
        )

    out = pd.DataFrame(rows)
    out["recommended_for_grouping"] = (
        (out["fill_rate"] >= 0.95)
        & (out["unique_non_missing"] > 1)
    )
    return out.sort_values(
        ["recommended_for_grouping", "fill_rate", "unique_non_missing"],
        ascending=[False, False, False],
    )


field_profile = profile_fields(labels_df)
print(field_profile)


# %% Cell 5 - Rank field combinations for categorization
def candidate_groupings(
    df: pd.DataFrame,
    max_combo_size: int = 3,
) -> pd.DataFrame:
    """Find practical label combinations for categorization.

    Scoring logic favors combinations that:
    1) have low missingness
    2) produce meaningful number of groups
    3) avoid one-record-per-group over-fragmentation
    """
    candidates = []
    total = len(df)

    fields_for_combo = [c for c in HEADER_FIELDS if c != "Se ID"]

    for size in range(1, max_combo_size + 1):
        for combo in combinations(fields_for_combo, size):
            subset = df[list(combo)]
            complete_mask = ~subset.isna().any(axis=1)
            complete_n = int(complete_mask.sum())
            if complete_n == 0:
                continue

            complete_df = subset.loc[complete_mask]
            groups = complete_df.value_counts().rename("n").reset_index()
            n_groups = len(groups)
            mean_group_size = complete_n / n_groups if n_groups else 0
            singleton_pct = (
                (groups["n"].eq(1).mean() * 100)
                if n_groups
                else 100
            )
            coverage = complete_n / total if total else 0

            score = (
                coverage * 0.45
                + min(n_groups / max(10, total * 0.1), 1.0) * 0.25
                + min(mean_group_size / 5, 1.0) * 0.2
                + (1 - singleton_pct / 100) * 0.1
            )

            candidates.append(
                {
                    "combo": " + ".join(combo),
                    "size": size,
                    "coverage": round(coverage, 4),
                    "complete_records": complete_n,
                    "n_groups": n_groups,
                    "mean_group_size": round(mean_group_size, 2),
                    "singleton_pct": round(singleton_pct, 2),
                    "score": round(score, 4),
                }
            )

    out = pd.DataFrame(candidates)
    if out.empty:
        return out
    return out.sort_values(
        ["score", "coverage", "mean_group_size"],
        ascending=[False, False, False],
    )


combo_rank = candidate_groupings(labels_df, max_combo_size=3)
print(combo_rank.head(25))


# %% Cell 6 - Final recommendation table + export CSV outputs
print("Top single-field category candidates:")
top_field_cols = [
    "position",
    "field",
    "fill_rate",
    "unique_non_missing",
    "uniqueness_rate",
    "recommended_for_grouping",
]
print(field_profile[top_field_cols].head(10).to_string(index=False))

print("\nTop multi-field combination candidates:")
if combo_rank.empty:
    print("No valid combinations found.")
else:
    print(combo_rank.head(10).to_string(index=False))

# Non-combo summary: Sampling Year only
single_field = "Sampling Year"
single_series = labels_df[single_field]
single_complete = single_series.dropna()
single_complete_records = int(single_complete.shape[0])
single_n_groups = int(single_complete.nunique())
single_coverage = (
    single_complete_records / len(labels_df)
    if len(labels_df)
    else 0.0
)
single_mean_group_size = (
    single_complete_records / single_n_groups if single_n_groups else 0.0
)

print("\nNon-combo summary (Sampling Year only):")
print(f"coverage={single_coverage:.4f}")
print(f"complete_records={single_complete_records}")
print(f"n_groups={single_n_groups}")
print(f"mean_group_size={single_mean_group_size:.2f}")

sampling_year_counts = (
    single_complete.value_counts()
    .sort_index()
    .rename("n_sequences")
    .to_frame()
)
print("\nGroups per Sampling Year:")
print(sampling_year_counts.to_string())

# Keep outputs out of legacy whole_genome_multialigned
output_dir = Path("data/whole_genome_unaligned/processed")
output_dir.mkdir(parents=True, exist_ok=True)

field_profile.to_csv(output_dir / "composition_field_profile.csv", index=False)
if not combo_rank.empty:
    combo_rank.to_csv(output_dir / "composition_combo_rank.csv", index=False)
sampling_year_counts.to_csv(output_dir / "sampling_year_group_counts.csv")

print(f"\nSaved reports in: {output_dir}")

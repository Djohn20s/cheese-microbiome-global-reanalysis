from pathlib import Path
import pandas as pd
import re

# Repository root (portable; independent of current working directory)
REPO_DIR = Path(__file__).resolve().parents[2]
PROJECT = REPO_DIR
PROJECT_DIR = REPO_DIR
BASE = PROJECT / "data/03_primary/ncbi/cheese_metagenome/v2"

WGS = BASE / "stage5/runinfo_WGS.csv"
AMP = BASE / "stage5/runinfo_AMPLICON.csv"
CENTER_MAP = BASE / "center_country_map_curated.csv"

ORIGIN = (
    BASE /
    "biosample_metadata/audit/geographic_origin/"
    "sra_origin_country_audit.tsv"
)

OUTDIR = (
    BASE /
    "biosample_metadata/audit/geographic_origin/"
    "center_vs_sample_origin"
)
OUTDIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# Helpers
# ============================================================

def clean(x):
    if pd.isna(x):
        return ""
    return re.sub(r"\s+", " ", str(x).strip())


def norm_center(x):
    return clean(x).upper()


def norm_country(x):
    x = clean(x)

    aliases = {
        "United States": "USA",
        "United States of America": "USA",
        "US": "USA",
        "U.S.A.": "USA",
        "UK": "United Kingdom",
        "Great Britain": "United Kingdom",
        "Russian Federation": "Russia",
        "Türkiye": "Turkey",
    }

    return aliases.get(x, x)


# ============================================================
# Load final SRA metadata
# ============================================================

frames = []

for path, strategy in [
    (WGS, "WGS"),
    (AMP, "AMPLICON"),
]:
    x = pd.read_csv(path, dtype=str)

    required = ["Run", "BioSample", "BioProject", "CenterName"]

    missing = [c for c in required if c not in x.columns]

    if missing:
        raise RuntimeError(
            f"{path} missing columns: {missing}"
        )

    x = x[required].copy()
    x["SequencingStrategy"] = strategy

    frames.append(x)

runs = pd.concat(frames, ignore_index=True)

print("Stage5 SRA:", runs["Run"].nunique())
print("Stage5 BioSamples:", runs["BioSample"].nunique())
print("Stage5 BioProjects:", runs["BioProject"].nunique())


# ============================================================
# Load CenterName -> Country mapping
# ============================================================

centers = pd.read_csv(CENTER_MAP, dtype=str)

required = [
    "CenterName",
    "Country",
    "Continent"
]

missing = [c for c in required if c not in centers.columns]

if missing:
    raise RuntimeError(
        f"{CENTER_MAP} missing columns: {missing}"
    )

centers["CenterName_key"] = centers["CenterName"].map(norm_center)

# Check ambiguity in curated mapping
amb = (
    centers.groupby("CenterName_key")["Country"]
    .nunique()
)

amb = amb[amb > 1]

if len(amb):
    print("\nWARNING: ambiguous CenterName -> Country mappings:")
    print(amb.to_string())

centers_small = (
    centers[
        [
            "CenterName_key",
            "CenterName_standardized",
            "Country",
            "Continent"
        ]
    ]
    .drop_duplicates("CenterName_key")
    .rename(
        columns={
            "Country": "CenterCountry",
            "Continent": "CenterContinent"
        }
    )
)


# ============================================================
# Assign sequencing-center country
# ============================================================

runs["CenterName_key"] = runs["CenterName"].map(norm_center)

runs = runs.merge(
    centers_small,
    on="CenterName_key",
    how="left"
)

runs["CenterCountry"] = (
    runs["CenterCountry"]
    .fillna("")
    .map(norm_country)
)

runs["CenterContinent"] = (
    runs["CenterContinent"]
    .fillna("")
)


# ============================================================
# Load BioSample-derived sample origin
# ============================================================

origin = pd.read_csv(ORIGIN, sep="\t", dtype=str)

needed = [
    "Run",
    "SampleOriginCountry",
    "OriginStatus",
    "OriginEvidence"
]

missing = [c for c in needed if c not in origin.columns]

if missing:
    raise RuntimeError(
        f"{ORIGIN} missing columns: {missing}"
    )

origin = (
    origin[needed]
    .drop_duplicates("Run")
)

origin["SampleOriginCountry"] = (
    origin["SampleOriginCountry"]
    .fillna("")
    .map(norm_country)
)

origin["OriginStatus"] = (
    origin["OriginStatus"]
    .fillna("unknown")
)

origin["OriginEvidence"] = (
    origin["OriginEvidence"]
    .fillna("")
)


# ============================================================
# Join
# ============================================================

audit = runs.merge(
    origin,
    on="Run",
    how="left"
)

audit["OriginStatus"] = (
    audit["OriginStatus"]
    .fillna("unknown")
)

audit["SampleOriginCountry"] = (
    audit["SampleOriginCountry"]
    .fillna("")
)

audit["OriginEvidence"] = (
    audit["OriginEvidence"]
    .fillna("")
)


# ============================================================
# Comparison
# ============================================================

def compare(row):

    origin_status = row["OriginStatus"]
    sample_country = clean(row["SampleOriginCountry"])
    center_country = clean(row["CenterCountry"])

    if origin_status != "confirmed" or not sample_country:
        return "no_confirmed_sample_origin"

    if not center_country:
        return "confirmed_origin_no_center_country"

    if sample_country == center_country:
        return "agreement"

    return "disagreement"


audit["Comparison"] = audit.apply(compare, axis=1)


# ============================================================
# Run-level summary
# ============================================================

run_summary = (
    audit["Comparison"]
    .value_counts(dropna=False)
    .rename_axis("Comparison")
    .reset_index(name="n_SRA")
)

run_summary["pct_all_SRA"] = (
    100 * run_summary["n_SRA"] / audit["Run"].nunique()
)


# ============================================================
# Confirmed-origin-only comparison
# ============================================================

confirmed = audit[
    audit["Comparison"].isin(
        ["agreement", "disagreement"]
    )
].copy()

confirmed_summary = (
    confirmed["Comparison"]
    .value_counts()
    .rename_axis("Comparison")
    .reset_index(name="n_SRA")
)

if len(confirmed):
    confirmed_summary["pct_confirmed_comparable_SRA"] = (
        100 *
        confirmed_summary["n_SRA"] /
        confirmed["Run"].nunique()
    )
else:
    confirmed_summary["pct_confirmed_comparable_SRA"] = []


# ============================================================
# Country-pair disagreement table
# ============================================================

disagreement = audit[
    audit["Comparison"] == "disagreement"
].copy()

pair_summary = (
    disagreement
    .groupby(
        ["CenterCountry", "SampleOriginCountry"],
        dropna=False
    )
    .agg(
        n_SRA=("Run", "nunique"),
        n_BioSamples=("BioSample", "nunique"),
        n_BioProjects=("BioProject", "nunique")
    )
    .reset_index()
    .sort_values(
        ["n_SRA", "n_BioProjects"],
        ascending=False
    )
)


# ============================================================
# BioProject-level comparison
# ============================================================

def summarize_project(g):

    confirmed_g = g[
        g["OriginStatus"] == "confirmed"
    ]

    sample_countries = sorted(
        {
            clean(x)
            for x in confirmed_g["SampleOriginCountry"]
            if clean(x)
        }
    )

    center_countries = sorted(
        {
            clean(x)
            for x in g["CenterCountry"]
            if clean(x)
        }
    )

    n_agree = int((g["Comparison"] == "agreement").sum())
    n_disagree = int((g["Comparison"] == "disagreement").sum())

    if len(sample_countries) == 0:
        origin_structure = "no_confirmed_origin"
    elif len(sample_countries) == 1:
        origin_structure = "single_country"
    else:
        origin_structure = "multi_country"

    if n_agree > 0 and n_disagree == 0:
        comparison_status = "all_agreement"
    elif n_disagree > 0 and n_agree == 0:
        comparison_status = "all_disagreement"
    elif n_agree > 0 and n_disagree > 0:
        comparison_status = "mixed"
    else:
        comparison_status = "not_comparable"

    return pd.Series({
        "n_SRA": g["Run"].nunique(),
        "n_BioSamples": g["BioSample"].nunique(),

        "n_confirmed_origin_SRA":
            confirmed_g["Run"].nunique(),

        "SampleOriginStructure":
            origin_structure,

        "SampleOriginCountries":
            ";".join(sample_countries),

        "CenterCountries":
            ";".join(center_countries),

        "n_agreement_SRA":
            n_agree,

        "n_disagreement_SRA":
            n_disagree,

        "ComparisonStatus":
            comparison_status
    })


project_summary = (
    audit
    .groupby("BioProject", group_keys=False)
    .apply(summarize_project)
    .reset_index()
)


# ============================================================
# CenterName audit
# ============================================================

center_unmapped = (
    audit[
        audit["CenterCountry"].eq("")
    ]
    .groupby("CenterName", dropna=False)
    .agg(
        n_SRA=("Run", "nunique"),
        n_BioProjects=("BioProject", "nunique")
    )
    .reset_index()
    .sort_values("n_SRA", ascending=False)
)


# ============================================================
# Save
# ============================================================

audit.to_csv(
    OUTDIR / "sra_center_vs_sample_origin.tsv",
    sep="\t",
    index=False
)

run_summary.to_csv(
    OUTDIR / "comparison_summary_all_sra.tsv",
    sep="\t",
    index=False
)

confirmed_summary.to_csv(
    OUTDIR / "comparison_summary_confirmed.tsv",
    sep="\t",
    index=False
)

pair_summary.to_csv(
    OUTDIR / "country_disagreement_pairs.tsv",
    sep="\t",
    index=False
)

project_summary.to_csv(
    OUTDIR / "bioproject_center_vs_origin.tsv",
    sep="\t",
    index=False
)

center_unmapped.to_csv(
    OUTDIR / "unmapped_centers.tsv",
    sep="\t",
    index=False
)


# ============================================================
# Console
# ============================================================

print("\n==============================================")
print("CENTER COUNTRY vs SAMPLE ORIGIN")
print("==============================================")

print("\nAll SRA:")
print(run_summary.to_string(index=False))

print("\nConfirmed origin + mapped center:")
print(confirmed_summary.to_string(index=False))

print("\nBioProject comparison:")
print(
    project_summary["ComparisonStatus"]
    .value_counts()
    .to_string()
)

print("\nSample-origin structure:")
print(
    project_summary["SampleOriginStructure"]
    .value_counts()
    .to_string()
)

print("\nTop disagreement pairs:")
if len(pair_summary):
    print(pair_summary.head(30).to_string(index=False))
else:
    print("No disagreements.")

print("\nUnmapped CenterName values:")
if len(center_unmapped):
    print(center_unmapped.head(30).to_string(index=False))
else:
    print("None")

print("\nOutputs:", OUTDIR)

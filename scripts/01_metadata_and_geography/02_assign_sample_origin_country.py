from pathlib import Path
import pandas as pd
import re

# Repository root (portable; independent of current working directory)
REPO_DIR = Path(__file__).resolve().parents[2]
PROJECT = REPO_DIR
PROJECT_DIR = REPO_DIR
BASE = PROJECT / "data/03_primary/ncbi/cheese_metagenome/v2"

ATTR = BASE / "biosample_metadata/audit/biosample_attributes_long.tsv"
RUNINFO = BASE / "stage5/runinfo_combined.xlsx"

OUTDIR = BASE / "biosample_metadata/audit/geographic_origin"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Helpers
# ============================================================

MISSING = {
    "", "na", "n/a", "nan", "none", "null",
    "missing", "not applicable", "'not applicable'",
    "not collected", "unknown", "unspecified",
    "restricted access"
}

COUNTRY_ALIASES = {
    "united states": "USA",
    "united states of america": "USA",
    "usa": "USA",
    "u.s.a.": "USA",
    "us": "USA",

    "united kingdom": "United Kingdom",
    "uk": "United Kingdom",
    "great britain": "United Kingdom",

    "france": "France",
    "italy": "Italy",
    "canada": "Canada",
    "norway": "Norway",
    "brazil": "Brazil",
    "switzerland": "Switzerland",
    "greece": "Greece",
    "china": "China",
    "mexico": "Mexico",
    "poland": "Poland",
    "spain": "Spain",
    "australia": "Australia",
    "turkey": "Turkey",
    "türkiye": "Turkey",
    "portugal": "Portugal",
    "russia": "Russia",
    "russian federation": "Russia",
    "iran": "Iran",
    "ecuador": "Ecuador",
    "ireland": "Ireland",
    "argentina": "Argentina",
    "bulgaria": "Bulgaria",
    "cyprus": "Cyprus",
}

def clean(x):
    if pd.isna(x):
        return ""
    return re.sub(r"\s+", " ", str(x).strip())

def norm_attr_name(x):
    x = clean(x).lower()
    x = x.replace("_", " ")
    x = re.sub(r"\s+", " ", x)
    return x

def informative(x):
    x = clean(x)
    return bool(x) and x.lower() not in MISSING

def country_from_geo(value):
    """
    Conservative parser.

    Intended for explicit geographic-location fields such as:
      France
      Italy: Vicenza, Veneto
      USA: Oregon
      Canada: Edmonton

    It does NOT try to infer countries from arbitrary city/region names.
    """
    value = clean(value)

    if not informative(value):
        return None

    # Country is normally before ":" in NCBI geo_loc_name.
    first = re.split(r"[:;,]", value, maxsplit=1)[0].strip()
    key = first.lower().rstrip(".")

    if key in COUNTRY_ALIASES:
        return COUNTRY_ALIASES[key]

    # exact full-value fallback
    key2 = value.lower().rstrip(".")
    if key2 in COUNTRY_ALIASES:
        return COUNTRY_ALIASES[key2]

    return None

# ============================================================
# Load BioSample attributes
# ============================================================

attr = pd.read_csv(ATTR, sep="\t", dtype=str, keep_default_na=False)

attr["attribute_norm"] = attr["attribute_name"].map(norm_attr_name)
attr["value_clean"] = attr["value"].map(clean)

# Only explicit geographic-origin attributes.
geo_names = {
    "geo loc name",
    "geographic location",
    "geographic location (country and/or sea)",
}

geo = attr[
    attr["attribute_norm"].isin(geo_names)
].copy()

geo = geo[geo["value_clean"].map(informative)].copy()

geo["SampleOriginCountry"] = geo["value_clean"].map(country_from_geo)

# ============================================================
# Audit unparsed geographic values
# ============================================================

unparsed = (
    geo[geo["SampleOriginCountry"].isna()]
    .groupby(
        ["attribute_name", "value_clean"],
        dropna=False
    )
    .agg(
        n_BioSamples=("BioSample", "nunique")
    )
    .reset_index()
    .sort_values(
        ["n_BioSamples", "value_clean"],
        ascending=[False, True]
    )
)

unparsed.to_csv(
    OUTDIR / "unparsed_geographic_values.tsv",
    sep="\t",
    index=False
)

# Keep only explicit values that could be converted conservatively.
geo_ok = geo[geo["SampleOriginCountry"].notna()].copy()

# ============================================================
# Resolve multiple geographic attributes within BioSample
# ============================================================

def resolve_biosample(group):
    countries = sorted(set(group["SampleOriginCountry"].dropna()))

    raw_evidence = " | ".join(
        sorted(
            set(
                f"{a}={v}"
                for a, v in zip(
                    group["attribute_name"],
                    group["value_clean"]
                )
            )
        )
    )

    if len(countries) == 1:
        status = "confirmed"
        country = countries[0]
    elif len(countries) > 1:
        status = "conflicting_biosample_metadata"
        country = ""
    else:
        status = "unknown"
        country = ""

    return pd.Series({
        "SampleOriginCountry": country,
        "OriginStatus": status,
        "OriginEvidence": raw_evidence,
        "n_distinct_origin_countries": len(countries),
        "all_origin_countries": ";".join(countries)
    })

if len(geo_ok):
    bio_origin = (
        geo_ok.groupby("BioSample", group_keys=False)
        .apply(resolve_biosample)
        .reset_index()
    )
else:
    bio_origin = pd.DataFrame(columns=[
        "BioSample",
        "SampleOriginCountry",
        "OriginStatus",
        "OriginEvidence",
        "n_distinct_origin_countries",
        "all_origin_countries"
    ])

# ============================================================
# Load final stage5 run metadata
# ============================================================

frames = []

xls = pd.ExcelFile(RUNINFO)

for sheet in ["WGS", "AMPLICON"]:
    d = pd.read_excel(
        RUNINFO,
        sheet_name=sheet,
        dtype=str
    )

    d["SequencingStrategy"] = sheet

    wanted = [
        "Run",
        "BioSample",
        "BioProject",
        "CenterName",
        "Continent",
        "SequencingStrategy"
    ]

    # Country may or may not be present in runinfo_combined.xlsx
    if "Country" in d.columns:
        wanted.insert(4, "Country")

    frames.append(d[wanted].copy())

runs = pd.concat(frames, ignore_index=True)

# ============================================================
# Add old curated/project Country if available elsewhere
# ============================================================

country_candidates = [
    PROJECT / "data/08_reporting/Figure2_BioProject_metadata_curated.csv",
    PROJECT / "data/08_reporting/zenodo_release_final/01_metadata/BioProject_SRA_counts_metadata_summary.csv",
]

project_country = None

for f in country_candidates:
    if not f.exists():
        continue

    try:
        x = pd.read_csv(f, dtype=str)
    except Exception:
        continue

    if "BioProject" in x.columns and "Country" in x.columns:
        project_country = (
            x[["BioProject", "Country"]]
            .dropna()
            .drop_duplicates()
        )
        print("Old Country source:", f)
        break

if "Country" in runs.columns:
    runs = runs.rename(columns={"Country": "Country_original"})
elif project_country is not None:
    runs = runs.merge(
        project_country.rename(
            columns={"Country": "Country_original"}
        ),
        on="BioProject",
        how="left"
    )
else:
    runs["Country_original"] = ""

# ============================================================
# Join confirmed BioSample origin
# ============================================================

run_origin = runs.merge(
    bio_origin,
    on="BioSample",
    how="left"
)

run_origin["OriginStatus"] = (
    run_origin["OriginStatus"]
    .fillna("unknown")
)

run_origin["SampleOriginCountry"] = (
    run_origin["SampleOriginCountry"]
    .fillna("")
)

run_origin["OriginEvidence"] = (
    run_origin["OriginEvidence"]
    .fillna("")
)

run_origin["Country_original"] = (
    run_origin["Country_original"]
    .fillna("")
)

# ============================================================
# Compare with old Country
# ============================================================

def norm_country_for_compare(x):
    x = clean(x)
    if not x:
        return ""

    key = x.lower().rstrip(".")

    if key in COUNTRY_ALIASES:
        return COUNTRY_ALIASES[key]

    return x

run_origin["Country_original_norm"] = (
    run_origin["Country_original"]
    .map(norm_country_for_compare)
)

def comparison_status(row):
    new = clean(row["SampleOriginCountry"])
    old = clean(row["Country_original_norm"])

    if row["OriginStatus"] != "confirmed":
        return row["OriginStatus"]

    if not old:
        return "confirmed_no_old_country"

    if new == old:
        return "confirmed_agreement"

    return "confirmed_disagreement"

run_origin["CountryComparison"] = run_origin.apply(
    comparison_status,
    axis=1
)

# ============================================================
# BioProject-level audit
# ============================================================

def project_summary(group):
    confirmed = group[
        group["OriginStatus"] == "confirmed"
    ]

    countries = sorted(
        set(
            x for x in confirmed["SampleOriginCountry"]
            if clean(x)
        )
    )

    if len(countries) == 0:
        project_origin_status = "no_confirmed_origin"
    elif len(countries) == 1:
        project_origin_status = "single_country"
    else:
        project_origin_status = "multi_country"

    old = sorted(
        set(
            x for x in group["Country_original_norm"]
            if clean(x)
        )
    )

    return pd.Series({
        "n_SRA": group["Run"].nunique(),
        "n_BioSamples": group["BioSample"].nunique(),
        "n_confirmed_origin_SRA": confirmed["Run"].nunique(),
        "n_confirmed_origin_BioSamples":
            confirmed["BioSample"].nunique(),
        "pct_SRA_confirmed_origin":
            100 * confirmed["Run"].nunique()
            / group["Run"].nunique(),
        "n_confirmed_countries": len(countries),
        "confirmed_countries": ";".join(countries),
        "ProjectOriginStatus": project_origin_status,
        "Country_original": ";".join(old),
        "n_agreement_SRA":
            (group["CountryComparison"]
             == "confirmed_agreement").sum(),
        "n_disagreement_SRA":
            (group["CountryComparison"]
             == "confirmed_disagreement").sum(),
    })

project_audit = (
    run_origin
    .groupby("BioProject", group_keys=False)
    .apply(project_summary)
    .reset_index()
)

# ============================================================
# BioProject × confirmed country table
# ============================================================

confirmed = run_origin[
    run_origin["OriginStatus"] == "confirmed"
].copy()

if len(confirmed):
    project_country_groups = (
        confirmed
        .groupby(
            ["BioProject", "SampleOriginCountry"],
            dropna=False
        )
        .agg(
            n_SRA=("Run", "nunique"),
            n_BioSamples=("BioSample", "nunique")
        )
        .reset_index()
        .sort_values(
            ["BioProject", "n_SRA"],
            ascending=[True, False]
        )
    )
else:
    project_country_groups = pd.DataFrame()

# ============================================================
# Summaries
# ============================================================

comparison_summary = (
    run_origin["CountryComparison"]
    .value_counts(dropna=False)
    .rename_axis("CountryComparison")
    .reset_index(name="n_SRA")
)

comparison_summary["pct_SRA"] = (
    100 * comparison_summary["n_SRA"] / len(run_origin)
)

project_status_summary = (
    project_audit["ProjectOriginStatus"]
    .value_counts()
    .rename_axis("ProjectOriginStatus")
    .reset_index(name="n_BioProjects")
)

# ============================================================
# Save
# ============================================================

bio_origin.to_csv(
    OUTDIR / "biosample_origin_country.tsv",
    sep="\t",
    index=False
)

run_origin.to_csv(
    OUTDIR / "sra_origin_country_audit.tsv",
    sep="\t",
    index=False
)

project_audit.to_csv(
    OUTDIR / "bioproject_origin_country_audit.tsv",
    sep="\t",
    index=False
)

project_country_groups.to_csv(
    OUTDIR / "bioproject_country_groups.tsv",
    sep="\t",
    index=False
)

comparison_summary.to_csv(
    OUTDIR / "country_comparison_summary.tsv",
    sep="\t",
    index=False
)

project_status_summary.to_csv(
    OUTDIR / "bioproject_origin_status_summary.tsv",
    sep="\t",
    index=False
)

# ============================================================
# Console report
# ============================================================

n_runs = run_origin["Run"].nunique()
n_bio = run_origin["BioSample"].nunique()
n_proj = run_origin["BioProject"].nunique()

n_confirmed_runs = (
    run_origin.loc[
        run_origin["OriginStatus"] == "confirmed",
        "Run"
    ].nunique()
)

n_confirmed_bio = (
    run_origin.loc[
        run_origin["OriginStatus"] == "confirmed",
        "BioSample"
    ].nunique()
)

print("\n==============================================")
print("GEOGRAPHIC ORIGIN AUDIT")
print("==============================================")
print(f"SRA                    : {n_runs}")
print(f"BioSamples             : {n_bio}")
print(f"BioProjects            : {n_proj}")
print(
    f"Confirmed-origin SRA   : {n_confirmed_runs} "
    f"({100*n_confirmed_runs/n_runs:.2f}%)"
)
print(
    f"Confirmed BioSamples   : {n_confirmed_bio} "
    f"({100*n_confirmed_bio/n_bio:.2f}%)"
)

print("\nSRA comparison with old Country:")
print(comparison_summary.to_string(index=False))

print("\nBioProject origin status:")
print(project_status_summary.to_string(index=False))

multi = project_audit[
    project_audit["ProjectOriginStatus"] == "multi_country"
]

print("\nMulti-country BioProjects:")
if len(multi):
    print(
        multi[
            [
                "BioProject",
                "n_SRA",
                "n_confirmed_origin_SRA",
                "n_confirmed_countries",
                "confirmed_countries",
                "Country_original"
            ]
        ].to_string(index=False)
    )
else:
    print("None")

print("\nUnparsed geographic values:", len(unparsed))
print("\nOutputs:")
for f in sorted(OUTDIR.glob("*.tsv")):
    print(f)

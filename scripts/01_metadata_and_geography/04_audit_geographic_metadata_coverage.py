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
OUTDIR = BASE / "biosample_metadata/audit"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# BioSample attributes
# ------------------------------------------------------------

attr = pd.read_csv(ATTR, sep="\t", dtype=str)

# normalize only the ATTRIBUTE NAME for grouping equivalent spelling/case
def norm_name(x):
    x = str(x).strip().lower()
    x = re.sub(r"[\s_-]+", " ", x)
    return x

attr["attribute_name_normalized"] = attr["attribute_name"].map(norm_name)

# ------------------------------------------------------------
# Final SRA metadata
# ------------------------------------------------------------

frames = []

xls = pd.ExcelFile(RUNINFO)

for sheet in ["WGS", "AMPLICON"]:
    x = pd.read_excel(
        RUNINFO,
        sheet_name=sheet,
        dtype=str
    )

    x["SequencingStrategy"] = sheet

    frames.append(x)

runs = pd.concat(frames, ignore_index=True)

print("Final SRA rows:", len(runs))
print("Unique runs:", runs["Run"].nunique())
print("Unique BioSamples:", runs["BioSample"].nunique())
print("Unique BioProjects:", runs["BioProject"].nunique())

# ------------------------------------------------------------
# Join
# ------------------------------------------------------------

m = attr.merge(
    runs[
        ["Run", "BioSample", "BioProject", "SequencingStrategy"]
    ],
    on="BioSample",
    how="inner"
)

print("Matched BioSamples:", m["BioSample"].nunique())

# ------------------------------------------------------------
# Coverage per normalized attribute
# ------------------------------------------------------------

summary = (
    m.groupby("attribute_name_normalized")
    .agg(
        n_BioSamples=("BioSample", "nunique"),
        n_SRA=("Run", "nunique"),
        n_BioProjects=("BioProject", "nunique"),
        n_attribute_records=("attribute_name", "size"),
    )
    .reset_index()
)

summary["pct_BioSamples"] = (
    100 * summary["n_BioSamples"] / runs["BioSample"].nunique()
)

summary["pct_SRA"] = (
    100 * summary["n_SRA"] / runs["Run"].nunique()
)

summary["pct_BioProjects"] = (
    100 * summary["n_BioProjects"] / runs["BioProject"].nunique()
)

# raw aliases belonging to normalized field
aliases = (
    m.groupby("attribute_name_normalized")["attribute_name"]
    .apply(lambda x: " | ".join(sorted(set(x.dropna()))))
    .rename("raw_attribute_names")
    .reset_index()
)

summary = summary.merge(
    aliases,
    on="attribute_name_normalized",
    how="left"
)

summary = summary.sort_values(
    ["n_BioProjects", "n_BioSamples"],
    ascending=[False, False]
)

summary.to_csv(
    OUTDIR / "attribute_coverage_by_project.tsv",
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# Exact BioProjects contributing each attribute
# ------------------------------------------------------------

projects = (
    m.groupby("attribute_name_normalized")["BioProject"]
    .apply(lambda x: ";".join(sorted(set(x.dropna()))))
    .rename("BioProjects")
    .reset_index()
)

projects.to_csv(
    OUTDIR / "attribute_bioproject_membership.tsv",
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# Coverage split WGS / AMPLICON
# ------------------------------------------------------------

strategy = (
    m.groupby(
        ["attribute_name_normalized", "SequencingStrategy"]
    )
    .agg(
        n_SRA=("Run", "nunique"),
        n_BioSamples=("BioSample", "nunique"),
        n_BioProjects=("BioProject", "nunique")
    )
    .reset_index()
)

strategy.to_csv(
    OUTDIR / "attribute_coverage_by_strategy.tsv",
    sep="\t",
    index=False
)

print()
print("Written:")
print(OUTDIR / "attribute_coverage_by_project.tsv")
print(OUTDIR / "attribute_bioproject_membership.tsv")
print(OUTDIR / "attribute_coverage_by_strategy.tsv")

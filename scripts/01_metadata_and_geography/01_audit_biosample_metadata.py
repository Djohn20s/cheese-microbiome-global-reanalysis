from pathlib import Path
import xml.etree.ElementTree as ET
import pandas as pd
from collections import Counter

# Repository root (portable; independent of current working directory)
REPO_DIR = Path(__file__).resolve().parents[2]
PROJECT = REPO_DIR
PROJECT_DIR = REPO_DIR
BASE = PROJECT / "data/03_primary/ncbi/cheese_metagenome/v2/biosample_metadata"
XMLDIR = BASE / "xml"
OUTDIR = BASE / "audit"
OUTDIR.mkdir(parents=True, exist_ok=True)

rows = []
errors = []

files = sorted(XMLDIR.glob("*.xml"))
print(f"XML files: {len(files)}")

for i, path in enumerate(files, 1):
    try:
        root = ET.parse(path).getroot()
        bs = root.find(".//BioSample")

        if bs is None:
            errors.append((path.name, "BioSample element not found"))
            continue

        accession = bs.attrib.get("accession", "")
        biosample_id = bs.attrib.get("id", "")

        title_el = bs.find("./Description/Title")
        title = title_el.text.strip() if title_el is not None and title_el.text else ""

        organism_el = bs.find("./Description/Organism")
        organism = ""
        taxid = ""

        if organism_el is not None:
            organism = organism_el.attrib.get("taxonomy_name", "")
            taxid = organism_el.attrib.get("taxonomy_id", "")

        package_el = bs.find("./Package")
        package = package_el.text.strip() if package_el is not None and package_el.text else ""

        for attr in bs.findall("./Attributes/Attribute"):
            value = (attr.text or "").strip()

            rows.append({
                "BioSample": accession,
                "BioSample_ID": biosample_id,
                "Title": title,
                "Organism": organism,
                "TaxID": taxid,
                "Package": package,
                "attribute_name": attr.attrib.get("attribute_name", ""),
                "harmonized_name": attr.attrib.get("harmonized_name", ""),
                "display_name": attr.attrib.get("display_name", ""),
                "value": value
            })

    except Exception as e:
        errors.append((path.name, str(e)))

    if i % 500 == 0:
        print(f"Processed {i}/{len(files)}")

df = pd.DataFrame(rows)

# pełna tabela long
df.to_csv(
    OUTDIR / "biosample_attributes_long.tsv",
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# AUDIT NAZW ATRYBUTÓW
# ------------------------------------------------------------

audit_rows = []

for name, g in df.groupby("attribute_name", dropna=False):

    values = (
        g["value"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    values = values[values != ""]

    counts = Counter(values)

    examples = [
        f"{v} (n={n})"
        for v, n in counts.most_common(10)
    ]

    harmonized = sorted(
        set(
            x for x in g["harmonized_name"].dropna().astype(str)
            if x.strip()
        )
    )

    display = sorted(
        set(
            x for x in g["display_name"].dropna().astype(str)
            if x.strip()
        )
    )

    audit_rows.append({
        "attribute_name": name,
        "n_BioSamples": g["BioSample"].nunique(),
        "percent_BioSamples": 100 * g["BioSample"].nunique() / len(files),
        "n_records": len(g),
        "n_unique_values": values.nunique(),
        "harmonized_names": " | ".join(harmonized),
        "display_names": " | ".join(display),
        "top_values": " | ".join(examples)
    })

audit = pd.DataFrame(audit_rows)

audit = audit.sort_values(
    ["n_BioSamples", "attribute_name"],
    ascending=[False, True]
)

audit.to_csv(
    OUTDIR / "attribute_name_audit.tsv",
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# PODSTAWOWA TABELA BIOSAMPLE
# ------------------------------------------------------------

basic = (
    df[
        ["BioSample", "BioSample_ID", "Title",
         "Organism", "TaxID", "Package"]
    ]
    .drop_duplicates()
)

basic.to_csv(
    OUTDIR / "biosample_basic_metadata.tsv",
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# ERRORS
# ------------------------------------------------------------

pd.DataFrame(
    errors,
    columns=["file", "error"]
).to_csv(
    OUTDIR / "parse_errors.tsv",
    sep="\t",
    index=False
)

print()
print("==============================================")
print("AUDIT COMPLETE")
print("==============================================")
print(f"XML files            : {len(files)}")
print(f"Parsed BioSamples    : {df['BioSample'].nunique()}")
print(f"Attribute records    : {len(df)}")
print(f"Unique attribute names: {df['attribute_name'].nunique()}")
print(f"Parse errors         : {len(errors)}")
print()
print("Outputs:")
print(OUTDIR / "biosample_attributes_long.tsv")
print(OUTDIR / "attribute_name_audit.tsv")
print(OUTDIR / "biosample_basic_metadata.tsv")
print(OUTDIR / "parse_errors.tsv")

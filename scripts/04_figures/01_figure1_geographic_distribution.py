#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Figure 1 maps for final curated cheese microbiome dataset, Stage 5 only.

Outputs:
- world map with continent-level pies: total BioProjects by continent, pie = AMPLICON vs WGS
- Europe + Turkey map without Russia, country-level pies: total BioProjects by country, pie = AMPLICON vs WGS
- combined world + Europe panel

Important metadata logic:
- runinfo_WGS.csv and runinfo_AMPLICON.csv do not need to contain Country.
- Country and curated Continent are added from center_country_map_curated.csv by normalized CenterName.
- If curated Continent is missing, fallback to runinfo Continent.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, List, Optional

import re
import unicodedata
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle

try:
    import geodatasets
except ImportError as exc:
    raise SystemExit(
        "Missing package: geodatasets. Install with: mamba install -c conda-forge geodatasets geopandas"
    ) from exc

try:
    import cartopy.io.shapereader as shpreader
except ImportError:
    shpreader = None


# ============================================================
# Paths
# ============================================================

# Repository root (portable; independent of current working directory)
REPO_DIR = Path(__file__).resolve().parents[2]
PROJECT = REPO_DIR
PROJECT_DIR = REPO_DIR
WGS_FILE = PROJECT_DIR / "metadata/curation/runinfo_WGS.csv"
AMP_FILE = PROJECT_DIR / "metadata/curation/runinfo_AMPLICON.csv"
ORIGIN_AUDIT = PROJECT_DIR / "metadata/geography/sra_origin_country_audit.tsv"

OUT_DIR = PROJECT_DIR / "results/figures/figure1"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# Configuration
# ============================================================

CONTINENT_CENTROIDS: Dict[str, Tuple[float, float]] = {
    "North America": (-100, 45),
    "South America": (-60, -18),
    "Europe": (15, 52),
    "Africa": (20, 5),
    "Asia": (90, 40),
    "Australia/Oceania": (135, -25),
    "Oceania": (135, -25),
}

# Approximate country label/pie positions in Europe panel.
# Only used for countries present in final data.
EUROPE_COUNTRY_CENTROIDS: Dict[str, Tuple[float, float]] = {
    "France": (2.5, 46.5),
    "Italy": (12.5, 42.8),
    "Spain": (-3.5, 40.2),
    "Portugal": (-8.0, 39.5),
    "Belgium": (4.7, 50.8),
    "Switzerland": (8.2, 46.8),
    "Ireland": (-8.0, 53.3),
    "United Kingdom": (-2.5, 54.5),
    "Norway": (10.0, 62.0),
    "Denmark": (10.0, 56.0),
    "Poland": (19.0, 52.0),
    "Greece": (22.0, 39.0),
    "Turkey": (35.0, 39.0),
    "Austria": (14.0, 47.6),
    "Slovakia": (19.5, 48.8),
    "Cyprus": (33.0, 35.0),
    "Bulgaria": (25.0, 42.8),
    "Unknown": (30.0, 56.0),
}

COUNTRY_NAME_MAP = {
    "United States of America": "USA",
    "United States": "USA",
    "Russian Federation": "Russia",
    "Türkiye": "Turkey",
    "Turkiye": "Turkey",
    "Republic of Türkiye": "Turkey",
    "United Kingdom": "United Kingdom",
}

WORLD_EXCLUDE = {"Antarctica"}
EUROPE_EXCLUDE_COUNTRIES = {"Russia", "Russian Federation"}
EUROPE_INCLUDE_COUNTRIES = {
    "France", "Italy", "Spain", "Portugal", "Belgium", "Switzerland", "Ireland",
    "United Kingdom", "Norway", "Denmark", "Poland", "Greece", "Turkey",
    "Austria", "Slovakia", "Cyprus", "Bulgaria", "Netherlands", "Germany",
    "Czechia", "Croatia", "Serbia", "Romania", "Hungary", "Slovenia",
    "Bosnia and Herzegovina", "Montenegro", "Albania", "North Macedonia",
    "Kosovo", "Moldova", "Ukraine", "Belarus", "Lithuania", "Latvia",
    "Estonia", "Finland", "Sweden", "Iceland", "Luxembourg", "Liechtenstein",
    "Monaco", "Andorra", "Malta", "San Marino", "Vatican", "Czech Republic",
}


@dataclass(frozen=True)
class PieMapConfig:
    strategies: Tuple[str, str] = ("AMPLICON", "WGS")
    colors: Tuple[str, str] = ("#4E79A7", "#F28E2B")
    continent_radius_scale: float = 0.95
    continent_radius_min: float = 2.0
    continent_radius_max: float = 10.5
    country_radius_scale: float = 0.45
    country_radius_min: float = 0.45
    country_radius_max: float = 2.8
    min_total_to_plot: int = 1
    output_dpi: int = 600


# ============================================================
# Helpers
# ============================================================

def normalize_center(x: object) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    s = re.sub(r"[\"'\.,\(\)]", "", s)
    s = " ".join(s.split())
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii")
    return s.upper().strip()


def standardize_unknown(s: pd.Series) -> pd.Series:
    out = s.astype("string").str.strip()
    out = out.mask(out.isna() | (out == ""), "Unknown")
    return out.astype(str)


def find_country_name_column(world: gpd.GeoDataFrame) -> str:
    candidates = [
        "NAME", "name", "ADMIN", "admin", "NAME_EN", "name_en", "SOVEREIGNT",
        "sovereignt", "NAME_LONG", "name_long", "FORMAL_EN", "formal_en",
    ]
    for col in candidates:
        if col in world.columns:
            return col

    object_cols = [c for c in world.columns if world[c].dtype == "object" and c != "geometry"]
    if not object_cols:
        raise ValueError(f"Could not identify country-name column. Available columns: {list(world.columns)}")

    # Prefer the object column with most unique values; avoids using feature class when possible.
    object_cols_sorted = sorted(object_cols, key=lambda c: world[c].nunique(dropna=True), reverse=True)
    return object_cols_sorted[0]


def load_world(require_country_names: bool = True) -> gpd.GeoDataFrame:
    """Load Natural Earth country polygons.

    geodatasets.get_path("naturalearth.land") is only a land mask in recent
    versions and has no country-name column. Therefore we first try country
    datasets and then Cartopy Natural Earth admin_0 countries.
    """

    candidate_keys = [
        "naturalearth.countries",
        "naturalearth.admin_0_countries",
        "naturalearth_lowres",
    ]

    last_errors = []

    for key in candidate_keys:
        try:
            path = geodatasets.get_path(key)
            world = gpd.read_file(path)
            name_col = find_country_name_column(world)
            print(f"Natural Earth country polygons loaded from geodatasets: {key}")
            print(f"Country-name column used: {name_col}")
            world = world.copy()
            world["plot_name"] = world[name_col].astype(str).replace(COUNTRY_NAME_MAP)
            world = world[~world["plot_name"].isin(WORLD_EXCLUDE)].copy()
            return world
        except Exception as exc:
            last_errors.append(f"{key}: {exc}")

    if shpreader is not None:
        try:
            shp = shpreader.natural_earth(
                resolution="110m",
                category="cultural",
                name="admin_0_countries",
            )
            world = gpd.read_file(shp)
            name_col = find_country_name_column(world)
            print("Natural Earth country polygons loaded from Cartopy admin_0_countries")
            print(f"Country-name column used: {name_col}")
            world = world.copy()
            world["plot_name"] = world[name_col].astype(str).replace(COUNTRY_NAME_MAP)
            world = world[~world["plot_name"].isin(WORLD_EXCLUDE)].copy()
            return world
        except Exception as exc:
            last_errors.append(f"cartopy admin_0_countries: {exc}")

    # Last-resort fallback: land mask. This is OK for the world background,
    # but not sufficient for Europe without Russia.
    try:
        path = geodatasets.get_path("naturalearth.land")
        world = gpd.read_file(path).copy()
        world["plot_name"] = "Land"
        if require_country_names:
            raise ValueError(
                "Only naturalearth.land is available, which has no country names. "
                "Install country polygons with: mamba install -c conda-forge cartopy "
                "or use a geodatasets version containing naturalearth.countries. "
                f"Previous attempts: {' | '.join(last_errors)}"
            )
        print("WARNING: using naturalearth.land fallback without country names.")
        return world
    except Exception as exc:
        last_errors.append(f"naturalearth.land fallback: {exc}")

    raise RuntimeError(
        "Could not load Natural Earth country polygons. Install with:\n"
        "mamba install -c conda-forge cartopy geopandas geodatasets pyogrio shapely fiona\n"
        f"Details: {' | '.join(last_errors)}"
    )


def country_to_continent(country: str) -> str:
    """Map confirmed sample-origin country to continent."""

    mapping = {
        # Europe
        "Austria": "Europe",
        "Belgium": "Europe",
        "Bulgaria": "Europe",
        "Croatia": "Europe",
        "Cyprus": "Europe",
        "Czechia": "Europe",
        "Czech Republic": "Europe",
        "Denmark": "Europe",
        "Estonia": "Europe",
        "Finland": "Europe",
        "France": "Europe",
        "Germany": "Europe",
        "Greece": "Europe",
        "Hungary": "Europe",
        "Iceland": "Europe",
        "Ireland": "Europe",
        "Italy": "Europe",
        "Latvia": "Europe",
        "Lithuania": "Europe",
        "Luxembourg": "Europe",
        "Malta": "Europe",
        "Moldova": "Europe",
        "Montenegro": "Europe",
        "Netherlands": "Europe",
        "North Macedonia": "Europe",
        "Norway": "Europe",
        "Poland": "Europe",
        "Portugal": "Europe",
        "Romania": "Europe",
        "Russia": "Europe",
        "Russian Federation": "Europe",
        "Serbia": "Europe",
        "Slovakia": "Europe",
        "Slovenia": "Europe",
        "Spain": "Europe",
        "Sweden": "Europe",
        "Switzerland": "Europe",
        "Turkey": "Europe",
        "Türkiye": "Europe",
        "Ukraine": "Europe",
        "United Kingdom": "Europe",

        # North America
        "USA": "North America",
        "United States": "North America",
        "United States of America": "North America",
        "Canada": "North America",
        "Mexico": "North America",

        # South America
        "Argentina": "South America",
        "Brazil": "South America",
        "Ecuador": "South America",

        # Asia
        "China": "Asia",
        "Iran": "Asia",
        "Kazakhstan": "Asia",

        # Australia/Oceania
        "Australia": "Australia/Oceania",
        "New Zealand": "Australia/Oceania",

        # Africa
        "South Africa": "Africa",
    }

    return mapping.get(str(country).strip(), "Unknown")


def load_final_runinfo() -> pd.DataFrame:
    """
    Load the final Stage-5 SRA records and replace sequencing-center
    geography with confirmed BioSample sample-origin geography.

    Only SRA records with OriginStatus == 'confirmed' and a non-empty
    SampleOriginCountry are retained for geographic plotting.

    BioProjects containing confirmed samples from more than one country
    are preserved at the BioProject x SampleOriginCountry level rather
    than being collapsed to the sequencing-center country.
    """

    if not WGS_FILE.exists():
        raise FileNotFoundError(f"Missing WGS file: {WGS_FILE}")
    if not AMP_FILE.exists():
        raise FileNotFoundError(f"Missing AMPLICON file: {AMP_FILE}")
    if not ORIGIN_AUDIT.exists():
        raise FileNotFoundError(f"Missing origin audit: {ORIGIN_AUDIT}")

    wgs = pd.read_csv(WGS_FILE, low_memory=False)
    amp = pd.read_csv(AMP_FILE, low_memory=False)

    wgs["strategy"] = "WGS"
    amp["strategy"] = "AMPLICON"

    runinfo = pd.concat([wgs, amp], ignore_index=True, sort=False)

    required = {"Run", "BioProject"}
    missing = required - set(runinfo.columns)
    if missing:
        raise ValueError(
            f"Stage-5 runinfo missing required columns: {sorted(missing)}"
        )

    runinfo["Run"] = runinfo["Run"].astype(str).str.strip()
    runinfo["BioProject"] = runinfo["BioProject"].astype(str).str.strip()

    origin = pd.read_csv(
        ORIGIN_AUDIT,
        sep="\t",
        low_memory=False
    )

    required_origin = {
        "Run",
        "BioProject",
        "SampleOriginCountry",
        "OriginStatus",
    }

    missing = required_origin - set(origin.columns)
    if missing:
        raise ValueError(
            f"Origin audit missing required columns: {sorted(missing)}"
        )

    origin["Run"] = origin["Run"].astype(str).str.strip()
    origin["BioProject"] = origin["BioProject"].astype(str).str.strip()

    # One geographic assignment per SRA record.
    origin_small = (
        origin[
            [
                "Run",
                "BioProject",
                "SampleOriginCountry",
                "OriginStatus",
            ]
        ]
        .drop_duplicates()
    )

    # Merge by Run; BioProject from Stage5 remains authoritative.
    runinfo = runinfo.merge(
        origin_small.drop(columns=["BioProject"]),
        on="Run",
        how="left",
        validate="one_to_one",
    )

    # Confirmed sample-origin geography only.
    confirmed = (
        runinfo["OriginStatus"].astype(str).str.lower().eq("confirmed")
        & runinfo["SampleOriginCountry"].notna()
        & runinfo["SampleOriginCountry"].astype(str).str.strip().ne("")
    )

    runinfo = runinfo.loc[confirmed].copy()


    # ========================================================
    # Keep only successfully profiled BioProjects
    # ========================================================

    profile_dir = (
        PROJECT_DIR
        / "processed_data/bioproject_matrices"
    )

    # Amplicon: successfully profiled by QIIME2
    qiime2_bp = pd.read_csv(
        profile_dir / "QIIME2_Amplicon_genus_BioProject_taxa_wide.csv",
        usecols=["BioProject"],
    )

    profiled_amp = set(
        qiime2_bp["BioProject"]
        .dropna()
        .astype(str)
        .unique()
    )

    # WGS: successfully profiled by Kraken2 OR mOTUs
    kraken_bp = pd.read_csv(
        profile_dir / "Kraken2_WGS_genus_BioProject_taxa_wide.csv",
        usecols=["BioProject"],
    )

    motus_bp = pd.read_csv(
        profile_dir / "mOTUs_WGS_genus_BioProject_taxa_wide.csv",
        usecols=["BioProject"],
    )

    profiled_wgs = (
        set(kraken_bp["BioProject"].dropna().astype(str).unique())
        |
        set(motus_bp["BioProject"].dropna().astype(str).unique())
    )

    print("\nSuccessfully profiled BioProjects:")
    print(f"  Amplicon / QIIME2: {len(profiled_amp)}")
    print(f"  WGS / Kraken2 or mOTUs: {len(profiled_wgs)}")

    before_amp = runinfo.loc[
        runinfo["strategy"].eq("AMPLICON"), "BioProject"
    ].nunique()

    before_wgs = runinfo.loc[
        runinfo["strategy"].eq("WGS"), "BioProject"
    ].nunique()

    keep_profiled = (
        (
            runinfo["strategy"].eq("AMPLICON")
            & runinfo["BioProject"].astype(str).isin(profiled_amp)
        )
        |
        (
            runinfo["strategy"].eq("WGS")
            & runinfo["BioProject"].astype(str).isin(profiled_wgs)
        )
    )

    runinfo = runinfo.loc[keep_profiled].copy()

    after_amp = runinfo.loc[
        runinfo["strategy"].eq("AMPLICON"), "BioProject"
    ].nunique()

    after_wgs = runinfo.loc[
        runinfo["strategy"].eq("WGS"), "BioProject"
    ].nunique()

    print("\nConfirmed-origin BioProjects before -> after profiling filter:")
    print(f"  Amplicon: {before_amp} -> {after_amp}")
    print(f"  WGS:      {before_wgs} -> {after_wgs}")

    runinfo["Country"] = (
        runinfo["SampleOriginCountry"]
        .astype(str)
        .str.strip()
        .replace(COUNTRY_NAME_MAP)
    )

    runinfo["Continent"] = runinfo["Country"].map(country_to_continent)

    runinfo["strategy"] = (
        runinfo["strategy"]
        .astype(str)
        .str.upper()
    )

    # --------------------------------------------------------
    # Diagnostics
    # --------------------------------------------------------

    print("")
    print("Confirmed sample-origin dataset loaded")
    print("--------------------------------------")
    print(f"SRA records with confirmed origin: {runinfo['Run'].nunique():,}")
    print(f"BioProjects represented: {runinfo['BioProject'].nunique():,}")
    print(f"Countries represented: {runinfo['Country'].nunique():,}")

    print("")
    print("SRA by strategy:")
    print(
        runinfo.groupby("strategy")["Run"]
        .nunique()
        .to_string()
    )

    print("")
    print("BioProjects by strategy:")
    print(
        runinfo.groupby("strategy")["BioProject"]
        .nunique()
        .to_string()
    )

    # BioProject x country table: critical for multi-country projects.
    project_country = (
        runinfo[
            [
                "BioProject",
                "Country",
                "Continent",
                "strategy",
            ]
        ]
        .drop_duplicates()
        .sort_values(
            ["BioProject", "Country", "strategy"]
        )
    )

    project_country.to_csv(
        OUT_DIR / "Source_Data_Figure1_BioProject_sample_origin_country_strategy.csv",
        index=False
    )

    country_counts = (
        project_country
        .groupby(["Country", "Continent"], dropna=False)
        .agg(
            n_BioProjects=("BioProject", "nunique")
        )
        .reset_index()
        .sort_values(
            ["Continent", "n_BioProjects", "Country"],
            ascending=[True, False, True]
        )
    )

    country_counts.to_csv(
        OUT_DIR / "Source_Data_Figure1_confirmed_origin_country_counts.csv",
        index=False
    )

    # Multi-country projects
    multi = (
        project_country[
            ["BioProject", "Country"]
        ]
        .drop_duplicates()
        .groupby("BioProject")
        .agg(
            n_countries=("Country", "nunique"),
            countries=(
                "Country",
                lambda x: ";".join(sorted(set(x)))
            )
        )
        .reset_index()
    )

    multi = multi[multi["n_countries"] > 1].copy()

    multi.to_csv(
        OUT_DIR / "Source_Data_Figure1_multicountry_BioProjects.csv",
        index=False
    )

    print("")
    print("Multi-country BioProjects:")
    if multi.empty:
        print("None")
    else:
        print(multi.to_string(index=False))

    unknown = runinfo[runinfo["Continent"] == "Unknown"].copy()

    if not unknown.empty:
        unknown[
            ["Run", "BioProject", "Country", "strategy"]
        ].to_csv(
            OUT_DIR / "Figure1_UNKNOWN_continent_check.csv",
            index=False
        )

        print("")
        print("WARNING: confirmed countries without continent mapping:")
        print(sorted(unknown["Country"].unique()))
    else:
        print("")
        print("All confirmed sample-origin countries mapped to continents.")

    return runinfo


def summarize_bioprojects(runinfo: pd.DataFrame, group_col: str) -> pd.DataFrame:
    """
    Summarize BioProjects by geographic group and sequencing strategy.

    n_bioprojects_strategy:
        number of distinct BioProjects represented by a given strategy.

    n_unique_bioprojects:
        number of distinct BioProjects in the geographic group irrespective
        of sequencing strategy. This value is used for the displayed n and
        pie radius, preventing projects represented by both AMPLICON and WGS
        from being counted twice.
    """

    strategy_counts = (
        runinfo.groupby([group_col, "strategy"], dropna=False)
        .agg(
            n_sra=("Run", "nunique"),
            n_bioprojects_strategy=("BioProject", "nunique")
        )
        .reset_index()
    )

    unique_counts = (
        runinfo.groupby(group_col, dropna=False)
        .agg(
            n_unique_bioprojects=("BioProject", "nunique")
        )
        .reset_index()
    )

    df = strategy_counts.merge(
        unique_counts,
        on=group_col,
        how="left"
    )

    df = df.rename(columns={group_col: "group"})
    df["group"] = df["group"].fillna("Unknown").astype(str)

    return df

def pivot_for_pies(summary: pd.DataFrame, cfg: PieMapConfig) -> pd.DataFrame:
    """
    Prepare strategy slices for pie charts while keeping the displayed
    total as the number of UNIQUE BioProjects.

    Strategy counts are intentionally retained as overlapping indicators:
    a BioProject represented by both AMPLICON and WGS contributes to both
    strategy slices, but only once to n_unique_bioprojects.
    """

    piv = (
        summary.pivot_table(
            index="group",
            columns="strategy",
            values="n_bioprojects_strategy",
            aggfunc="sum"
        )
        .fillna(0)
    )

    for s in cfg.strategies:
        if s not in piv.columns:
            piv[s] = 0

    piv = piv[list(cfg.strategies)]

    unique_n = (
        summary.groupby("group")["n_unique_bioprojects"]
        .max()
    )

    piv["total"] = unique_n.reindex(piv.index)

    return piv

def radius_from_total(total: float, scale: float, rmin: float, rmax: float) -> float:
    r = scale * np.sqrt(max(total, 0.0))
    return float(np.clip(r, rmin, rmax))


def draw_pie(ax, x: float, y: float, values: np.ndarray, colors: List[str], r: float):
    total = float(np.nansum(values))
    if total <= 0:
        return
    angles = np.cumsum(values) / total * 360.0
    start = 0.0
    for v, end, c in zip(values, angles, colors):
        if v <= 0:
            continue
        ax.add_patch(Wedge((x, y), r, start, float(end), facecolor=c, edgecolor="white", linewidth=0.7, zorder=5))
        start = float(end)
    ax.add_patch(Circle((x, y), r, facecolor="none", edgecolor="white", linewidth=0.9, zorder=6))


def add_strategy_legend(ax, cfg: PieMapConfig, loc="upper left", bbox_to_anchor=(1.02, 1.0)):
    handles = [Circle((0, 0), radius=0.2, facecolor=c, edgecolor="none") for c in cfg.colors]
    ax.legend(handles, list(cfg.strategies), loc=loc, bbox_to_anchor=bbox_to_anchor, frameon=False, title="Strategy")


def plot_world_continent_pies(by_continent: pd.DataFrame, cfg: PieMapConfig, ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 6))

    world = load_world(require_country_names=False)
    world.plot(ax=ax, linewidth=0.45, edgecolor="0.80", color="0.95", zorder=1)

    piv = pivot_for_pies(by_continent, cfg)

    for cont, row in piv.iterrows():
        if cont not in CONTINENT_CENTROIDS:
            continue
        total = float(row["total"])
        if total < cfg.min_total_to_plot:
            continue
        x, y = CONTINENT_CENTROIDS[cont]
        r = radius_from_total(total, cfg.continent_radius_scale, cfg.continent_radius_min, cfg.continent_radius_max)
        vals = row[list(cfg.strategies)].to_numpy(dtype=float)
        draw_pie(ax, x, y, vals, list(cfg.colors), r)
        ax.text(x, y - r - 2.0, f"{cont}\n(n={int(round(total))})", ha="center", va="top", fontsize=8.5, zorder=10)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 85)
    ax.set_axis_off()
    ax.set_title("A", loc="left", fontweight="bold", fontsize=14)
    return ax


def plot_europe_country_pies(by_country: pd.DataFrame, cfg: PieMapConfig, ax=None):
    """
    Panel B: Europe at country-level resolution.

    Geography is based exclusively on confirmed SampleOriginCountry.
    Russia is excluded from the displayed map.
    Turkey is retained.
    """

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 7.5))

    world = load_world(require_country_names=True)
    world["plot_name"] = world["plot_name"].replace(COUNTRY_NAME_MAP)

    # --------------------------------------------------------
    # European map background
    # Russia deliberately excluded.
    # Turkey deliberately included.
    # --------------------------------------------------------

    europe = world[
        world["plot_name"].isin(EUROPE_INCLUDE_COUNTRIES)
        & ~world["plot_name"].isin(
            {"Russia", "Russian Federation"}
        )
    ].copy()

    europe.plot(
        ax=ax,
        linewidth=0.55,
        edgecolor="0.70",
        color="0.95",
        zorder=1
    )

    # --------------------------------------------------------
    # BioProject counts by confirmed sample-origin country
    # --------------------------------------------------------

    piv = pivot_for_pies(by_country, cfg)

    countries_plotted = []

    for country, row in piv.iterrows():

        # Russia is not shown in Panel B
        if country in {"Russia", "Russian Federation"}:
            continue

        # Only countries for which a European plotting position exists
        if country not in EUROPE_COUNTRY_CENTROIDS:
            continue

        total = float(row["total"])

        if total < cfg.min_total_to_plot:
            continue

        x, y = EUROPE_COUNTRY_CENTROIDS[country]

        r = radius_from_total(
            total,
            cfg.country_radius_scale,
            cfg.country_radius_min,
            cfg.country_radius_max
        )

        vals = row[list(cfg.strategies)].to_numpy(dtype=float)

        draw_pie(
            ax,
            x,
            y,
            vals,
            list(cfg.colors),
            r
        )

        ax.text(
            x,
            y - r - 0.55,
            f"{country}\n(n={int(round(total))})",
            ha="center",
            va="top",
            fontsize=6.6,
            zorder=10
        )

        countries_plotted.append(country)

    # --------------------------------------------------------
    # Geographic extent
    #
    # Western Europe -> Turkey
    # Scandinavia -> Mediterranean
    # Russia excluded
    # --------------------------------------------------------

    ax.set_xlim(-13, 43)
    ax.set_ylim(33, 72)

    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()

    ax.set_title(
        "B",
        loc="left",
        fontweight="bold",
        fontsize=14
    )

    print("")
    print("Panel B — European countries plotted:")
    print(", ".join(sorted(countries_plotted)))

    missing_europe = sorted(
        set(piv.index)
        .intersection(EUROPE_INCLUDE_COUNTRIES)
        .difference(countries_plotted)
        .difference({"Russia", "Russian Federation"})
    )

    if missing_europe:
        print("")
        print("WARNING: European countries present in data but not plotted:")
        print(", ".join(missing_europe))

    return ax


def main():
    cfg = PieMapConfig()
    runinfo = load_final_runinfo()

    by_continent = summarize_bioprojects(runinfo, "Continent")
    by_country = summarize_bioprojects(runinfo, "Country")

    by_continent.to_csv(OUT_DIR / "Source_Data_Figure1A_continent_strategy_BioProjects.csv", index=False)
    by_country.to_csv(OUT_DIR / "Source_Data_Figure1B_country_strategy_BioProjects.csv", index=False)

    # World only
    fig, ax = plt.subplots(figsize=(12, 6))
    plot_world_continent_pies(by_continent, cfg, ax=ax)
    add_strategy_legend(ax, cfg)
    fig.savefig(OUT_DIR / "Figure1A_confirmed_origin_world.png", dpi=cfg.output_dpi, bbox_inches="tight")
    fig.savefig(OUT_DIR / "Figure1A_confirmed_origin_world.pdf", bbox_inches="tight")
    plt.close(fig)

    # Europe + Turkey only, no Russia
    fig, ax = plt.subplots(figsize=(8, 7))
    plot_europe_country_pies(by_country, cfg, ax=ax)
    add_strategy_legend(ax, cfg)
    fig.savefig(OUT_DIR / "Figure1B_confirmed_origin_Europe.png", dpi=cfg.output_dpi, bbox_inches="tight")
    fig.savefig(OUT_DIR / "Figure1B_confirmed_origin_Europe.pdf", bbox_inches="tight")
    plt.close(fig)

    # Combined panel
    fig, axes = plt.subplots(1, 2, figsize=(17, 7), gridspec_kw={"width_ratios": [1.45, 1.0]})
    plot_world_continent_pies(by_continent, cfg, ax=axes[0])
    plot_europe_country_pies(by_country, cfg, ax=axes[1])
    add_strategy_legend(axes[1], cfg)
    fig.savefig(OUT_DIR / "Figure1_REVISED_confirmed_origin.png", dpi=cfg.output_dpi, bbox_inches="tight")
    fig.savefig(OUT_DIR / "Figure1_REVISED_confirmed_origin.pdf", bbox_inches="tight")
    plt.close(fig)

    print("Saved outputs to:")
    print(OUT_DIR)
    print("- Figure1A_confirmed_origin_world.png/pdf")
    print("- Figure1B_confirmed_origin_Europe.png/pdf")
    print("- Figure1_REVISED_confirmed_origin.png/pdf")
    print("- Source_Data_Figure1A_continent_strategy_BioProjects.csv")
    print("- Source_Data_Figure1B_country_strategy_BioProjects.csv")


if __name__ == "__main__":
    main()

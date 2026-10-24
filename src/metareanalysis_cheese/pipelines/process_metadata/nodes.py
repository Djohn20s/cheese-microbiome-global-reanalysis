from __future__ import annotations

import logging
import re
from dataclasses import dataclass
import pandas as pd
from pathlib import Path


logger = logging.getLogger(__name__)


def _norm_text(x) -> str:
    if x is None:
        return ""
    s = str(x).strip()
    s = re.sub(r"\s+", " ", s)
    return s


def _norm_lower(x) -> str:
    return _norm_text(x).lower()


def _existing_fields(df: pd.DataFrame, fields: list[str]) -> list[str]:
    return [f for f in fields if f in df.columns]


def _build_text_blob(df: pd.DataFrame, fields: list[str]) -> pd.Series:
    use_fields = _existing_fields(df, fields)
    if not use_fields:
        return pd.Series([""] * len(df), index=df.index)
    parts = [df[f].fillna("").astype(str) for f in use_fields]
    blob = parts[0]
    for p in parts[1:]:
        blob = blob + " | " + p
    return blob.str.lower()


@dataclass(frozen=True)
class Stage1Config:
    whitelist: list[str]
    blacklist_regex: str
    text_fields: list[str]
    keep_nonwhitelisted: bool = False
    drop_unknown_continent: bool = False


def stage1_runlevel_curation(
    runinfo: pd.DataFrame,
    config: dict,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Returns:
      kept_df, dropped_df(with drop_reason), drop_reason_counts_df
    """
    cfg = Stage1Config(
        whitelist=[_norm_lower(x) for x in config["scientificname_whitelist"]],
        blacklist_regex=config["blacklist_regex"],
        text_fields=list(config["text_fields"]),
        keep_nonwhitelisted=bool(config.get("keep_nonwhitelisted", False)),
        drop_unknown_continent=bool(config.get("drop_unknown_continent", False)),
    )

    required_cols = ["ScientificName"]
    missing = [c for c in required_cols if c not in runinfo.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # compile blacklist
    try:
        bl_re = re.compile(cfg.blacklist_regex, flags=re.IGNORECASE)
    except re.error as e:
        raise ValueError(f"Invalid blacklist regex: {e}\n{cfg.blacklist_regex}") from e

    df = runinfo.copy()
    n0 = len(df)
    drop_reason = pd.Series([""] * len(df), index=df.index)

    # 1) whitelist
    if cfg.keep_nonwhitelisted:
        df_wl = df
        logger.info("Stage1 whitelist: DISABLED (keep_nonwhitelisted=true)")
    else:
        wl_set = set(cfg.whitelist)
        sci = df["ScientificName"].fillna("").map(_norm_lower)
        mask_wl = sci.isin(wl_set)
        drop_reason.loc[~mask_wl] = "not_whitelisted_scientificname"
        df_wl = df[mask_wl].copy()
        logger.info("Stage1 whitelist: kept=%s dropped=%s (from %s)", len(df_wl), n0 - len(df_wl), n0)

    # 2) blacklist keyword scan
    blob = _build_text_blob(df_wl, cfg.text_fields)
    hit_bl = blob.apply(lambda s: bool(bl_re.search(s)))
    df_keep = df_wl[~hit_bl].copy()
    df_drop_bl = df_wl[hit_bl].copy()
    if len(df_drop_bl) > 0:
        drop_reason.loc[df_drop_bl.index] = "blacklist_keyword_match"
    logger.info("Stage1 blacklist: kept=%s dropped=%s", len(df_keep), len(df_drop_bl))

    # 3) optional continent drop
    if cfg.drop_unknown_continent:
        if "Continent" not in df_keep.columns:
            logger.warning("Stage1: drop_unknown_continent requested but Continent column missing")
        else:
            cont = df_keep["Continent"].fillna("").astype(str).str.strip().str.lower()
            is_unknown = cont.isin({"", "na", "nan", "unknown_or_missing", "unknown"})
            df_drop_cont = df_keep[is_unknown].copy()
            df_keep = df_keep[~is_unknown].copy()
            if len(df_drop_cont) > 0:
                drop_reason.loc[df_drop_cont.index] = "unknown_or_missing_continent"
            logger.info("Stage1 continent: dropped=%s", len(df_drop_cont))

    df_dropped = df[drop_reason != ""].copy()
    df_dropped["drop_reason"] = drop_reason.loc[df_dropped.index].values

    drop_counts = (
        df_dropped["drop_reason"]
        .fillna("")
        .astype(str)
        .value_counts()
        .rename_axis("drop_reason")
        .reset_index(name="count")
    )

    logger.info("Stage1 summary: input=%s kept=%s dropped=%s", n0, len(df_keep), len(df_dropped))
    return df_keep, df_dropped, drop_counts


def summarize_runinfo_markdown(
    runinfo: pd.DataFrame,
    params: dict,
) -> str:
    """
    Kedro-friendly summary to Markdown string (no file IO).
    """
    from collections import Counter
    from datetime import datetime, timezone

    def safe_int(x) -> int:
        try:
            return int(float(x))
        except Exception:
            return 0

    def md_table(counter: Counter, top: int = 20) -> str:
        lines = ["| Value | Count |", "|---|---:|"]
        for k, v in counter.most_common(top):
            k_disp = str(k) if str(k).strip() else "(empty)"
            lines.append(f"| {k_disp} | {v} |")
        return "\n".join(lines)

    def quantiles(values: list[int]) -> dict:
        if not values:
            return {}
        v = sorted(values)
        n = len(v)

        def q(p):
            idx = int(round((n - 1) * p))
            return v[idx]

        return {
            "min": v[0],
            "p25": q(0.25),
            "median": q(0.50),
            "p75": q(0.75),
            "p90": q(0.90),
            "max": v[-1],
            "n": n,
        }

    topn = int(params.get("top", 20))
    dataset_prefix = str(params.get("dataset_prefix", "cheese_metagenome"))
    version = str(params.get("version", "v?"))

    # Counters
    c_strategy = Counter()
    c_layout = Counter()
    c_platform = Counter()
    c_model = Counter()
    c_source = Counter()
    c_library_source = Counter()
    c_scientific = Counter()
    c_center = Counter()
    c_continent = Counter()

    runs_per_bp = Counter()
    spots_list: list[int] = []
    bases_list: list[int] = []

    has_continent_col = "Continent" in runinfo.columns

    for _, row in runinfo.iterrows():
        strategy = _norm_text(row.get("LibraryStrategy"))
        layout = _norm_text(row.get("LibraryLayout"))
        platform = _norm_text(row.get("Platform"))
        model = _norm_text(row.get("Model"))
        lib_source = _norm_text(row.get("LibrarySource"))
        source = _norm_text(row.get("source"))
        sci = _norm_text(row.get("ScientificName"))
        center = _norm_text(row.get("CenterName"))
        bp = _norm_text(row.get("BioProject"))

        c_strategy[strategy] += 1
        c_layout[layout] += 1
        c_platform[platform] += 1
        c_model[model] += 1
        c_library_source[lib_source] += 1
        c_source[source] += 1
        c_scientific[sci] += 1
        c_center[center] += 1

        if has_continent_col:
            cont = _norm_text(row.get("Continent")).strip()
            if cont and cont.lower() not in {"unknown", "unknown_or_missing"}:
                c_continent[cont] += 1
            else:
                c_continent["unknown_or_missing"] += 1

        if bp:
            runs_per_bp[bp] += 1

        spots_list.append(safe_int(row.get("spots")))
        bases_list.append(safe_int(row.get("bases")))

    bp_q = quantiles(list(runs_per_bp.values()))
    spots_q = quantiles([x for x in spots_list if x > 0])
    bases_q = quantiles([x for x in bases_list if x > 0])

    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    md: list[str] = []
    md.append(f"# Runinfo summary – {dataset_prefix} ({version})\n")
    md.append(f"**Generated:** {now_utc} UTC\n")

    md.append("## Overall counts\n")
    md.append(f"- Total runs (SRR rows): **{len(runinfo)}**\n")
    md.append(f"- Unique BioProjects: **{len(runs_per_bp)}**\n")
    if has_continent_col:
        md.append(f"- Unknown/missing `Continent`: **{c_continent.get('unknown_or_missing', 0)}**\n")
    else:
        md.append("- `Continent` column: **NOT FOUND**\n")

    md.append("\n## LibraryStrategy (top)\n")
    md.append(md_table(c_strategy, topn) + "\n")

    md.append("## LibraryLayout (top)\n")
    md.append(md_table(c_layout, topn) + "\n")

    md.append("## Platform (top)\n")
    md.append(md_table(c_platform, topn) + "\n")

    md.append("## Model (top)\n")
    md.append(md_table(c_model, topn) + "\n")

    md.append("## LibrarySource (top)\n")
    md.append(md_table(c_library_source, topn) + "\n")

    md.append("## `source` field (top)\n")
    md.append(md_table(c_source, topn) + "\n")

    md.append("## ScientificName (top)\n")
    md.append(md_table(c_scientific, topn) + "\n")

    md.append("## CenterName (top)\n")
    md.append(md_table(c_center, topn) + "\n")

    if has_continent_col:
        md.append("## Continent (top)\n")
        md.append(md_table(c_continent, topn) + "\n")

    md.append("## BioProject size (runs per BioProject)\n")
    if bp_q:
        md.append(
            "- Quantiles (runs per BioProject):\n"
            f"  - n: {bp_q['n']}\n"
            f"  - min: {bp_q['min']}\n"
            f"  - p25: {bp_q['p25']}\n"
            f"  - median: {bp_q['median']}\n"
            f"  - p75: {bp_q['p75']}\n"
            f"  - p90: {bp_q['p90']}\n"
            f"  - max: {bp_q['max']}\n"
        )

    md.append("\n## Sequencing depth (spots, bases)\n")
    if spots_q:
        md.append(
            "- `spots` quantiles (non-zero):\n"
            f"  - n: {spots_q['n']}\n"
            f"  - min: {spots_q['min']}\n"
            f"  - median: {spots_q['median']}\n"
            f"  - p90: {spots_q['p90']}\n"
            f"  - max: {spots_q['max']}\n"
        )
    if bases_q:
        md.append(
            "\n- `bases` quantiles (non-zero):\n"
            f"  - n: {bases_q['n']}\n"
            f"  - min: {bases_q['min']}\n"
            f"  - median: {bases_q['median']}\n"
            f"  - p90: {bases_q['p90']}\n"
            f"  - max: {bases_q['max']}\n"
        )

    md.append("\n## Notes\n")
    md.append("- Descriptive report to guide filtering decisions.\n")

    return "\n".join(md)


def _decide_project(n_runs: int, frac_whitelist: float, keep_thr: float, drop_thr: float, min_runs: int) -> str:
    if n_runs < min_runs:
        return "MANUAL"
    if frac_whitelist >= keep_thr:
        return "KEEP"
    if frac_whitelist <= drop_thr:
        return "DROP"
    return "MANUAL"


def stage3_bioproject_curation(
    stage1_kept: pd.DataFrame,
    params: dict,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Applies:
      - require LibrarySource
      - allowed strategies
      - per-BioProject decision based on frac whitelisted ScientificName
    Returns:
      decisions_df, kept_df, dropped_df, manual_df
    """
    required = ["BioProject", "ScientificName", "LibrarySource", "LibraryStrategy"]
    missing = [c for c in required if c not in stage1_kept.columns]
    if missing:
        raise ValueError(f"Missing required columns for stage3: {missing}")

    whitelist = [_norm_lower(x) for x in params["scientificname_whitelist"]]
    wl_set = set(whitelist)

    req_src = str(params.get("require_librarysource", "METAGENOMIC")).strip().upper()
    allowed_strategy = {str(x).strip().upper() for x in params.get("allowed_strategy", ["AMPLICON", "WGS"])}

    keep_thr = float(params.get("keep_thr", 0.80))
    drop_thr = float(params.get("drop_thr", 0.20))
    min_runs = int(params.get("min_runs", 5))

    df = stage1_kept.copy()

    src = df["LibrarySource"].fillna("").astype(str).str.strip().str.upper()
    strat = df["LibraryStrategy"].fillna("").astype(str).str.strip().str.upper()

    mask = (src == req_src) & (strat.isin(allowed_strategy))
    df = df[mask].copy()

    if len(df) == 0:
        raise ValueError("Stage3 prefilter produced 0 rows (check require_librarysource/allowed_strategy)")

    df["_is_whitelisted_sciname"] = df["ScientificName"].fillna("").map(_norm_lower).isin(wl_set)

    g = df.groupby("BioProject", dropna=False)
    decisions = g.agg(
        runs=("BioProject", "size"),
        whitelisted_runs=("_is_whitelisted_sciname", "sum"),
    ).reset_index()
    decisions["frac_whitelisted"] = decisions["whitelisted_runs"] / decisions["runs"]
    decisions["decision"] = decisions.apply(
        lambda r: _decide_project(
            n_runs=int(r["runs"]),
            frac_whitelist=float(r["frac_whitelisted"]),
            keep_thr=keep_thr,
            drop_thr=drop_thr,
            min_runs=min_runs,
        ),
        axis=1,
    )
    decisions["reason"] = decisions.apply(
        lambda r: (
            "too_few_runs"
            if int(r["runs"]) < min_runs
            else "high_frac_whitelisted"
            if float(r["frac_whitelisted"]) >= keep_thr
            else "low_frac_whitelisted"
            if float(r["frac_whitelisted"]) <= drop_thr
            else "ambiguous"
        ),
        axis=1,
    )

    decision_map = dict(zip(decisions["BioProject"], decisions["decision"]))
    df["bioproject_decision"] = df["BioProject"].map(decision_map).fillna("MANUAL")

    kept = df[df["bioproject_decision"] == "KEEP"].copy()
    dropped = df[df["bioproject_decision"] == "DROP"].copy()
    manual = df[df["bioproject_decision"] == "MANUAL"].copy()

    for dfx in (kept, dropped, manual):
        dfx.drop(columns=["_is_whitelisted_sciname"], inplace=True, errors="ignore")

    return decisions, kept, dropped, manual


def make_manual_bioproject_template(stage3_manual: pd.DataFrame) -> pd.DataFrame:
    """
    Returns TSV template with BioProject + decision (OK/NOPE).
    """
    if "BioProject" not in stage3_manual.columns:
        raise ValueError("BioProject missing in stage3_manual")
    uniq = (
        stage3_manual["BioProject"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .sort_values()
        .reset_index(drop=True)
    )
    return pd.DataFrame({"BioProject": uniq, "decision": [""] * len(uniq)})


def apply_manual_and_merge_kept(
    stage1_kept_full: pd.DataFrame,
    stage3_auto_kept: pd.DataFrame,
    manual_decisions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Outputs:
      final_kept_runs_df (from stage1 full, filtered by union(auto_kept + manual_ok))
      final_bioprojects_df
    """
    for df, name in [
        (stage1_kept_full, "stage1_kept_full"),
        (stage3_auto_kept, "stage3_auto_kept"),
        (manual_decisions, "manual_decisions"),
    ]:
        if "BioProject" not in df.columns:
            raise ValueError(f"BioProject missing in {name}")

    if not {"BioProject", "decision"}.issubset(manual_decisions.columns):
        raise ValueError("manual_decisions must include BioProject + decision")

    auto_bps = set(stage3_auto_kept["BioProject"].dropna().astype(str).unique())

    dec = manual_decisions.copy()
    dec["decision_norm"] = dec["decision"].fillna("").astype(str).str.strip().str.lower()
    ok_bps = set(dec[dec["decision_norm"] == "ok"]["BioProject"].dropna().astype(str).unique())

    final_bps = auto_bps | ok_bps

    df = stage1_kept_full.copy()
    df["_bp"] = df["BioProject"].fillna("").astype(str)
    df_final = df[df["_bp"].isin(final_bps)].copy()
    df_final.drop(columns=["_bp"], inplace=True)

    df_final["bioproject_final_decision"] = "KEEP"
    df_final["bioproject_decision_source"] = df_final["BioProject"].apply(
        lambda x: "AUTO_KEEP" if str(x) in auto_bps else "MANUAL_OK"
    )

    final_bps_df = pd.DataFrame({"BioProject": sorted(final_bps)})

    return df_final, final_bps_df


def split_strategy(final_kept: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Returns: wgs_df, amplicon_df, others_df, strategy_counts_df
    """
    if "LibraryStrategy" not in final_kept.columns:
        raise ValueError("LibraryStrategy missing in final_kept")

    strat = final_kept["LibraryStrategy"].fillna("").astype(str).str.strip().str.upper()
    wgs = final_kept[strat == "WGS"].copy()
    amp = final_kept[strat == "AMPLICON"].copy()
    oth = final_kept[~strat.isin({"WGS", "AMPLICON"})].copy()

    counts = (
        strat.replace("", "(empty)")
        .value_counts()
        .rename_axis("LibraryStrategy")
        .reset_index(name="count")
    )
    return wgs, amp, oth, counts


def _norm_center(x) -> str:
    if x is None:
        return ""
    return str(x).strip()


def add_continent_by_centername(
    runinfo: pd.DataFrame,
    center_map: pd.DataFrame,
) -> pd.DataFrame:
    """
    Left-join: runinfo.CenterName -> center_map.CenterName and adds/overwrites Continent.
    center_map must have columns: CenterName, Continent (Notes optional).
    """
    if "CenterName" not in runinfo.columns:
        raise ValueError("runinfo missing column: CenterName")

    required_map_cols = {"CenterName", "Continent"}
    if not required_map_cols.issubset(center_map.columns):
        missing = required_map_cols - set(center_map.columns)
        raise ValueError(f"center_map missing columns: {sorted(missing)}")

    df = runinfo.copy()
    m = center_map.copy()

    df["CenterName"] = df["CenterName"].map(_norm_center)
    m["CenterName"] = m["CenterName"].map(_norm_center)

    m["Continent"] = m["Continent"].fillna("").astype(str).str.strip()

    m["_is_nonempty"] = (m["Continent"].str.len() > 0).astype(int)
    m = (
        m.sort_values(["CenterName", "_is_nonempty"], ascending=[True, False])
        .drop_duplicates(subset=["CenterName"], keep="first")
        .drop(columns=["_is_nonempty"], errors="ignore")
    )

    out = df.merge(
        m[["CenterName", "Continent"]],
        on="CenterName",
        how="left",
        suffixes=("", "_mapped"),
    )

    if "Continent_mapped" in out.columns:
        if "Continent" in df.columns:
            old = out["Continent"].fillna("").astype(str).str.strip()
            mapped = out["Continent_mapped"].fillna("").astype(str).str.strip()
            out["Continent"] = mapped.where(mapped != "", old)
        else:
            out["Continent"] = out["Continent_mapped"]

        out.drop(columns=["Continent_mapped"], inplace=True, errors="ignore")

    cont = out["Continent"].fillna("").astype(str).str.strip()
    unknown = (cont == "").sum()
    logger.info("Continent mapping: total_runs=%s unknown_continent=%s", len(out), int(unknown))

    return out

def require_and_load_manual_bioprojects_filled(
    manual_template: pd.DataFrame,
    manual_filled_path: str,
    manual_template_path: str,
) -> pd.DataFrame:
    """
    Guard node:
    - TEMPLATE jest generowany i zapisywany jako dataset (manual_template DF jest tylko "anchor" w grafie).
    - Jeśli FILLED nie istnieje -> STOP z jasnym komunikatem.
    - Jeśli istnieje -> wczytaj i zwróć DF decyzji + walidacja wartości.
    """
    filled = Path(manual_filled_path)

    allowed_values = ["OK", "NOPE"]  # <- jawnie komunikujemy i walidujemy
    allowed_values_lower = {v.lower() for v in allowed_values}

    if not filled.exists():
        msg = (
            "\n[STOP] Missing manual BioProject decisions file.\n"
            f"Expected FILLED TSV at:\n  {filled}\n\n"
            "A TEMPLATE TSV has been generated. Please fill it and save as FILLED.\n"
            f"Template TSV location:\n  {manual_template_path}\n\n"
            "How to fill:\n"
            "- Open TEMPLATE TSV\n"
            "- For each BioProject set column 'decision' to one of:\n"
            f"    - {allowed_values[0]}  (keep BioProject)\n"
            f"    - {allowed_values[1]} (drop BioProject)\n"
            "- Leave no empty decisions for projects you want to decide now\n"
            "- Save as TSV with columns: BioProject, decision\n\n"
            "Then resume with:\n"
            "  kedro run --pipeline process_metadata --from-nodes stage3_apply_manual_and_merge\n"
        )
        raise FileNotFoundError(msg)

    df = pd.read_csv(filled, sep="\t", dtype=str)

    required = {"BioProject", "decision"}
    if not required.issubset(df.columns):
        raise ValueError(
            f"[STOP] FILLED TSV invalid. Required columns {sorted(required)} "
            f"but got {list(df.columns)}. Path: {filled}"
        )

    df = df[["BioProject", "decision"]].copy()
    df["BioProject"] = df["BioProject"].fillna("").astype(str).str.strip()
    df["decision"] = df["decision"].fillna("").astype(str).str.strip()

    # sanity
    empties_bp = (df["BioProject"] == "").sum()
    if empties_bp:
        raise ValueError(f"[STOP] FILLED TSV has empty BioProject rows: {int(empties_bp)}. Path: {filled}")

    # validate decision values (ignore empty lines if you want; tu: wymagamy poprawnych, niepustych)
    dnorm = df["decision"].str.lower()
    bad_mask = ~dnorm.isin(allowed_values_lower)
    if bad_mask.any():
        bad = (
            df.loc[bad_mask, ["BioProject", "decision"]]
            .head(20)
            .to_string(index=False)
        )
        raise ValueError(
            "\n[STOP] FILLED TSV contains invalid values in column 'decision'.\n"
            f"Allowed values: {', '.join(allowed_values)}\n"
            f"File: {filled}\n\n"
            "Examples of invalid rows (first 20):\n"
            f"{bad}\n"
        )

    return df
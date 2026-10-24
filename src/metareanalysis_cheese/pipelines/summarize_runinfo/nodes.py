from __future__ import annotations

import pandas as pd


def summarize_stage5_metadata(
    wgs: pd.DataFrame,
    amplicon: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    """
    Build summary tables:
    A) Continent x Strategy (unique BioProjects)
    B) Strategy totals (unique BioProjects)
    C) Year x Strategy (unique BioProjects)
    D) Year x Strategy x Platform (RUN COUNTS; Platform = sequencer platform)

    Returns dict[str, DataFrame] suitable for a PartitionedDataset.
    """

    def _prep_runs(df: pd.DataFrame, strategy: str) -> pd.DataFrame:
        """Keep run-level rows (do NOT collapse to BioProject)."""
        if df is None or df.empty:
            return pd.DataFrame(columns=["BioProject", "Continent", "ReleaseDate", "year", "strategy", "Platform"])

        out = df.copy()

        # Tolerate missing columns
        for col in ["BioProject", "Continent", "ReleaseDate", "Platform"]:
            if col not in out.columns:
                out[col] = pd.NA

        out["strategy"] = strategy
        out["BioProject"] = out["BioProject"].astype("string")
        out["Continent"] = out["Continent"].fillna("Unknown").astype("string")
        out["Platform"] = out["Platform"].fillna("Unknown").astype("string")

        out["ReleaseDate"] = pd.to_datetime(out["ReleaseDate"], errors="coerce")
        out["year"] = out["ReleaseDate"].dt.year

        return out[["BioProject", "Continent", "ReleaseDate", "year", "strategy", "Platform"]]

    def _prep_projects(run_df: pd.DataFrame) -> pd.DataFrame:
        """Collapse run-level rows to one row per BioProject for project-level summaries."""
        if run_df is None or run_df.empty:
            return pd.DataFrame(columns=["BioProject", "Continent", "year", "strategy"])

        def _mode_or_unknown(s: pd.Series) -> str:
            s = s.dropna().astype(str)
            if s.empty:
                return "Unknown"
            m = s.mode()
            return str(m.iloc[0]) if not m.empty else "Unknown"

        proj = (
            run_df.groupby("BioProject", as_index=False)
            .agg(
                Continent=("Continent", _mode_or_unknown),
                year=("year", "min"),
                strategy=("strategy", "first"),
            )
        )
        proj["Continent"] = proj["Continent"].fillna("Unknown").astype("string")
        proj["strategy"] = proj["strategy"].astype("string")
        return proj

    # ---- build run-level tables (keeps Platform) ----
    runs_wgs = _prep_runs(wgs, "WGS")
    runs_amp = _prep_runs(amplicon, "AMPLICON")

    # overlap safety (shouldn't happen, but keep consistent)
    overlap = set(runs_wgs["BioProject"].dropna()).intersection(set(runs_amp["BioProject"].dropna()))
    if overlap:
        runs_amp = runs_amp[~runs_amp["BioProject"].isin(overlap)].copy()

    run_all = pd.concat([runs_wgs, runs_amp], ignore_index=True)

    # ---- build project-level tables (unique BioProjects) ----
    proj_all = _prep_projects(run_all)

    # A) continent x strategy (BioProjects)
    by_continent = (
        proj_all.groupby(["Continent", "strategy"], as_index=False)
        .agg(n_bioprojects=("BioProject", "nunique"))
    )
    totals_cont = by_continent.groupby("Continent")["n_bioprojects"].transform("sum")
    by_continent["pct_within_continent"] = (by_continent["n_bioprojects"] / totals_cont) * 100
    by_continent = by_continent.sort_values(["Continent", "strategy"]).reset_index(drop=True)

    # B) strategy totals (BioProjects)
    by_strategy = (
        proj_all.groupby("strategy", as_index=False)
        .agg(n_bioprojects=("BioProject", "nunique"))
        .sort_values("n_bioprojects", ascending=False)
        .reset_index(drop=True)
    )
    total_all = float(by_strategy["n_bioprojects"].sum())
    by_strategy["pct_total"] = (by_strategy["n_bioprojects"] / total_all) * 100 if total_all else 0.0

    # C) year x strategy (BioProjects)
    by_year_strategy = (
        proj_all.dropna(subset=["year"])
        .groupby(["year", "strategy"], as_index=False)
        .agg(n_bioprojects=("BioProject", "nunique"))
        .sort_values(["year", "strategy"])
        .reset_index(drop=True)
    )
    year_totals = by_year_strategy.groupby("year")["n_bioprojects"].transform("sum")
    by_year_strategy["pct_within_year"] = (by_year_strategy["n_bioprojects"] / year_totals) * 100

    # D) year x strategy x Platform (RUN COUNTS)
    by_year_strategy_platform = (
        run_all.dropna(subset=["year"])
        .assign(Platform=lambda d: d["Platform"].fillna("Unknown").astype("string"))
        .groupby(["year", "strategy", "Platform"], as_index=False)
        .agg(n_runs=("Platform", "size"))
        .sort_values(["year", "strategy", "n_runs"], ascending=[True, True, False])
        .reset_index(drop=True)
    )

    return {
        "by_continent": by_continent,
        "by_strategy": by_strategy,
        "by_year_strategy": by_year_strategy,
        "by_year_strategy_platform": by_year_strategy_platform,
    }
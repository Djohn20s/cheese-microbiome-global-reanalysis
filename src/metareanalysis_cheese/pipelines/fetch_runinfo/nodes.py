from __future__ import annotations

import io
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Tuple

import pandas as pd


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _now_fs() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")


def _run_bash_capture(cmd: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["bash", "-lc", cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def _extract_esearch_count(esearch_xml: str) -> int:
    for line in esearch_xml.splitlines():
        line = line.strip()
        if line.startswith("<Count>") and line.endswith("</Count>"):
            inner = line.replace("<Count>", "").replace("</Count>", "").strip()
            try:
                return int(inner)
            except Exception:
                return -1
    return -1


def _safe_unique_centers(runinfo_df: pd.DataFrame) -> list[str]:
    if "CenterName" not in runinfo_df.columns:
        return []
    centers = (
        runinfo_df["CenterName"]
        .fillna("")
        .astype(str)
        .str.strip()
    )
    centers = centers[centers != ""].drop_duplicates().sort_values()
    return centers.tolist()


def fetch_runinfo(ncbi: dict) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Fetch SRA RunInfo via EDirect:
      esearch -db sra -query "<QUERY>" | efetch -format runinfo

    Outputs (Kedro):
      - runinfo_raw: DataFrame (saved by Catalog)
      - runinfo_download_summary: dict (JSON)

    Side files (for reproducibility / manual step):
      - data/01_raw/ncbi/runinfo_<prefix>_<version>_<timestamp>.csv (snapshot)
      - data/01_raw/ncbi/queries/esearch_<prefix>_<version>_<timestamp>.query
      - data/01_raw/ncbi/meta/esearch_<prefix>_<version>_<timestamp>.meta.tsv
      - data/01_raw/ncbi/center_continent_map_TEMPLATE.tsv (template)
    """
    prefix = str(ncbi["dataset_prefix"]).strip()
    version = str(ncbi["version"]).strip()
    query = str(ncbi["query"]).strip()

    ts_utc = _now_utc()
    ts_fs = _now_fs()

    # output locations for logs/snapshots
    raw_dir = Path("data/01_raw/ncbi")
    q_dir = raw_dir / "queries"
    meta_dir = raw_dir / "meta"
    raw_dir.mkdir(parents=True, exist_ok=True)
    q_dir.mkdir(parents=True, exist_ok=True)
    meta_dir.mkdir(parents=True, exist_ok=True)

    out_runinfo_snapshot = raw_dir / f"runinfo_{prefix}_{version}_{ts_fs}.csv"
    out_query = q_dir / f"esearch_{prefix}_{version}_{ts_fs}.query"
    out_meta = meta_dir / f"esearch_{prefix}_{version}_{ts_fs}.meta.tsv"

    out_query.write_text(query + "\n", encoding="utf-8")

    # 1) count
    cmd_esearch = "esearch -db sra -query " + subprocess.list2cmdline([query])
    res_count = _run_bash_capture(cmd_esearch)
    if res_count.returncode != 0:
        raise RuntimeError(f"esearch failed (exit={res_count.returncode}): {res_count.stderr[:4000]}")
    count = _extract_esearch_count(res_count.stdout)

    # 2) runinfo
    cmd_pipe = (
        "esearch -db sra -query "
        + subprocess.list2cmdline([query])
        + " | efetch -format runinfo"
    )
    res_runinfo = _run_bash_capture(cmd_pipe)
    if res_runinfo.returncode != 0:
        raise RuntimeError(f"efetch runinfo failed (exit={res_runinfo.returncode}): {res_runinfo.stderr[:4000]}")

    content = res_runinfo.stdout or ""
    lines = content.splitlines()

    # sanity: header + at least one record
    if len(lines) < 2:
        preview = content[:500].replace("\n", "\\n")
        raise RuntimeError(
            f"RunInfo empty or header-only. lines={len(lines)} preview={preview} stderr={res_runinfo.stderr[:1000]}"
        )

    header_preview = lines[0][:200]
    n_lines = len(lines)

    # Parse to DataFrame -> goes to `runinfo_raw` via Catalog
    runinfo_df = pd.read_csv(io.StringIO(content), low_memory=False)

    # Snapshot for reproducibility
    out_runinfo_snapshot.write_text(content, encoding="utf-8")

    # Manual step: center -> continent map
    center_map_path = Path("data/01_raw/ncbi/center_continent_map.tsv")
    center_map_template_path = Path("data/01_raw/ncbi/center_continent_map_TEMPLATE.tsv")

    centers = _safe_unique_centers(runinfo_df)
    template_created = False
    if centers:
        tmpl = pd.DataFrame(
            {"CenterName": centers, "Continent": [""] * len(centers), "Notes": [""] * len(centers)}
        )
        tmpl.to_csv(center_map_template_path, sep="\t", index=False)
        template_created = True

    # meta tsv
    out_meta.write_text(
        "timestamp_utc\tsra_count\truninfo_snapshot\tquery_file\theader_preview\tlines\tcenter_map\tcenter_map_template\n"
        f"{ts_utc}\t{count}\t{out_runinfo_snapshot}\t{out_query}\t{header_preview}\t{n_lines}\t{center_map_path}\t{center_map_template_path}\n",
        encoding="utf-8",
    )

    summary: Dict[str, Any] = {
        "status": "OK",
        "timestamp_utc": ts_utc,
        "sra_count": int(count),
        "runinfo_rows": int(len(runinfo_df)),
        "runinfo_cols_n": int(len(runinfo_df.columns)),
        "header_preview": header_preview,
        "lines": int(n_lines),
        "runinfo_snapshot": str(out_runinfo_snapshot),
        "query_file": str(out_query),
        "meta_file": str(out_meta),
        "center_continent_map": str(center_map_path),
        "center_continent_map_template": str(center_map_template_path),
        "center_map_template_created": bool(template_created),
        "action_required": (
            "Uzupełnij plik center_continent_map.tsv: przypisz Continent do CenterName. "
            "Możesz skopiować center_continent_map_TEMPLATE.tsv -> center_continent_map.tsv "
            "i wypełnić kolumnę Continent."
        ),
        "why_required": (
            "Bez uzupełnienia mapy, add_continent_by_centername nie przypisze kontynentu "
            "i Continent będzie unknown_or_missing w downstream (summary/mapy/filtry)."
        ),
    }

    return runinfo_df, summary
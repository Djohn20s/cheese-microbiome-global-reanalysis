from __future__ import annotations

import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd


# ----------------------------
# shell helpers
# ----------------------------
def _run(cmd: str, cwd: Optional[Path] = None, timeout_sec: Optional[int] = None) -> subprocess.CompletedProcess:
    """
    Runs bash -lc <cmd>, captures stdout/stderr, optional timeout.
    If timeout triggers, raises subprocess.TimeoutExpired (caught upstream).
    """
    return subprocess.run(
        ["bash", "-lc", cmd],
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=timeout_sec,
    )


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _norm(x: Any) -> str:
    return re.sub(r"\s+", "", str(x).strip()) if x is not None else ""


def _pick_col(df: pd.DataFrame, candidates: List[str], df_name: str) -> str:
    col = next((c for c in candidates if c in df.columns), None)
    if col is None:
        raise ValueError(f"{df_name}: missing required column. Candidates={candidates}. Found={list(df.columns)}")
    return col


def _copytree_merge(src: Path, dst: Path) -> None:
    _ensure_dir(dst)
    for item in src.iterdir():
        s = item
        d = dst / item.name
        if s.is_dir():
            _copytree_merge(s, d)
        else:
            shutil.copy2(s, d)


def _fastq_files(fastq_dir: Path, run: str) -> List[Path]:
    return sorted(fastq_dir.glob(f"{run}*.fastq")) + sorted(fastq_dir.glob(f"{run}*.fastq.gz"))


def _validate_fastq_files(fqs: List[Path]) -> None:
    if not fqs:
        raise RuntimeError("No FASTQ files produced")
    for p in fqs:
        if not p.exists():
            raise RuntimeError(f"FASTQ missing: {p}")
        if p.stat().st_size == 0:
            raise RuntimeError(f"FASTQ empty: {p}")


def _looks_like_internal_mirror_error(text: str) -> bool:
    t = (text or "")
    needles = [
        "sra-download-internal.ncbi.nlm.nih.gov",
        "Could not resolve host",
        "rcConnection",
        "rcNotFound",
        "Cannot KNSManagerMakeClientRequest",
        "connection not found while validating within network system module",
    ]
    return any(n in t for n in needles)


# ----------------------------
# ENA HTTPS downloader
# ----------------------------
def _ena_fastq_https_urls(run: str, user_agent: str = "metareanalysis-cheese/1.0", timeout_sec: int = 60) -> List[str]:
    """
    Query ENA portal for FASTQ FTP paths and convert to HTTPS URLs.

    Returns list of URLs pointing to *.fastq.gz
    """
    cmd = (
        f"curl -s -A '{user_agent}' "
        f"\"https://www.ebi.ac.uk/ena/portal/api/filereport?"
        f"accession={run}&result=read_run&fields=fastq_ftp&format=tsv\""
    )
    res = _run(cmd, timeout_sec=timeout_sec)
    if res.returncode != 0:
        raise RuntimeError(f"ENA filereport failed: {res.stderr or res.stdout}")

    lines = [ln.strip() for ln in (res.stdout or "").splitlines() if ln.strip()]
    if len(lines) < 2:
        raise RuntimeError(f"ENA filereport empty for run={run}")

    # expected: header + one data line
    parts = lines[1].split("\t")
    if len(parts) < 2 or not parts[1].strip():
        raise RuntimeError(f"ENA filereport missing fastq_ftp for run={run}: {lines[1]}")

    ftp_field = parts[1].strip()
    urls: List[str] = []
    for item in ftp_field.split(";"):
        item = item.strip()
        if not item:
            continue

        # Output often looks like: ftp.sra.ebi.ac.uk/vol1/fastq/...
        if item.startswith("ftp.sra.ebi.ac.uk/"):
            urls.append("https://" + item)
        elif item.startswith("ftp://"):
            urls.append(item.replace("ftp://", "https://", 1))
        elif item.startswith("https://"):
            urls.append(item)
        else:
            # fallback: assume host/path without scheme
            urls.append("https://" + item)

    # de-dup, preserve order
    seen = set()
    out: List[str] = []
    for u in urls:
        if u not in seen:
            seen.add(u)
            out.append(u)
    return out


def _download_ena_https(
    run: str,
    fastq_dir: Path,
    retries: int = 5,
    timeout_sec: int = 600,
    user_agent: str = "metareanalysis-cheese/1.0",
) -> None:
    """
    Download .fastq.gz files via ENA HTTPS into fastq_dir.
    Uses wget -c for resume, per-file retries with backoff.
    """
    _ensure_dir(fastq_dir)
    urls = _ena_fastq_https_urls(run, user_agent=user_agent)
    if not urls:
        raise RuntimeError(f"No ENA FASTQ urls for run={run}")

    for url in urls:
        out_path = fastq_dir / Path(url).name
        ok = False
        last_err = ""

        for i in range(1, retries + 1):
            # -c resume, -T per-connection timeout, -t 1 single try per loop attempt
            cmd = f"wget -c -T 30 -t 1 -O '{out_path}' '{url}'"
            try:
                res = _run(cmd, timeout_sec=timeout_sec)
            except subprocess.TimeoutExpired:
                last_err = f"timeout after {timeout_sec}s"
                time.sleep(min(60 * i, 300))
                continue

            if res.returncode == 0 and out_path.exists() and out_path.stat().st_size > 0:
                ok = True
                break

            last_err = (res.stderr or res.stdout or "").strip()
            time.sleep(min(60 * i, 300))

        if not ok:
            raise RuntimeError(f"ENA download failed for {url}: {last_err[:500]}")


# ----------------------------
# public API
# ----------------------------
def build_runs_table(final_wgs: pd.DataFrame, final_amplicon: pd.DataFrame) -> pd.DataFrame:
    """Tworzy jedną tabelę runów z kolumną strategy (WGS/AMPLICON)."""

    def prep(df: pd.DataFrame, strategy: str) -> pd.DataFrame:
        if df is None or df.empty:
            return pd.DataFrame(columns=["Run", "strategy"])
        run_col = _pick_col(df, ["Run", "run", "SRR", "srr", "accession", "run_accession"], f"final_{strategy}")
        out = pd.DataFrame({"Run": df[run_col].astype(str), "strategy": strategy})
        # opcjonalnie dołącz BioProject/Continent jeśli istnieją (przyda się do logów)
        for col in ["BioProject", "Continent"]:
            if col in df.columns:
                out[col] = df[col].astype("string")
        return out.dropna(subset=["Run"]).drop_duplicates(subset=["Run"])

    wgs = prep(final_wgs, "WGS")
    amp = prep(final_amplicon, "AMPLICON")
    return pd.concat([wgs, amp], ignore_index=True)


def fetch_unpack_stage(
    runs_table: pd.DataFrame,
    params: Dict[str, Any],
) -> pd.DataFrame:
    """
    Pobiera SRA + robi FASTQ w RAM, potem przenosi FASTQ na HDD (staging) i zwalnia RAM.

    Output: manifest z kolumnami:
      Run, strategy, status, message, staged_fastq_dir, r1, r2, se,
      prefetch_sec, fasterq_sec, total_sec, ...
    """
    cfg = params

    max_runs_wgs = int(cfg.get("max_runs_wgs", 0))  # 0 = bez limitu
    max_runs_amplicon = int(cfg.get("max_runs_amplicon", 0))
    threads = int(cfg.get("threads", 8))

    ram_work_dir = cfg.get("ram_work_dir")  # np. /dev/shm/cheese_sra_work
    hdd_stage_dir = cfg["hdd_stage_dir"]  # e.g. work/sra_stage_fastq
    cleanup_ram = bool(cfg.get("cleanup_ram", True))

    # backend selection
    download_backend = str(cfg.get("download_backend", "sratools")).strip().lower()  # "sratools" | "ena_https"
    fallback_ena_https = bool(cfg.get("fallback_ena_https", False))

    # tool commands (only used in sratools mode)
    prefetch_cmd = str(cfg.get("prefetch_cmd", "prefetch {run} --output-directory '{sra_dir}'"))
    fasterq_cmd = str(
        cfg.get(
            "fasterq_cmd",
            "fasterq-dump {run} --split-files --threads {threads} --outdir '{fastq_dir}' -t '{tmp_dir}'",
        )
    )

    # timeouts (avoid hanging forever)
    prefetch_timeout_sec = int(cfg.get("prefetch_timeout_sec", 1800))  # 30 min
    fasterq_timeout_sec = int(cfg.get("fasterq_timeout_sec", 1800))  # 30 min
    ena_timeout_sec = int(cfg.get("ena_timeout_sec", 600))  # per file
    ena_retries = int(cfg.get("ena_retries", 5))
    ena_user_agent = str(cfg.get("ena_user_agent", "metareanalysis-cheese/1.0"))

    # sanity: tools (only if needed)
    if download_backend == "sratools" or fallback_ena_https:
        for tool in ["prefetch", "fasterq-dump"]:
            chk = _run(f"command -v {tool} >/dev/null 2>&1 && echo OK || echo MISSING")
            if "MISSING" in (chk.stdout or ""):
                raise RuntimeError(f"{tool} not found in PATH (activate env with sra-tools)")

    df = runs_table.copy()
    if df.empty:
        raise ValueError("runs_table is empty (nothing to fetch/unpack).")

    # apply per-strategy limits
    out_rows: List[Dict[str, Any]] = []
    wgs_df = df[df["strategy"] == "WGS"].copy()
    amp_df = df[df["strategy"] == "AMPLICON"].copy()

    if max_runs_wgs > 0:
        wgs_df = wgs_df.head(max_runs_wgs)
    if max_runs_amplicon > 0:
        amp_df = amp_df.head(max_runs_amplicon)

    df = pd.concat([wgs_df, amp_df], ignore_index=True)
    if df.empty:
        raise ValueError("After applying limits, runs_table is empty.")

    ram_base = Path(ram_work_dir) if ram_work_dir else None
    hdd_base = Path(hdd_stage_dir)
    _ensure_dir(hdd_base)

    for _, r in df.iterrows():
        run = str(r["Run"])
        strategy = str(r["strategy"])
        bioproject = str(r.get("BioProject", "")) if "BioProject" in r else ""
        continent = str(r.get("Continent", "")) if "Continent" in r else ""

        t0 = time.perf_counter()
        status = "OK"
        msg = ""
        prefetch_sec = 0.0
        fasterq_sec = 0.0

        staged_dir = hdd_base / strategy / run
        staged_fastq_dir = staged_dir / "fastq"
        _ensure_dir(staged_fastq_dir)

        # RAM dirs
        if ram_base:
            ram_run_dir = ram_base / strategy / run
        else:
            ram_run_dir = Path(cfg.get("fallback_work_dir", "data/07_model_output/sra_fetch_unpack/_work")) / strategy / run

        sra_dir = ram_run_dir / "sra"
        fastq_dir = ram_run_dir / "fastq"
        tmp_dir = ram_run_dir / "tmp"
        _ensure_dir(sra_dir)
        _ensure_dir(fastq_dir)
        _ensure_dir(tmp_dir)

        r1_out = ""
        r2_out = ""
        se_out = ""

        try:
            # ----------------------------
            # DOWNLOAD / UNPACK
            # ----------------------------
            if download_backend == "ena_https":
                # ENA delivers gz FASTQ directly
                t = time.perf_counter()
                _download_ena_https(
                    run=run,
                    fastq_dir=fastq_dir,
                    retries=ena_retries,
                    timeout_sec=ena_timeout_sec,
                    user_agent=ena_user_agent,
                )
                fasterq_sec = time.perf_counter() - t

            else:
                # 1) prefetch (RAM)
                cmd = prefetch_cmd.format(run=run, sra_dir=sra_dir)
                t = time.perf_counter()
                try:
                    res = _run(cmd, timeout_sec=prefetch_timeout_sec)
                except subprocess.TimeoutExpired:
                    raise RuntimeError(f"prefetch timeout after {prefetch_timeout_sec}s")
                prefetch_sec = time.perf_counter() - t
                if res.returncode != 0:
                    err = res.stderr or res.stdout or ""
                    if fallback_ena_https and _looks_like_internal_mirror_error(err):
                        _download_ena_https(
                            run=run,
                            fastq_dir=fastq_dir,
                            retries=ena_retries,
                            timeout_sec=ena_timeout_sec,
                            user_agent=ena_user_agent,
                        )
                    else:
                        raise RuntimeError(f"prefetch failed: {err}")

                else:
                    # 2) fasterq-dump (RAM)
                    cmd = fasterq_cmd.format(run=run, fastq_dir=fastq_dir, tmp_dir=tmp_dir, threads=threads)
                    t = time.perf_counter()
                    try:
                        res = _run(cmd, timeout_sec=fasterq_timeout_sec)
                    except subprocess.TimeoutExpired:
                        raise RuntimeError(f"fasterq-dump timeout after {fasterq_timeout_sec}s")
                    fasterq_sec = time.perf_counter() - t
                    if res.returncode != 0:
                        err = res.stderr or res.stdout or ""
                        if fallback_ena_https and _looks_like_internal_mirror_error(err):
                            _download_ena_https(
                                run=run,
                                fastq_dir=fastq_dir,
                                retries=ena_retries,
                                timeout_sec=ena_timeout_sec,
                                user_agent=ena_user_agent,
                            )
                        else:
                            raise RuntimeError(f"fasterq-dump failed: {err}")

            # ----------------------------
            # VALIDATE + STAGE TO HDD
            # ----------------------------
            fqs = _fastq_files(fastq_dir, run)
            _validate_fastq_files(fqs)

            # 3) stage FASTQ -> HDD
            _copytree_merge(fastq_dir, staged_fastq_dir)

            # detect R1/R2 in staged dir
            r1 = sorted(staged_fastq_dir.glob(f"{run}*_1.fastq")) + sorted(staged_fastq_dir.glob(f"{run}*_1.fastq.gz"))
            r2 = sorted(staged_fastq_dir.glob(f"{run}*_2.fastq")) + sorted(staged_fastq_dir.glob(f"{run}*_2.fastq.gz"))
            r1_out = str(r1[0]) if r1 else ""
            r2_out = str(r2[0]) if r2 else ""

            # detect SE (ignore _1/_2)
            se = sorted(staged_fastq_dir.glob(f"{run}*.fastq")) + sorted(staged_fastq_dir.glob(f"{run}*.fastq.gz"))
            se = [p for p in se if not re.search(r"_[12]\.fastq(\.gz)?$", p.name)]
            se_out = str(se[0]) if se else ""

            # 4) cleanup RAM to free space
            if cleanup_ram:
                shutil.rmtree(ram_run_dir, ignore_errors=True)

        except Exception as e:
            status = "FAILED"
            msg = str(e).replace("\n", " ")[:2000]

        total_sec = time.perf_counter() - t0

        out = {
            "Run": run,
            "strategy": strategy,
            "status": status,
            "message": msg,
            "staged_fastq_dir": str(staged_fastq_dir),
            "r1": r1_out,
            "r2": r2_out,
            "se": se_out,
            "prefetch_sec": round(prefetch_sec, 3),
            "fasterq_sec": round(fasterq_sec, 3),
            "total_sec": round(total_sec, 3),
            "BioProject": bioproject,
            "Continent": continent,
        }
        out_rows.append(out)

    return pd.DataFrame(out_rows)
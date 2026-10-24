from __future__ import annotations

import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd


def _run(cmd: str, cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["bash", "-lc", cmd],
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _write_text(path: Path, text: str) -> None:
    _ensure_dir(path.parent)
    path.write_text(text or "", encoding="utf-8", errors="replace")


def _copy_if_exists(src: Path, dst: Path) -> None:
    try:
        if src.exists():
            _ensure_dir(dst.parent)
            shutil.copy2(src, dst)
    except Exception:
        pass


def _pick_fastq_mode_and_paths(row: pd.Series) -> Tuple[str, str, str, str]:
    fastq_dir = Path(str(row["staged_fastq_dir"]))

    r1 = str(row.get("r1") or "").strip()
    r2 = str(row.get("r2") or "").strip()
    se = str(row.get("se") or "").strip()

    if r1 and r2:
        return "PE", r1, r2, ""
    if se:
        return "SE", "", "", se

    run = str(row["Run"])
    r1_c = sorted(fastq_dir.glob(f"{run}*_1.fastq")) + sorted(fastq_dir.glob(f"{run}*_1.fastq.gz"))
    r2_c = sorted(fastq_dir.glob(f"{run}*_2.fastq")) + sorted(fastq_dir.glob(f"{run}*_2.fastq.gz"))
    if r1_c and r2_c:
        return "PE", str(r1_c[0]), str(r2_c[0]), ""

    se_c = sorted(fastq_dir.glob(f"{run}*.fastq")) + sorted(fastq_dir.glob(f"{run}*.fastq.gz"))
    se_c = [p for p in se_c if not re.search(r"_[12]\.fastq(\.gz)?$", p.name)]
    if se_c:
        return "SE", "", "", str(se_c[0])

    raise RuntimeError(f"Cannot detect FASTQ inputs for run={run} in {fastq_dir}")


def parse_kraken_report(report_path: Path, run: str) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    rows = []
    total_reads = 0
    unclassified = 0

    with report_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                parts = re.split(r"\s+", line.strip(), maxsplit=5)
            if len(parts) < 6:
                continue

            pct = float(parts[0])
            clade_reads = int(parts[1])
            direct_reads = int(parts[2])
            rank = parts[3]
            taxid = parts[4]
            depth = len(re.match(r"^\s*", parts[5]).group(0)) if parts[5] else 0
            name = parts[5].strip()

            rows.append(
                dict(
                    Run=run,
                    percent=pct,
                    clade_reads=clade_reads,
                    direct_reads=direct_reads,
                    rank=rank,
                    taxid=str(taxid),
                    name=name,
                    depth=depth,
                )
            )

            if rank == "U" and str(taxid) == "0":
                unclassified = clade_reads
            if rank == "R" and str(taxid) == "1":
                total_reads = clade_reads

    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError(f"Kraken report parsed empty: {report_path}")

    if total_reads == 0:
        total_reads = int(df["clade_reads"].max())

    classified = max(total_reads - unclassified, 0)
    summary = dict(
        kraken_report=str(report_path),
        kraken_total_reads=int(total_reads),
        kraken_unclassified_reads=int(unclassified),
        kraken_classified_reads=int(classified),
    )
    return df, summary


def parse_motus_tsv(motus_path: Path, run: str) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    df = pd.read_csv(motus_path, sep="\t", comment="#", dtype=str, engine="python")
    if df.empty:
        raise RuntimeError(f"mOTUs TSV empty: {motus_path}")

    cols = list(df.columns)
    if "mOTU" not in cols:
        raise RuntimeError(f"mOTUs TSV missing mOTU column: {motus_path} (cols={cols})")
    if "Taxonomy" not in cols:
        raise RuntimeError(f"mOTUs TSV missing Taxonomy column: {motus_path} (cols={cols})")

    count_col = run if run in cols else cols[-1]
    df = df.rename(columns={count_col: "count"}).copy()
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0).astype(float)

    long_df = df[["mOTU", "Taxonomy", "count"]].copy()
    long_df.insert(0, "Run", run)

    summary = dict(
        motus_tsv=str(motus_path),
        motus_nonzero_motus=int((long_df["count"] > 0).sum()),
        motus_total_count=float(long_df["count"].sum()),
    )
    return long_df, summary


def _best_kraken_report(out_dir: Path) -> Optional[Path]:
    cands = list(out_dir.glob("*.kreport")) + list(out_dir.rglob("*.kreport"))
    if not cands:
        cands = list(out_dir.rglob("*.report")) + list(out_dir.rglob("*report*.txt"))
    if not cands:
        return None
    cands.sort(key=lambda p: p.stat().st_size if p.exists() else 0, reverse=True)
    return cands[0]


def _load_matrix(path: Path, key_col: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=[key_col])
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return pd.DataFrame(columns=[key_col])
        if key_col not in df.columns:
            df = df.rename(columns={df.columns[0]: key_col})
        return df
    except Exception:
        return pd.DataFrame(columns=[key_col])


def _merge_matrix(existing: pd.DataFrame, new_wide: pd.DataFrame, key_col: str) -> pd.DataFrame:
    if existing is None or existing.empty:
        out = new_wide
    else:
        out = existing.merge(new_wide, on=key_col, how="outer")
    # NaN -> 0 dla kolumn Run
    for c in out.columns:
        if c != key_col:
            out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0.0)
    return out


def analyze_wgs_from_manifest(
    sra_manifest_staged: pd.DataFrame,
    params: Dict[str, Any],
) -> Dict[str, pd.DataFrame]:
    cfg = params

    threads = int(cfg.get("threads", 8))
    kraken_db = str(cfg["kraken_db"])

    results_base = Path(str(cfg["results_dir"]))
    matrices_dir = Path(str(cfg["matrices_dir"]))
    _ensure_dir(results_base)
    _ensure_dir(matrices_dir)

    kraken_long_dir = matrices_dir / "kraken_taxa_long"
    motus_long_dir = matrices_dir / "motus_taxa_long"
    _ensure_dir(kraken_long_dir)
    _ensure_dir(motus_long_dir)

    # macierze wynikowe (Tax x Run)
    kraken_matrix_path = matrices_dir / "kraken_matrix.tsv"
    motus_matrix_path = matrices_dir / "motus_matrix.tsv"

    cleanup_staged_fastq = bool(cfg.get("cleanup_staged_fastq", False))

    # kasowanie results/<RUN>
    cleanup_results_dir = bool(cfg.get("cleanup_results_dir", False))
    cleanup_results_dir_force = bool(cfg.get("cleanup_results_dir_force", False))

    # hub logów
    log_hub_dir = Path(str(cfg.get("log_hub_dir", "data/08_reporting/sra_orchestrator/logs")))
    _ensure_dir(log_hub_dir)

    df = sra_manifest_staged.copy()
    df = df[(df["strategy"] == "WGS") & (df["status"] == "OK")].copy()
    if df.empty:
        raise ValueError("No WGS runs with status=OK in sra_manifest_staged.")

    run_rows: List[Dict[str, Any]] = []
    kraken_long_all: List[pd.DataFrame] = []
    motus_long_all: List[pd.DataFrame] = []

    for _, r in df.iterrows():
        run = str(r["Run"])
        fastq_dir = Path(str(r["staged_fastq_dir"]))
        out_dir = results_base / run
        _ensure_dir(out_dir)

        run_hub = log_hub_dir / run
        _ensure_dir(run_hub)

        kraken_long_path = kraken_long_dir / f"{run}.tsv"
        motus_long_path = motus_long_dir / f"{run}.tsv"

        # logi w results
        kraken_stdout = out_dir / f"{run}.kraken.stdout.log"
        kraken_stderr = out_dir / f"{run}.kraken.stderr.log"
        motus_stdout = out_dir / f"{run}.motus.stdout.log"
        motus_stderr = out_dir / f"{run}.motus.stderr.log"

        t0 = time.perf_counter()
        status = "OK"
        msg = ""
        analysis_sec = 0.0

        fastq_mode, r1, r2, se = _pick_fastq_mode_and_paths(r)

        extra: Dict[str, Any] = {
            "fastq_mode": fastq_mode,
            "r1": r1,
            "r2": r2,
            "se": se,
            "kraken_long_tsv": str(kraken_long_path),
            "motus_long_tsv": str(motus_long_path),
            "kraken_report": "",
            "motus_tsv": "",
        }

        try:
            t = time.perf_counter()

            # ---------- KRAKEN2 ----------
            kraken_report = out_dir / f"{run}.kreport"
            kraken_out = out_dir / f"{run}.kraken"

            if fastq_mode == "PE":
                cmd_k = (
                    f"mkdir -p '{out_dir}'; "
                    f"kraken2 --db '{kraken_db}' --threads {threads} --paired "
                    f"'{r1}' '{r2}' "
                    f"--report '{kraken_report}' "
                    f"--output '{kraken_out}'"
                )
            else:
                cmd_k = (
                    f"mkdir -p '{out_dir}'; "
                    f"kraken2 --db '{kraken_db}' --threads {threads} "
                    f"'{se}' "
                    f"--report '{kraken_report}' "
                    f"--output '{kraken_out}'"
                )

            resk = _run(cmd_k)
            _write_text(kraken_stdout, resk.stdout)
            _write_text(kraken_stderr, resk.stderr)
            if resk.returncode != 0:
                raise RuntimeError(f"Kraken2 failed: {resk.stderr or resk.stdout}")

            # ---------- mOTUs4 ----------
            motus_tsv = out_dir / f"{run}.motus.tsv"

            if fastq_mode == "PE":
                cmd_m = (
                    f"mkdir -p '{out_dir}'; "
                    f"mamba run -n motus4 motus profile "
                    f"-f '{r1}' -r '{r2}' "
                    f"-n '{run}' -t {threads} "
                    f"-o '{motus_tsv}' "
                    f"1> '{motus_stdout}' "
                    f"2> '{motus_stderr}'"
                )
            else:
                cmd_m = (
                    f"mkdir -p '{out_dir}'; "
                    f"mamba run -n motus4 motus profile "
                    f"-s '{se}' "
                    f"-n '{run}' -t {threads} "
                    f"-o '{motus_tsv}' "
                    f"1> '{motus_stdout}' "
                    f"2> '{motus_stderr}'"
                )

            resm = _run(cmd_m)
            analysis_sec = time.perf_counter() - t
            if resm.returncode != 0:
                raise RuntimeError(f"mOTUs failed: {resm.stderr or resm.stdout}")

            # parse outputs
            kraken_report_path = _best_kraken_report(out_dir)
            if kraken_report_path is None:
                raise RuntimeError("Kraken output not found (no .kreport in out_dir)")

            if not motus_tsv.exists():
                raise RuntimeError(f"mOTUs TSV not found: {motus_tsv}")

            # Kraken long -> zapisz (pełny long)
            kdf, ksum = parse_kraken_report(kraken_report_path, run)
            kdf.to_csv(kraken_long_path, sep="\t", index=False)
            extra.update(ksum)

            # MOTUs long -> zapisz
            mdf, msum = parse_motus_tsv(motus_tsv, run)
            mdf.to_csv(motus_long_path, sep="\t", index=False)
            extra.update(msum)

            # przygotuj long do pivot:
            # Kraken: użyj "name" jako Tax key + clade_reads jako count
            k_long = kdf[["Run", "name", "clade_reads"]].copy()
            k_long = k_long.rename(columns={"name": "Tax", "clade_reads": "count"})
            k_long["count"] = pd.to_numeric(k_long["count"], errors="coerce").fillna(0.0).astype(float)
            kraken_long_all.append(k_long)

            # mOTUs: użyj "Taxonomy" jako Tax key + count
            m_long = mdf[["Run", "Taxonomy", "count"]].copy()
            m_long = m_long.rename(columns={"Taxonomy": "Tax"})
            m_long["count"] = pd.to_numeric(m_long["count"], errors="coerce").fillna(0.0).astype(float)
            motus_long_all.append(m_long)

        except Exception as e:
            status = "FAILED"
            msg = str(e).replace("\n", " ")[:2000]

        total_sec = time.perf_counter() - t0

        # --- DUPLIKACJA logów do HUB (zostają też w results) ---
        _copy_if_exists(kraken_stdout, run_hub / kraken_stdout.name)
        _copy_if_exists(kraken_stderr, run_hub / kraken_stderr.name)
        _copy_if_exists(motus_stdout, run_hub / motus_stdout.name)
        _copy_if_exists(motus_stderr, run_hub / motus_stderr.name)

        # dodatkowo reporty
        _copy_if_exists(out_dir / f"{run}.kreport", run_hub / f"{run}.kreport")

        # cleanup staged fastq AFTER analysis (tylko jeśli OK)
        if cleanup_staged_fastq and status == "OK":
            shutil.rmtree(fastq_dir, ignore_errors=True)

        # cleanup results/<RUN> AFTER log duplication + AFTER parsing
        if cleanup_results_dir and (status == "OK" or cleanup_results_dir_force):
            shutil.rmtree(out_dir, ignore_errors=True)

        run_rows.append(
            dict(
                Run=run,
                strategy="WGS",
                status=status,
                message=msg,
                results_dir=str(out_dir),
                analysis_sec=round(analysis_sec, 3),
                total_sec=round(total_sec, 3),
                staged_fastq_dir=str(fastq_dir),
                log_hub_dir=str(run_hub),
                **extra,
            )
        )

    # -------------------------
    # RUN MATRIX
    # -------------------------
    run_matrix = pd.DataFrame(run_rows)

    # -------------------------
    # KRAKEN MATRIX (Tax x Run) incremental
    # -------------------------
    if kraken_long_all:
        k_all = pd.concat(kraken_long_all, ignore_index=True)
        k_wide = (
            k_all.pivot_table(index="Tax", columns="Run", values="count", aggfunc="sum", fill_value=0.0)
            .reset_index()
        )
        existing = _load_matrix(kraken_matrix_path, key_col="Tax")
        kraken_matrix = _merge_matrix(existing, k_wide, key_col="Tax")
        kraken_matrix.to_csv(kraken_matrix_path, sep="\t", index=False)
    else:
        kraken_matrix = _load_matrix(kraken_matrix_path, key_col="Tax")

    # -------------------------
    # MOTUS MATRIX (Tax x Run) incremental
    # -------------------------
    if motus_long_all:
        m_all = pd.concat(motus_long_all, ignore_index=True)
        m_wide = (
            m_all.pivot_table(index="Tax", columns="Run", values="count", aggfunc="sum", fill_value=0.0)
            .reset_index()
        )
        existing = _load_matrix(motus_matrix_path, key_col="Tax")
        motus_matrix = _merge_matrix(existing, m_wide, key_col="Tax")
        motus_matrix.to_csv(motus_matrix_path, sep="\t", index=False)
    else:
        motus_matrix = _load_matrix(motus_matrix_path, key_col="Tax")

    return {
        "sra_profile_wgs_run_matrix": run_matrix,
        "sra_profile_wgs_kraken_matrix": kraken_matrix,
        "sra_profile_wgs_motus_matrix": motus_matrix,
    }
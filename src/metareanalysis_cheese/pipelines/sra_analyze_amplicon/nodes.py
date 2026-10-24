from __future__ import annotations

import csv
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd


# ----------------------------
# helpers
# ----------------------------
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


def _assert_exists(p: Path, label: str) -> None:
    if not p.exists():
        raise FileNotFoundError(f"{label} does not exist: {p}")


def _has_pigz() -> bool:
    return _run("command -v pigz >/dev/null 2>&1").returncode == 0


def _gzip_copy(src: Path, dst: Path, threads: int) -> None:
    """Create gz copy of src at dst WITHOUT deleting src."""
    _ensure_dir(dst.parent)
    if _has_pigz():
        cmd = f"pigz -p {threads} -c '{src}' > '{dst}'"
    else:
        cmd = f"gzip -c '{src}' > '{dst}'"
    res = _run(cmd)
    if res.returncode != 0:
        raise RuntimeError(f"gzip copy failed: {cmd}\nstderr={res.stderr[:2000]}")
    _assert_exists(dst, f"gz copy of {src.name}")


def _detect_fastq_mode(r1: Optional[str], r2: Optional[str], se: Optional[str]) -> str:
    if r1 and r2:
        return "PE"
    if se:
        return "SE"
    raise ValueError("Cannot detect fastq mode: need (r1+r2) or (se).")


def _pick_fastq_mode_and_paths(row: pd.Series) -> Tuple[str, str, str, str]:
    """
    Try row.r1/r2/se first; if missing -> autodetect in staged_fastq_dir.
    Returns: (mode, r1, r2, se)
    """
    r1 = str(row.get("r1") or "").strip()
    r2 = str(row.get("r2") or "").strip()
    se = str(row.get("se") or "").strip()

    if r1 and r2:
        return "PE", r1, r2, ""
    if se:
        return "SE", "", "", se

    fastq_dir = Path(str(row.get("staged_fastq_dir") or "")).resolve()
    run = str(row.get("Run") or "").strip()

    if not run:
        raise RuntimeError("Missing Run in manifest row (cannot autodetect fastq filenames).")
    if not fastq_dir.exists():
        raise RuntimeError(f"staged_fastq_dir does not exist for run={run}: {fastq_dir}")

    # paired candidates: SRRxxxx_1.fastq(.gz) etc.
    r1_c = sorted(fastq_dir.glob(f"{run}*_1.fastq")) + sorted(fastq_dir.glob(f"{run}*_1.fastq.gz"))
    r2_c = sorted(fastq_dir.glob(f"{run}*_2.fastq")) + sorted(fastq_dir.glob(f"{run}*_2.fastq.gz"))
    if r1_c and r2_c:
        return "PE", str(r1_c[0]), str(r2_c[0]), ""

    # single-end: SRRxxxx*.fastq(.gz) excluding _1/_2
    se_c = sorted(fastq_dir.glob(f"{run}*.fastq")) + sorted(fastq_dir.glob(f"{run}*.fastq.gz"))
    se_c = [p for p in se_c if not re.search(r"_[12]\.fastq(\.gz)?$", p.name)]
    if se_c:
        return "SE", "", "", str(se_c[0])

    # last resort: if directory has exactly 2 fastqs -> assume paired, 1 -> single
    all_fq = sorted(fastq_dir.glob("*.fastq")) + sorted(fastq_dir.glob("*.fastq.gz"))
    if len(all_fq) == 2:
        return "PE", str(all_fq[0]), str(all_fq[1]), ""
    if len(all_fq) == 1:
        return "SE", "", "", str(all_fq[0])

    raise RuntimeError(f"Cannot detect FASTQ inputs for run={run} in staged_fastq_dir={fastq_dir}")


def _write_qiime2_manifest_pe_tsv(manifest_path: Path, run: str, r1_gz: Path, r2_gz: Path) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    if manifest_path.exists():
        manifest_path.unlink()
    with manifest_path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"])
        w.writerow([run, str(r1_gz.resolve()), str(r2_gz.resolve())])


def _write_qiime2_manifest_se_tsv(manifest_path: Path, run: str, se_gz: Path) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    if manifest_path.exists():
        manifest_path.unlink()
    with manifest_path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample-id", "absolute-filepath"])
        w.writerow([run, str(se_gz.resolve())])


def _read_biom_tsv_single_sample(feature_table_tsv: Path) -> pd.DataFrame:
    """
    biom convert --to-tsv produces:
      # Constructed from biom file
      #OTU ID <tab> SAMPLE
      FeatureOrTax <tab> count
    """
    rows: List[List[str]] = []
    with feature_table_tsv.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(line.split("\t"))

    if not rows:
        raise RuntimeError(f"Empty feature table TSV after removing comments: {feature_table_tsv}")

    header = rows[0]
    data = rows[1:]
    if len(header) < 2:
        raise RuntimeError(f"Unexpected feature table header: {header}")

    sample_id = header[1]
    out: List[Tuple[str, str]] = []
    for r in data:
        if len(r) < 2:
            continue
        out.append((r[0], r[1]))

    df = pd.DataFrame(out, columns=["Tax", sample_id])
    df[sample_id] = pd.to_numeric(df[sample_id], errors="coerce").fillna(0).astype(float)
    return df


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
    for c in out.columns:
        if c != key_col:
            out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0.0)
    return out


def _matrix_has_run_col(matrix_path: Path, run: str) -> bool:
    """Fast check if run column already present (read only header)."""
    if not matrix_path.exists():
        return False
    try:
        with matrix_path.open("r", encoding="utf-8", errors="replace") as f:
            header = f.readline().rstrip("\n")
        if not header:
            return False
        cols = header.split("\t")
        return run in cols[1:]  # after Tax
    except Exception:
        return False


# ----------------------------
# MAIN
# ----------------------------
def analyze_amplicon_from_manifest(
    sra_manifest_staged: pd.DataFrame,
    params: Dict[str, Any],
) -> Dict[str, pd.DataFrame]:
    """
    WGS-style AMPLICON:
      - per-run long TSV: matrices/amplicon_taxa_long/<RUN>.tsv
      - incremental wide matrix: matrices/tax_matrix.tsv (Tax x Run)
      - run_matrix with status + metrics
      - logs in results/<RUN> + copied to log_hub_dir/<RUN>/

    Params supported:
      - target_run: "" or "SRR..."
      - skip_if_done: bool (skip if column exists in tax_matrix.tsv)
      - log_hub_dir: path
      - cleanup_staged_fastq: bool
      - cleanup_results_dir: bool
      - cleanup_results_dir_force: bool
      - add_total_reads_row: bool
      - add_unassigned_row: bool
    """
    cfg = params

    threads = int(cfg.get("threads", 8))
    qiime2_env = str(cfg.get("qiime2_env", "qiime2-amplicon-2026.1"))

    results_base = Path(str(cfg["results_dir"]))
    matrices_dir = Path(str(cfg["matrices_dir"]))
    _ensure_dir(results_base)
    _ensure_dir(matrices_dir)

    taxa_long_dir = matrices_dir / "amplicon_taxa_long"
    _ensure_dir(taxa_long_dir)

    tax_matrix_path = matrices_dir / "tax_matrix.tsv"

    target_run = str(cfg.get("target_run", "") or "").strip()
    skip_if_done = bool(cfg.get("skip_if_done", False))

    cleanup_staged_fastq = bool(cfg.get("cleanup_staged_fastq", False))
    cleanup_results_dir = bool(cfg.get("cleanup_results_dir", False))
    cleanup_results_dir_force = bool(cfg.get("cleanup_results_dir_force", False))

    add_total_reads_row = bool(cfg.get("add_total_reads_row", False))
    add_unassigned_row = bool(cfg.get("add_unassigned_row", False))

    log_hub_dir = Path(str(cfg.get("log_hub_dir", "data/08_reporting/sra_orchestrator/logs")))
    _ensure_dir(log_hub_dir)

    silva_classifier_qza = Path(str(cfg["silva_classifier_qza"]))
    _assert_exists(silva_classifier_qza, "silva_classifier_qza")

    trunc_len_f = int(cfg.get("trunc_len_f", 0))
    trunc_len_r = int(cfg.get("trunc_len_r", 0))
    trunc_len_s = int(cfg.get("trunc_len_s", 0))

    tax_level = str(cfg.get("tax_level", "genus")).lower()
    level_map = {
        "domain": 1,
        "phylum": 2,
        "class": 3,
        "order": 4,
        "family": 5,
        "genus": 6,
        "species": 7,
    }
    collapse_level = int(level_map.get(tax_level, 6))

    df = sra_manifest_staged.copy()
    df = df[(df["strategy"] == "AMPLICON") & (df["status"] == "OK")].copy()
    if target_run:
        df = df[df["Run"].astype(str) == target_run].copy()

    if df.empty:
        raise ValueError("No AMPLICON runs with status=OK in sra_manifest_staged (after filtering).")

    run_rows: List[Dict[str, Any]] = []
    long_all: List[pd.DataFrame] = []

    totals_map: Dict[str, float] = {}
    unassigned_map: Dict[str, float] = {}

    for _, row in df.iterrows():
        run = str(row["Run"])

        out_dir = results_base / run
        _ensure_dir(out_dir)

        run_hub = log_hub_dir / run
        _ensure_dir(run_hub)

        staged_dir = Path(str(row.get("staged_fastq_dir") or ""))

        taxa_long_path = taxa_long_dir / f"{run}.tsv"

        stdout_log = out_dir / f"{run}.qiime2.stdout.log"
        stderr_log = out_dir / f"{run}.qiime2.stderr.log"

        tmp_fastq_dir = out_dir / "tmp_fastq"

        # skip if already present in matrix header
        if skip_if_done and _matrix_has_run_col(tax_matrix_path, run):
            run_rows.append(
                dict(
                    Run=run,
                    strategy="AMPLICON",
                    status="SKIPPED",
                    message="skip_if_done: already present in tax_matrix",
                    results_dir=str(out_dir),
                    analysis_sec=0.0,
                    total_sec=0.0,
                    staged_fastq_dir=str(staged_dir),
                    log_hub_dir=str(run_hub),
                    fastq_mode="",
                    r1=str(row.get("r1") or ""),
                    r2=str(row.get("r2") or ""),
                    se=str(row.get("se") or ""),
                    qiime2_env=qiime2_env,
                    silva_classifier_qza=str(silva_classifier_qza),
                    tax_level=tax_level,
                    collapse_level=collapse_level,
                    taxa_long_tsv=str(taxa_long_path),
                )
            )
            continue

        # qiime2 paths
        manifest_tsv = out_dir / f"{run}.manifest.tsv"
        demux_qza = out_dir / f"{run}.demux.qza"

        dada2_dir = out_dir / "dada2"
        rep_seqs_qza = dada2_dir / "representative_sequences.qza"
        table_qza = dada2_dir / "table.qza"

        taxonomy_qza = out_dir / f"{run}.taxonomy.qza"

        collapsed_qza = out_dir / f"{run}.table_L{collapse_level}.qza"
        export_dir = out_dir / f"export_L{collapse_level}"
        feature_table_biom = export_dir / "feature-table.biom"
        feature_table_tsv = out_dir / f"{run}.table_L{collapse_level}.tsv"

        t0 = time.perf_counter()
        status = "OK"
        msg = ""
        analysis_sec = 0.0

        total_reads_retained: Optional[float] = None
        unassigned_reads: Optional[float] = None

        # IMPORTANT: always pick inputs robustly (fallback to staged_fastq_dir)
        fastq_mode = ""
        r1 = ""
        r2 = ""
        se = ""

        try:
            fastq_mode, r1, r2, se = _pick_fastq_mode_and_paths(row)

            t = time.perf_counter()

            # We accept gz already (from staged_fastq_dir); but we still create gz copies into results/<RUN>/tmp_fastq
            # to ensure QIIME2 manifest points to stable paths under out_dir.
            if fastq_mode == "PE":
                r1p = Path(r1)
                r2p = Path(r2)
                _assert_exists(r1p, f"{run} R1")
                _assert_exists(r2p, f"{run} R2")

                r1gz = tmp_fastq_dir / f"{run}_1.fastq.gz"
                r2gz = tmp_fastq_dir / f"{run}_2.fastq.gz"

                # if inputs already .gz, just copy; else gzip copy
                _ensure_dir(tmp_fastq_dir)
                if r1p.suffixes[-2:] == [".fastq", ".gz"] or r1p.name.endswith(".fastq.gz"):
                    shutil.copy2(r1p, r1gz)
                else:
                    _gzip_copy(r1p, r1gz, threads=threads)

                if r2p.suffixes[-2:] == [".fastq", ".gz"] or r2p.name.endswith(".fastq.gz"):
                    shutil.copy2(r2p, r2gz)
                else:
                    _gzip_copy(r2p, r2gz, threads=threads)

                _write_qiime2_manifest_pe_tsv(manifest_tsv, run, r1gz, r2gz)

                import_cmd = (
                    "qiime tools import "
                    "--type 'SampleData[PairedEndSequencesWithQuality]' "
                    f"--input-path '{manifest_tsv}' "
                    f"--output-path '{demux_qza}' "
                    "--input-format PairedEndFastqManifestPhred33V2"
                )
                denoise_cmd = (
                    "qiime dada2 denoise-paired "
                    f"--i-demultiplexed-seqs '{demux_qza}' "
                    f"--p-trunc-len-f {trunc_len_f} "
                    f"--p-trunc-len-r {trunc_len_r} "
                    f"--p-n-threads {threads} "
                    f"--output-dir '{dada2_dir}'"
                )
            else:
                sep = Path(se)
                _assert_exists(sep, f"{run} SE")

                segz = tmp_fastq_dir / f"{run}.fastq.gz"
                _ensure_dir(tmp_fastq_dir)
                if sep.suffixes[-2:] == [".fastq", ".gz"] or sep.name.endswith(".fastq.gz"):
                    shutil.copy2(sep, segz)
                else:
                    _gzip_copy(sep, segz, threads=threads)

                _write_qiime2_manifest_se_tsv(manifest_tsv, run, segz)

                import_cmd = (
                    "qiime tools import "
                    "--type 'SampleData[SequencesWithQuality]' "
                    f"--input-path '{manifest_tsv}' "
                    f"--output-path '{demux_qza}' "
                    "--input-format SingleEndFastqManifestPhred33V2"
                )
                denoise_cmd = (
                    "qiime dada2 denoise-single "
                    f"--i-demultiplexed-seqs '{demux_qza}' "
                    f"--p-trunc-len {trunc_len_s} "
                    f"--p-n-threads {threads} "
                    f"--output-dir '{dada2_dir}'"
                )

            classify_cmd = (
                "qiime feature-classifier classify-sklearn "
                f"--i-classifier '{silva_classifier_qza}' "
                f"--i-reads '{rep_seqs_qza}' "
                f"--p-n-jobs {threads} "
                f"--o-classification '{taxonomy_qza}'"
            )

            collapse_cmd = (
                "qiime taxa collapse "
                f"--i-table '{table_qza}' "
                f"--i-taxonomy '{taxonomy_qza}' "
                f"--p-level {collapse_level} "
                f"--o-collapsed-table '{collapsed_qza}'"
            )

            export_cmd = (
                f"rm -rf '{export_dir}' && mkdir -p '{export_dir}'; "
                f"qiime tools export --input-path '{collapsed_qza}' --output-path '{export_dir}'; "
                f"biom convert -i '{feature_table_biom}' -o '{feature_table_tsv}' --to-tsv"
            )

            bash_script = (
                "set -euo pipefail; "
                f"mamba run -n {qiime2_env} bash -lc "
                + "'"
                + f"mkdir -p '{out_dir}'; "
                + import_cmd
                + "; "
                + denoise_cmd
                + "; "
                + classify_cmd
                + "; "
                + collapse_cmd
                + "; "
                + export_cmd
                + "'"
            )

            res = _run(f"{bash_script} 1> '{stdout_log}' 2> '{stderr_log}'")
            analysis_sec = time.perf_counter() - t

            if res.returncode != 0:
                tail = ""
                if stderr_log.exists():
                    tail = "\n".join(stderr_log.read_text(errors="replace").splitlines()[-200:])
                raise RuntimeError(f"QIIME2 failed. stderr_tail={tail[:4000]}")

            _assert_exists(feature_table_tsv, "collapsed feature-table TSV")

            ft = _read_biom_tsv_single_sample(feature_table_tsv)
            sample_cols = [c for c in ft.columns if c != "Tax"]
            if not sample_cols:
                raise RuntimeError(f"No sample column found in {feature_table_tsv}")
            sample_col = sample_cols[0]

            ft["Tax"] = ft["Tax"].fillna("Unassigned").astype(str)

            total_reads_retained = float(ft[sample_col].sum())
            unassigned_reads = float(ft.loc[ft["Tax"] == "Unassigned", sample_col].sum())

            long_df = ft.rename(columns={sample_col: "count"}).copy()
            long_df.insert(0, "Run", run)
            long_df = long_df[["Run", "Tax", "count"]].copy()
            long_df["count"] = pd.to_numeric(long_df["count"], errors="coerce").fillna(0.0).astype(float)

            _ensure_dir(taxa_long_path.parent)
            long_df.to_csv(taxa_long_path, sep="\t", index=False)
            _assert_exists(taxa_long_path, "amplicon taxa_long TSV")

            long_all.append(long_df)
            totals_map[run] = float(total_reads_retained)
            unassigned_map[run] = float(unassigned_reads)

        except Exception as e:
            status = "FAILED"
            msg = str(e).replace("\n", " ")[:4000]

        total_sec = time.perf_counter() - t0

        # copy logs to hub (even if failed)
        _copy_if_exists(stdout_log, run_hub / stdout_log.name)
        _copy_if_exists(stderr_log, run_hub / stderr_log.name)

        # cleanup tmp gz only on OK
        if status == "OK" and tmp_fastq_dir.exists():
            shutil.rmtree(tmp_fastq_dir, ignore_errors=True)

        # optional cleanup staged fastq after OK
        if cleanup_staged_fastq and status == "OK" and staged_dir.exists():
            shutil.rmtree(staged_dir, ignore_errors=True)

        # optional cleanup results/<RUN>
        if cleanup_results_dir and (status == "OK" or cleanup_results_dir_force):
            shutil.rmtree(out_dir, ignore_errors=True)

        run_rows.append(
            dict(
                Run=run,
                strategy="AMPLICON",
                status=status,
                message=msg,
                results_dir=str(out_dir),
                analysis_sec=round(analysis_sec, 3),
                total_sec=round(total_sec, 3),
                staged_fastq_dir=str(staged_dir),
                log_hub_dir=str(run_hub),
                fastq_mode=fastq_mode,
                r1=r1,
                r2=r2,
                se=se,
                qiime2_env=qiime2_env,
                silva_classifier_qza=str(silva_classifier_qza),
                tax_level=tax_level,
                collapse_level=collapse_level,
                trunc_len_f=trunc_len_f,
                trunc_len_r=trunc_len_r,
                trunc_len_s=trunc_len_s,
                taxa_long_tsv=str(taxa_long_path),
                total_reads_retained=total_reads_retained if status == "OK" else None,
                unassigned_reads=unassigned_reads if status == "OK" else None,
            )
        )

    run_matrix = pd.DataFrame(run_rows)

    # -------------------------
    # TAX MATRIX incremental
    # -------------------------
    if long_all:
        all_long = pd.concat(long_all, ignore_index=True)
        new_wide = (
            all_long.pivot_table(index="Tax", columns="Run", values="count", aggfunc="sum", fill_value=0.0)
            .reset_index()
        )

        existing = _load_matrix(tax_matrix_path, key_col="Tax")
        tax_matrix = _merge_matrix(existing, new_wide, key_col="Tax")

        if add_total_reads_row:
            row = {"Tax": "__TOTAL_READS_RETAINED__"}
            for rr, val in totals_map.items():
                row[rr] = float(val)
            tax_matrix = _merge_matrix(tax_matrix, pd.DataFrame([row]), key_col="Tax")

        if add_unassigned_row:
            row = {"Tax": "__UNASSIGNED__"}
            for rr, val in unassigned_map.items():
                row[rr] = float(val)
            tax_matrix = _merge_matrix(tax_matrix, pd.DataFrame([row]), key_col="Tax")

        tax_matrix.to_csv(tax_matrix_path, sep="\t", index=False)
    else:
        tax_matrix = _load_matrix(tax_matrix_path, key_col="Tax")

    return {
        "sra_profile_amplicon_run_matrix": run_matrix,
        "sra_profile_amplicon_tax_matrix": tax_matrix,
    }
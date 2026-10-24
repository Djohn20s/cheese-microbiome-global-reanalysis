# src/metareanalysis_cheese/pipelines/sra_orchestrator/nodes.py
from __future__ import annotations

import csv
import faulthandler
import shutil
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

import pandas as pd

from metareanalysis_cheese.pipelines.sra_fetch_unpack.nodes import build_runs_table, fetch_unpack_stage
from metareanalysis_cheese.pipelines.sra_analyze_wgs.nodes import analyze_wgs_from_manifest
from metareanalysis_cheese.pipelines.sra_analyze_amplicon.nodes import analyze_amplicon_from_manifest


# ----------------------------
# small helpers
# ----------------------------
STATE_BASE_COLS = ["Run", "strategy", "status", "step", "message", "updated_at"]


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _now_iso() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _as_list(x: Any) -> List[str]:
    if x is None:
        return []
    if isinstance(x, list):
        return [str(v).strip() for v in x if str(v).strip()]
    s = str(x).strip()
    if not s:
        return []
    return [v.strip() for v in s.split(",") if v.strip()]


def _normalize_strategy(s: str) -> str:
    s = (s or "").strip().upper()
    if s in {"WGS", "AMPLICON"}:
        return s
    if s in {"AMP", "AMPLICONS", "16S"}:
        return "AMPLICON"
    if s in {"BOTH", "ALL", "*", ""}:
        return "BOTH"
    return s


def _allowed_strategies(cfg: Dict[str, Any]) -> Set[str]:
    """
    New options (choose one):
      - only_strategy: "WGS" | "AMPLICON" | "BOTH"  (recommended)
      - run_wgs: bool, run_amplicon: bool (legacy-style)
      - strategies: ["WGS","AMPLICON"] (list)
    Precedence: only_strategy > strategies > (run_wgs/run_amplicon) > default BOTH
    """
    only_strategy = _normalize_strategy(str(cfg.get("only_strategy", "")).strip())
    if only_strategy in {"WGS", "AMPLICON"}:
        return {only_strategy}
    if only_strategy == "BOTH":
        return {"WGS", "AMPLICON"}

    strategies = cfg.get("strategies", None)
    if strategies is not None:
        st = {_normalize_strategy(x) for x in _as_list(strategies) if _normalize_strategy(x) in {"WGS", "AMPLICON"}}
        if st:
            return st

    run_wgs = bool(cfg.get("run_wgs", True))
    run_amp = bool(cfg.get("run_amplicon", True))
    out: Set[str] = set()
    if run_wgs:
        out.add("WGS")
    if run_amp:
        out.add("AMPLICON")
    return out or {"WGS", "AMPLICON"}


def _read_header_columns(path: Path) -> List[str]:
    if not path.exists():
        return []
    try:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            header = f.readline().rstrip("\n")
        return header.split("\t") if header else []
    except Exception:
        return []


def _done_runs_from_matrices(wgs_matrices_dir: Path, amp_matrices_dir: Path) -> Tuple[Set[str], Set[str]]:
    """
    Reads Run columns already present in matrix TSVs.
    WGS: kraken_matrix.tsv -> cols[1:]
    AMP: tax_matrix.tsv -> cols[1:]
    """
    done_wgs: Set[str] = set()
    done_amp: Set[str] = set()

    wgs_path = wgs_matrices_dir / "kraken_matrix.tsv"
    cols = _read_header_columns(wgs_path)
    if cols:
        done_wgs.update([str(c).strip() for c in cols[1:] if str(c).strip()])

    amp_path = amp_matrices_dir / "tax_matrix.tsv"
    cols = _read_header_columns(amp_path)
    if cols:
        done_amp.update([str(c).strip() for c in cols[1:] if str(c).strip()])

    return done_wgs, done_amp


def _safe_copy_glob(src_dir: Path, patterns: List[str], dst_dir: Path) -> None:
    if not src_dir.exists():
        return
    _ensure_dir(dst_dir)
    for pat in patterns:
        for p in src_dir.glob(pat):
            try:
                if p.is_file():
                    shutil.copy2(p, dst_dir / p.name)
            except Exception:
                pass


# ----------------------------
# STATE (NO pandas): Run -> row dict, written as TSV
# ----------------------------
def _load_state_map(state_path: Path) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    """
    Returns:
      - state_map: {Run: {col: value_as_str}}
      - fieldnames: columns in file (if exists) else STATE_BASE_COLS
    """
    if not state_path.exists():
        return {}, STATE_BASE_COLS.copy()

    try:
        with state_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            fieldnames = list(reader.fieldnames or [])
            if not fieldnames:
                return {}, STATE_BASE_COLS.copy()

            for c in STATE_BASE_COLS:
                if c not in fieldnames:
                    fieldnames.append(c)

            state_map: Dict[str, Dict[str, str]] = {}
            for row in reader:
                run = str(row.get("Run", "") or "").strip()
                if not run:
                    continue
                norm = {k: "" if row.get(k) is None else str(row.get(k)) for k in fieldnames}
                norm["Run"] = run
                state_map[run] = norm
            return state_map, fieldnames
    except Exception:
        return {}, STATE_BASE_COLS.copy()


def _write_state_map(state_path: Path, state_map: Dict[str, Dict[str, str]], fieldnames: List[str]) -> None:
    _ensure_dir(state_path.parent)

    base = [c for c in STATE_BASE_COLS if c in fieldnames]
    extra = [c for c in fieldnames if c not in base]
    fns = base + extra

    runs_sorted = sorted(state_map.keys())

    tmp = state_path.with_suffix(".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fns, delimiter="\t")
        w.writeheader()
        for run in runs_sorted:
            row = state_map[run]
            out = {k: "" for k in fns}
            for k, v in row.items():
                if k not in out:
                    continue
                out[k] = "" if v is None else str(v)
            out["Run"] = run
            w.writerow(out)
    tmp.replace(state_path)


def _upsert_state(
    state_path: Path,
    state_map: Dict[str, Dict[str, str]],
    fieldnames: List[str],
    row: Dict[str, Any],
) -> None:
    run = str(row.get("Run", "") or "").strip()
    if not run:
        return

    if run not in state_map:
        state_map[run] = {c: "" for c in fieldnames}
        state_map[run]["Run"] = run

    for k in row.keys():
        if k not in fieldnames:
            fieldnames.append(k)
            for r in state_map.values():
                r[k] = ""

    for k, v in row.items():
        state_map[run][k] = "" if v is None else str(v)

    state_map[run]["updated_at"] = _now_iso()
    _write_state_map(state_path, state_map, fieldnames)


def _should_skip_by_state(run: str, state_map: Dict[str, Dict[str, str]], continue_mode: bool, rerun_failed: bool) -> bool:
    if not continue_mode:
        return False
    if run not in state_map:
        return False
    status = str(state_map[run].get("status", "")).upper()
    if status == "OK":
        return True
    if status == "FAILED" and not rerun_failed:
        return True
    return False


# ----------------------------
# MAIN
# ----------------------------
def run_sra_orchestrator(
    final_wgs: pd.DataFrame,
    final_amplicon: pd.DataFrame,
    params_orchestrator: Dict[str, Any],
    params_fetch: Dict[str, Any],
    params_wgs: Dict[str, Any],
    params_amplicon: Dict[str, Any],
) -> pd.DataFrame:
    faulthandler.enable(all_threads=True, file=sys.stderr)

    cfg = params_orchestrator

    base_dir = Path(str(cfg.get("orchestrator_dir", "data/08_reporting/sra_orchestrator")))
    state_path = base_dir / "state.tsv"
    logs_dir = base_dir / "logs"
    _ensure_dir(base_dir)
    _ensure_dir(logs_dir)

    continue_mode = bool(cfg.get("continue_mode", True))
    rerun_failed = bool(cfg.get("rerun_failed", False))
    stop_on_fail = bool(cfg.get("stop_on_fail", False))

    start_after = str(cfg.get("start_after_run", "")).strip()
    skip_runs = set(_as_list(cfg.get("skip_runs", [])))

    max_total = int(cfg.get("max_total", 0))
    max_wgs = int(cfg.get("max_wgs", 0))
    max_amp = int(cfg.get("max_amplicon", 0))

    skip_if_done_in_matrices = bool(cfg.get("skip_if_done_in_matrices", True))

    allowed = _allowed_strategies(cfg)  # <- NEW: choose only WGS / only AMPLICON / both

    # matrices "done"
    wgs_matrices_dir = Path(str(params_wgs["matrices_dir"]))
    amp_matrices_dir = Path(str(params_amplicon["matrices_dir"]))
    done_wgs, done_amp = _done_runs_from_matrices(wgs_matrices_dir, amp_matrices_dir)

    # state map (NO pandas)
    state_map, fieldnames = _load_state_map(state_path)

    # run list
    runs_table = build_runs_table(final_wgs=final_wgs, final_amplicon=final_amplicon)
    runs_table = runs_table.dropna(subset=["Run"]).drop_duplicates(subset=["Run"]).reset_index(drop=True)

    # filter by allowed strategies (NEW)
    if "strategy" in runs_table.columns and allowed:
        runs_table = runs_table[runs_table["strategy"].astype(str).str.upper().isin(allowed)].copy()
        runs_table = runs_table.reset_index(drop=True)

    started = (start_after == "")
    total_done = 0
    wgs_done = 0
    amp_done = 0

    for row in runs_table.itertuples(index=False):
        run = str(getattr(row, "Run"))
        strategy = str(getattr(row, "strategy")).upper().strip()

        # extra safety: if build_runs_table gave something unexpected
        if allowed and strategy not in allowed:
            continue

        if not started:
            if run == start_after:
                started = True
            continue
        if run == start_after:
            continue

        if skip_runs and run in skip_runs:
            _upsert_state(
                state_path,
                state_map,
                fieldnames,
                dict(Run=run, strategy=strategy, status="SKIPPED", step="USER_SKIP", message="skip_runs"),
            )
            continue

        if max_total > 0 and total_done >= max_total:
            break
        if strategy == "WGS" and max_wgs > 0 and wgs_done >= max_wgs:
            continue
        if strategy == "AMPLICON" and max_amp > 0 and amp_done >= max_amp:
            continue

        if _should_skip_by_state(run, state_map, continue_mode=continue_mode, rerun_failed=rerun_failed):
            continue

        if skip_if_done_in_matrices:
            if strategy == "WGS" and run in done_wgs:
                _upsert_state(
                    state_path,
                    state_map,
                    fieldnames,
                    dict(Run=run, strategy=strategy, status="SKIPPED", step="DONE_IN_MATRIX", message="already in kraken_matrix"),
                )
                continue
            if strategy == "AMPLICON" and run in done_amp:
                _upsert_state(
                    state_path,
                    state_map,
                    fieldnames,
                    dict(Run=run, strategy=strategy, status="SKIPPED", step="DONE_IN_MATRIX", message="already in tax_matrix"),
                )
                continue

        run_log_dir = logs_dir / run
        _ensure_dir(run_log_dir)

        _upsert_state(state_path, state_map, fieldnames, dict(Run=run, strategy=strategy, status="RUNNING", step="START", message=""))

        # 1) FETCH
        try:
            _upsert_state(state_path, state_map, fieldnames, dict(Run=run, strategy=strategy, status="RUNNING", step="FETCH", message=""))

            one = pd.DataFrame([{"Run": run, "strategy": strategy}])
            manifest = fetch_unpack_stage(runs_table=one, params=params_fetch)

            if manifest.empty:
                raise RuntimeError("fetch_unpack_stage returned empty manifest")
            m0 = manifest.iloc[0]
            if str(m0.get("status", "")).upper() != "OK":
                raise RuntimeError(f"fetch failed: {m0.get('message','')}")

            _upsert_state(state_path, state_map, fieldnames, dict(Run=run, strategy=strategy, status="RUNNING", step="FETCH_OK", message=""))

        except Exception as e:
            _upsert_state(state_path, state_map, fieldnames, dict(Run=run, strategy=strategy, status="FAILED", step="FETCH", message=str(e)[:2000]))
            if stop_on_fail:
                break
            continue

        # 2) ANALYZE
        try:
            _upsert_state(state_path, state_map, fieldnames, dict(Run=run, strategy=strategy, status="RUNNING", step="ANALYZE", message=""))

            if strategy == "WGS":
                # WGS node already does incremental matrices (kraken/motus) itself
                analyze_wgs_from_manifest(sra_manifest_staged=manifest, params=params_wgs)

                done_wgs.add(run)
                wgs_done += 1

                out_dir = Path(str(params_wgs["results_dir"])) / run
                _safe_copy_glob(out_dir, ["*.log", "*.kreport", "*.kraken", "*.tsv"], run_log_dir)

            elif strategy == "AMPLICON":
                # AMPLICON node (updated) does incremental tax_matrix itself
                analyze_amplicon_from_manifest(sra_manifest_staged=manifest, params=params_amplicon)

                done_amp.add(run)
                amp_done += 1

                out_dir = Path(str(params_amplicon["results_dir"])) / run
                _safe_copy_glob(out_dir, ["*.qiime2.*.log", "*.log", "*.tsv"], run_log_dir)

            else:
                raise ValueError(f"Unknown strategy={strategy}")

            _upsert_state(state_path, state_map, fieldnames, dict(Run=run, strategy=strategy, status="OK", step="DONE", message=""))
            total_done += 1

        except Exception as e:
            _upsert_state(state_path, state_map, fieldnames, dict(Run=run, strategy=strategy, status="FAILED", step="ANALYZE", message=str(e)[:2000]))
            if stop_on_fail:
                break
            continue

    return pd.DataFrame(list(state_map.values()))
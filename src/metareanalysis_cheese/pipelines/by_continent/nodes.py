from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple, List

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle
import geodatasets

# Centroidy kontynentów (możesz przesunąć, jeśli koła nachodzą).
CONTINENT_CENTROIDS: Dict[str, Tuple[float, float]] = {
    "North America": (-100, 45),
    "South America": (-60, -15),
    "Europe": (15, 52),
    "Africa": (20, 5),
    "Asia": (90, 40),
    "Oceania": (135, -25),
}

@dataclass(frozen=True)
class ContinentPieConfig:
    strategies: Tuple[str, str] = ("AMPLICON", "WGS")
    radius_scale: float = 0.85     # skala promienia
    radius_min: float = 2.0        # min promień (stopnie)
    radius_max: float = 10.0       # max promień
    min_total_to_plot: int = 1
    output_dpi: int = 300


def _radius_from_total(total: float, cfg: ContinentPieConfig) -> float:
    r = cfg.radius_scale * np.sqrt(max(total, 0.0))
    return float(np.clip(r, cfg.radius_min, cfg.radius_max))


def _draw_pie(ax, x: float, y: float, values: np.ndarray, colors: List, r: float):
    total = float(values.sum())
    if total <= 0:
        return

    angles = np.cumsum(values) / total * 360.0
    start = 0.0
    for v, end, c in zip(values, angles, colors):
        if v <= 0:
            continue
        ax.add_patch(Wedge((x, y), r, start, float(end), facecolor=c, edgecolor="white", linewidth=0.7))
        start = float(end)

    ax.add_patch(Circle((x, y), r, facecolor="none", edgecolor="white", linewidth=0.9))


def plot_by_continent_strategy_pies(
    by_continent: pd.DataFrame,
    cfg: ContinentPieConfig = ContinentPieConfig(),
) -> plt.Figure:
    """
    Wejście (long):
      Continent | strategy | n_bioprojects | pct_within_continent
    Rysuje mapę świata i koła na kontynentach: promień ~ suma n_bioprojects, wycinki = strategie.
    """

    # Walidacja kolumn (Twoje nazwy)
    required = {"Continent", "strategy", "n_bioprojects", "pct_within_continent"}
    missing = required - set(by_continent.columns)
    if missing:
        raise ValueError(f"by_continent missing columns: {sorted(missing)}")

    df = by_continent.copy()
    df["n_bioprojects"] = pd.to_numeric(df["n_bioprojects"], errors="coerce").fillna(0)

    # Pivot: Continent x strategy -> n_bioprojects
    piv = (
        df.pivot_table(index="Continent", columns="strategy", values="n_bioprojects", aggfunc="sum")
        .fillna(0)
    )

    # Dociągnij brakujące strategie jako kolumny 0
    for s in cfg.strategies:
        if s not in piv.columns:
            piv[s] = 0
    piv = piv[list(cfg.strategies)]

    # Mapa świata (Natural Earth - wbudowane w GeoPandas)
    world = gpd.read_file(geodatasets.get_path("naturalearth.land"))

    fig = plt.figure(figsize=(12, 6))
    ax = plt.gca()
    world.plot(ax=ax, linewidth=0.5, edgecolor="0.8", color="0.95")

    # Kolory (domyślna paleta matplotlib – nie ustawiamy ręcznie)
    cmap = plt.get_cmap("tab10")
    colors = [cmap(0), cmap(1)]  # AMPLICON, WGS

    # Koła
    for cont in piv.index:
        if cont not in CONTINENT_CENTROIDS:
            continue
        total = float(piv.loc[cont].sum())
        if total < cfg.min_total_to_plot:
            continue

        x, y = CONTINENT_CENTROIDS[cont]
        r = _radius_from_total(total, cfg)
        vals = piv.loc[cont].to_numpy(dtype=float)

        _draw_pie(ax, x, y, vals, colors, r)
        ax.text(x, y - r - 2.0, f"{cont}\n(n={int(round(total))})", ha="center", va="top", fontsize=9)

    ax.set_axis_off()
    ax.set_title("Bioprojects by continent: size = total count, pie = AMPLICON vs WGS")

    # Legenda
    handles = [
        Circle((0, 0), radius=0.2, facecolor=colors[0], edgecolor="none"),
        Circle((0, 0), radius=0.2, facecolor=colors[1], edgecolor="none"),
    ]
    ax.legend(handles, list(cfg.strategies), loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False)

    # Legenda rozmiaru (promień)
    totals = piv.sum(axis=1).values
    totals = totals[totals > 0]
    size_vals = [10, 25, 50] if len(totals) == 0 else [int(round(v)) for v in np.quantile(totals, [0.25, 0.5, 0.75])]

    x0, y0 = -170, -55
    ax.text(x0, y0 + 15, "total n (radius)", ha="left", va="bottom", fontsize=9)
    for i, sv in enumerate(size_vals):
        rr = _radius_from_total(float(sv), cfg)
        ax.add_patch(Circle((x0 + i * 25, y0), rr, facecolor="none", edgecolor="0.5", linewidth=0.8))
        ax.text(x0 + i * 25, y0 - rr - 2.0, str(sv), ha="center", va="top", fontsize=8)

    plt.tight_layout()
    return fig


def save_by_continent_strategy_pies_map(
    by_continent: pd.DataFrame,
    output_path: str,
    cfg: ContinentPieConfig = ContinentPieConfig(),
) -> str:
    fig = plot_by_continent_strategy_pies(by_continent, cfg)
    fig.savefig(output_path, dpi=cfg.output_dpi, bbox_inches="tight")
    plt.close(fig)
    return output_path
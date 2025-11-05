# src/brainmaps/result.py

"""
result.py — generate a boxplot of spin-test null correlations with empirical r overlaid.

Assumes you have:
  - out/correlations.csv (or correlations_fdr.csv)
  - out/nulls/<map>.npy  (each is a 1D array of r under spins for that map)

If null .npy files are missing, this script will raise a friendly error explaining
how to save them from step_stats (see patch in section 2).
"""


import os
from pathlib import Path
from typing import Union, Iterable, List, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


def _out_dir(path: Union[str, os.PathLike] = "out") -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def load_df(out_dir: Union[str, os.PathLike] = "out", use_fdr: bool = True, alpha: float = 0.05) -> pd.DataFrame:
    out = _out_dir(out_dir)
    csv = out / ("correlations_fdr.csv" if use_fdr else "correlations.csv")
    if not csv.exists():
        raise FileNotFoundError(f"Results CSV not found: {csv}. Run `python run.py stats` or `all` first.")
    df = pd.read_csv(csv)
    df["sig"] = df["p_spin"] < alpha

    return df.sort_values("map").reset_index(drop=True)


def load_nulls(out_dir: Union[str, os.PathLike], maps: Iterable[str]) -> Dict[str, np.ndarray]:
    null_dir = _out_dir(out_dir) / "nulls"
    if not null_dir.exists():
        raise FileNotFoundError(
            f"Nulls directory not found: {null_dir}\n"
            "Make sure your stats step saves null arrays (see step_stats patch)."
        )
    out: Dict[str, np.ndarray] = {}
    for m in maps:
        fp = null_dir / f"{m}.npy"
        if not fp.exists():
            raise FileNotFoundError(
                f"Missing nulls for '{m}': {fp}\n"
                "Re-run stats after enabling null saving."
            )
        out[m] = np.load(fp)
    return out


def _as_1d_nulls(arr, map_name: str) -> np.ndarray:
    """Return a 1-D null-correlation vector or raise with a precise message."""
    a = np.asarray(arr)
    a = np.squeeze(a)

    # Accept 1-D
    if a.ndim == 1:
        return a.astype(float)

    # Accept degenerate 2-D (n_perm,1) or (1,n_perm)
    if a.ndim == 2 and 1 in a.shape:
        return a.reshape(-1).astype(float)

    # Looks like raw spins matrix, not correlations
    raise ValueError(
        f"Nulls for '{map_name}' must be a 1-D array of correlation values per "
        f"permutation, but got shape {a.shape}. This usually means you saved the "
        f"raw spun maps instead of the per-permutation correlations. "
        f"Delete out/nulls/{map_name}.npy and re-run stats to save r_null (1-D)."
    )

def boxplot_nulls_vs_empirical(
    df: pd.DataFrame,
    null_corrs: Dict[str, np.ndarray],
    outpath: Union[str, os.PathLike] = "out/figs/boxplots.png",
    ylim=(-0.8, 0.8),
):
    maps_in_order: List[str] = list(df["map"])
    null_arrays: List[np.ndarray] = [_as_1d_nulls(null_corrs[m], m) for m in maps_in_order]
    xpos = np.arange(len(maps_in_order))

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.boxplot(
        null_arrays,
        positions=xpos,
        widths=0.45,
        patch_artist=False,
        boxprops=dict(color="gray", linewidth=1),
        whiskerprops=dict(color="gray", linewidth=1),
        capprops=dict(color="gray", linewidth=1),
        medianprops=dict(color="gray", linewidth=1),
        flierprops=dict(
            marker="o", markerfacecolor="none", markeredgecolor="gray",
            markersize=3, linestyle="none"
        ),
    )

    for xi, row in df.iterrows():
        ax.scatter(
            xi,
            float(row["r"]),
            s=80,
            color=("red" if bool(row["sig"]) else "#d58f6c"),
            zorder=3,
        )

    ax.set_ylabel("Pearson's r")
    ax.set_xticks(xpos)
    ax.set_xticklabels(maps_in_order, rotation=45, ha="right")
    ax.set_ylim(*ylim)

    emp_ns = mlines.Line2D([], [], color="#d58f6c", marker="o", linestyle="None", markersize=8,
                           label="Empirical (P_spin ≥ 0.05)")
    emp_sig = mlines.Line2D([], [], color="red", marker="o", linestyle="None", markersize=8,
                            label="Empirical (P_spin < 0.05)")
    null_leg = mlines.Line2D([], [], color="gray", marker="s", fillstyle="none", linestyle="-", markersize=8,
                             label="Spatial null")
    ax.legend(handles=[emp_ns, emp_sig, null_leg],
              loc="upper left", bbox_to_anchor=(0, -0.3), ncol=3, frameon=False)

    Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"[result] wrote {outpath}")

# src/brainmaps/boxplot.py

'''
boxplot.py — generate a boxplot of spin-test null correlations with empirical r overlaid.

Assumes you have:
  - out/correlations.csv (or correlations_fdr.csv)
  - out/nulls/<map>.npy  (each is a 1D array of r under spins for that map)

If null .npy files are missing, this script will raise a friendly error explaining
how to save them from step_stats (see patch in section 2).
'''


import os
from pathlib import Path
from typing import Union, Iterable, List, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches


def _out_dir(path: Union[str, os.PathLike] = "out") -> Path:
    '''
    Ensures an output directory exists and returns it as a Path object.

    :param path: Directory path to create or verify (default: "out").
    :returns: A pathlib.Path object pointing to the ensured directory.
    '''
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def load_df(out_dir: Union[str, os.PathLike] = "out", use_fdr: bool = True, alpha: float = 0.05) -> pd.DataFrame:
    '''
    Loads the correlation results from CSV and marks significant entries.

    :param out_dir: Base output directory where CSVs are stored (default: "out").
    :param use_fdr: If True, loads 'correlations_fdr.csv'; otherwise, 'correlations.csv'.
    :param alpha: Significance threshold for marking p_spin < alpha (default: 0.05).
    :returns: A pandas DataFrame sorted by 'map' with an added 'sig' column.
    :raises FileNotFoundError: If the required CSV file does not exist.
    '''
    out = _out_dir(out_dir)
    csv = out / ("correlations_fdr.csv" if use_fdr else "correlations.csv")
    if not csv.exists():
        raise FileNotFoundError(f"Results CSV not found: {csv}. Run `python run.py stats` or `all` first.")
    df = pd.read_csv(csv)
    df["sig"] = df["p_spin"] < alpha

    return df.sort_values("map").reset_index(drop=True)


def load_nulls(out_dir: Union[str, os.PathLike], maps: Iterable[str]) -> Dict[str, np.ndarray]:
    '''
    Loads per-map null correlation arrays from disk.

    :param out_dir: Base output directory containing the 'nulls' subdirectory.
    :param maps: Iterable of map names for which to load null correlation arrays.
    :returns: A dictionary mapping each map name to a 1D NumPy array of null correlations.
    :raises FileNotFoundError: If the 'nulls' directory or any required file is missing.
    '''
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
    '''
    Validates and coerces a null correlation array to a 1D float vector.

    :param arr: The array to validate; may be 1D or a degenerate 2D shape.
    :param map_name: The name of the map (used in error messages).
    :returns: A 1D NumPy array of float correlation values.
    :raises ValueError: If the array does not resemble a vector of correlations
                        (e.g., raw spins matrix instead of per-permutation r).
    '''
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
    ylim=(-1.0, 1.0),
):
    """
    Creates a boxplot comparing spatial spin-test null distributions
    with empirical Pearson correlations for each target map, and
    overlays past-work ranges/points.

    :param df: DataFrame with columns ['map', 'r', 'p_spin', 'p_fdr', 'sig' (optional)].
    :param null_corrs: Dict mapping technical map name to 1D array of null correlation values.
    :param outpath: Path to save the output PNG (default: "out/figs/boxplots.png").
    :param ylim: Tuple defining y-axis limits for the boxplot (default: (-1.0, 1.0)).
    :returns: A tuple of (fig, ax) — the matplotlib Figure and Axes.
    """

    alpha = 0.05

    # Mapping of short map names 
    map_labels = {
        "genepc1": "PC1 gene Expression",
        "myelin": "T1w/T2w ratio",
        "devexp": "Developmental expansion",
        "evoexp": "Evolutionary expansion",
        "gradient_pc1": "Functional gradient",
        "isv": "Intersubject variability",
        "cbf": "Cerebral blood flow",
        "cbv": "Cerebral blood volume",
        "cmro2": "Oxygen metabolism",
        "cmrglc": "Glucose metabolism",
        "scalingnih": "Allometric scaling (NIH)",
        "scalingpnc": "Allometric scaling (PNC)",
    }

    # Desired publication order (display labels)
    desired_order = [
        "PC1 gene expression",
        "T1w/T2w ratio",
        "Oxygen metabolism",
        "Functional gradient",
        "Glucose metabolism",
        "Intersubject variability",
        "Cerebral blood flow",
        "Evolutionary expansion",
        "Allometric scaling (NIH)",
        "Cerebral blood volume",
        "Allometric scaling (PNC)",
        "Developmental expansion",
    ]
    order_index = {name.lower(): i for i, name in enumerate(desired_order)}

    # Work on a copy so we don't mutate the caller's df in-place
    df = df.copy()

    # Choose which p to use for significance
    if "p_fdr" not in df.columns:
        raise ValueError("DataFrame must contain a 'p_fdr' column for FDR-corrected p-values.")
    df["sig"] = df["p_fdr"] < alpha

    maps_in_order: List[str] = list(df["map"])

    # First pass: get display labels for each row
    display_labels_pre = [
        map_labels.get(m.replace("_fsLR32k", ""), m) for m in maps_in_order
    ]
    display_labels_pre_lower = [lab.lower() for lab in display_labels_pre]

    # Reorder to match desired publication order
    sort_idx = sorted(
        range(len(maps_in_order)),
        key=lambda i: order_index.get(
            display_labels_pre_lower[i],
            len(desired_order) + i,  # fallback: keep relative order if not found
        ),
    )

    maps_in_order = [maps_in_order[i] for i in sort_idx]
    df = df.iloc[sort_idx].reset_index(drop=True)

    display_labels = [
        map_labels.get(m.replace("_fsLR32k", ""), m) for m in maps_in_order
    ]
    display_labels_lower = [lab.lower() for lab in display_labels]

    # Prepare null arrays in the same (reordered) order
    null_arrays: List[np.ndarray] = [
        _as_1d_nulls(null_corrs[m], m) for m in maps_in_order
    ]
    xpos = np.arange(len(maps_in_order))

    # --- Plot setup ---
    fig, ax = plt.subplots(figsize=(10, 4))

    # gray boxplots for null distributions
    ax.boxplot(
        null_arrays,
        positions=xpos,
        widths=0.4,
        patch_artist=False,
        boxprops=dict(color="gray", linewidth=1),
        whiskerprops=dict(color="gray", linewidth=1),
        capprops=dict(color="gray", linewidth=1),
        medianprops=dict(color="gray", linewidth=1),
        flierprops=dict(
            marker="o",
            markerfacecolor="none",
            markeredgecolor="gray",
            markersize=3,
            linestyle="none",
        ),
    )

    # overlay empirical r as dots
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
    ax.set_xticklabels(display_labels, rotation=45, ha="right")
    ax.set_ylim(*ylim)

    # current y limits (after setting ylim)
    ymin, ymax = ax.get_ylim()

    # shade negative region for PC1 gene expression (sign only)
    pc1_label = "PC1 gene expression"  # match ignoring case
    if pc1_label.lower() in display_labels_lower:
        pc1_idx = display_labels_lower.index(pc1_label.lower())
        shade_pc1 = mpatches.Rectangle(
            (pc1_idx - 0.4, ymin),   # bottom left (x, y)
            width=0.8,               # span of that category
            height=0 - ymin,         # up to zero
            facecolor="gray",
            alpha=0.2,
            edgecolor="none",
            zorder=0,                # behind boxplots and dots
        )
        ax.add_patch(shade_pc1)

    # shade past work range for oxygen and glucose metabolism, etc. 
    past_ranges = {
        "Oxygen metabolism": (-1.0, 0),
        "Glucose metabolism": (-0.71, -0.24),
        "T1w/T2w ratio": (-1.0, 0),
        "Cerebral blood flow": (-1.0, 0),
        "Cerebral blood volume": (-1.0, 0),
    }

    for map_name, (y0, y1) in past_ranges.items():
        key = map_name.lower()
        if key in display_labels_lower:
            xi = display_labels_lower.index(key)
            shade_range = mpatches.Rectangle(
                (xi - 0.35, y0),
                width=0.7,
                height=y1 - y0,
                facecolor="gray",
                alpha=0.2,
                edgecolor="none",
                zorder=0.5,          # behind dots but above background
            )
            ax.add_patch(shade_range)

    # dictionary of past work correlations (black diamonds)
    past_work_r = {
        "Intersubject variability": 0.05,
    }

    for map_name, r_prev in past_work_r.items():
        key = map_name.lower()
        if key in display_labels_lower:
            xi = display_labels_lower.index(key)
            ax.scatter(
                xi,
                r_prev,
                s=90,
                marker="D",   # diamond
                color="gray",
                zorder=4,
            )

    # past work with no relationship (new marker: x)
    no_relation_maps = [
        "Developmental expansion",
        "Evolutionary expansion",
        "Functional gradient",
        "Allometric scaling (NIH)",
        "Allometric scaling (PNC)",
    ]

    for map_name in no_relation_maps:
        key = map_name.lower()
        if key in display_labels_lower:
            xi = display_labels_lower.index(key)
            ax.scatter(
                xi,
                0.0,                 
                s=90,
                marker="x",
                color="gray",
                zorder=4,
            )

    # Legend 
    emp_ns = mlines.Line2D(
        [], [], color="#d58f6c", marker="o", linestyle="None", markersize=8,
        label="Empirical (p_spin ≥ 0.05)",
    )
    emp_sig = mlines.Line2D(
        [], [], color="red", marker="o", linestyle="None", markersize=8,
        label="Empirical (p_spin < 0.05)",
    )
    null_leg = mpatches.Patch(
        facecolor="none", edgecolor="gray",
        label="Spatial null",
    )
    past_leg = mpatches.Patch(
        facecolor="gray", edgecolor="none", alpha=0.2,
        label="Past work range",
    )
    prev_leg = mlines.Line2D(
        [], [], color="gray", marker="D", linestyle="None", markersize=8,
        label="Past work",
    )
    no_rel_leg = mlines.Line2D(
        [], [], color="gray", marker="x", linestyle="None", markersize=8,
        label="Past work (no relationship)",
    )

    ax.legend(
        handles=[emp_ns, emp_sig, null_leg, past_leg, prev_leg, no_rel_leg],
        loc="upper center",
        bbox_to_anchor=(0.5, 1.28),
        ncol=3,
        frameon=False,
    )

    # --- Save ---
    plt.subplots_adjust(top=0.8, bottom=0.3)
    Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"[result] wrote {outpath}")

    return fig, ax

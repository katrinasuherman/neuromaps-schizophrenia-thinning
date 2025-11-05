'''
stats.py contains utilities for computing correlations, spin-test nulls,
and FDR correction for brain-map comparisons.
'''

import numpy as np
import pandas as pd
from neuromaps.stats import compare_images
from neuromaps.nulls import alexander_bloch
from statsmodels.stats.multitest import multipletests


def pearson_nan_safe(a, b):
    '''
    Compute the Pearson correlation (r) between two brain maps using
    neuromaps.compare_images for consistent masking/NaN handling.

    :param a: 1D NumPy array (or image-like supported by compare_images) for map A.
    :param b: 1D NumPy array (or image-like supported by compare_images) for map B.
    :return: Float Pearson correlation coefficient r (rounded to 6 decimals).
    '''
    r = compare_images(a, b, metric='pearsonr', ignore_zero=True)
    return round(float(r), 6)


# def spin_test(src_vec, tgt_vec, atlas="fsLR", density="32k", n_perm=1000, seed=42):
#     '''
#     Compute observed Pearson r and spin-test p-value using neuromaps.
#     :param src_vec: 1D NumPy array for source map.
#     :param tgt_vec: 1D NumPy array for target map.
#     :param atlas: Surface atlas name.
#     :param density: Surface mesh density.
#     :param n_perm: Number of permutations.
#     :param seed: Random seed.
#     :return: (r_obs, p_spin, r_null_array)
#     '''
#     # Clean NaNs/Infs
#     src_vec = np.asarray(src_vec, dtype=float)
#     tgt_vec = np.asarray(tgt_vec, dtype=float)
#     src_vec[~np.isfinite(src_vec)] = np.nan
#     tgt_vec[~np.isfinite(tgt_vec)] = np.nan

#     # Generate null maps
#     nulls = alexander_bloch(
#         data=src_vec,
#         atlas=atlas,
#         density=density,
#         n_perm=n_perm,
#         seed=seed
#     )

#     # Observed correlation (with ignore_zero=True)
#     r_obs = compare_images(src_vec, tgt_vec, metric="pearsonr", ignore_zero=True)

#     # Try the fast path (if compare_images supports nulls)
#     try:
#         r_obs2, p_val = compare_images(src_vec, tgt_vec, metric="pearsonr", nulls=nulls, ignore_zero=True)
#         r_obs = round(float(r_obs2), 6)
#         p_spin = round(float(p_val), 6)
#         r_null = nulls
#     except TypeError:
#         # Fallback: manual computation
#         if nulls.shape[0] == nulls.shape[1]:  # (n_vertices, n_perm)
#             r_null = np.array([
#                 compare_images(nulls[:, i], tgt_vec, metric="pearsonr", ignore_zero=True)
#                 for i in range(n_perm)
#             ])
#         else:  # (n_perm, n_vertices)
#             r_null = np.array([
#                 compare_images(nulls[i, :], tgt_vec, metric="pearsonr", ignore_zero=True)
#                 for i in range(n_perm)
#             ])
#         p_spin = (np.sum(np.abs(r_null) >= abs(r_obs)) + 1) / (n_perm + 1)
#         r_obs = round(float(r_obs), 6)
#         p_spin = round(float(p_spin), 6)

#     return r_obs, p_spin, r_null.astype(float)

def spin_test(src_vec, tgt_vec, atlas="fsLR", density="32k", n_perm=1000, seed=42):
    """
    Compute observed Pearson r and a spin-test p-value.

    Returns:
        r_obs (float): observed Pearson correlation (rounded to 6 d.p.)
        p_spin (float): two-tailed spin-test p-value (rounded to 6 d.p.)
        r_null (np.ndarray): 1-D array of null correlations, shape (n_perm,)
    """
    # Clean NaNs/Infs
    src_vec = np.asarray(src_vec, dtype=float)
    tgt_vec = np.asarray(tgt_vec, dtype=float)
    src_vec[~np.isfinite(src_vec)] = np.nan
    tgt_vec[~np.isfinite(tgt_vec)] = np.nan

    # Generate raw spun maps (vertices x permutations OR permutations x vertices)
    nulls = alexander_bloch(
        data=src_vec, atlas=atlas, density=density,
        n_perm=n_perm, seed=seed
    )

    # Observed correlation
    r_obs = compare_images(src_vec, tgt_vec, metric="pearsonr", ignore_zero=True)

    # Convert raw spins -> per-permutation correlation vector (always 1-D)
    n_vert = src_vec.size
    if nulls.ndim != 2:
        raise ValueError(f"Unexpected nulls ndim={nulls.ndim}; expected 2-D spin matrix.")
    if nulls.shape[0] == n_vert:
        # (n_vertices, n_perm)
        it = (nulls[:, i] for i in range(nulls.shape[1]))
    elif nulls.shape[1] == n_vert:
        # (n_perm, n_vertices)
        it = (nulls[i, :] for i in range(nulls.shape[0]))
    else:
        raise ValueError(
            f"Nulls shape {nulls.shape} is incompatible with source length {n_vert}."
        )

    r_null = np.array(
        [compare_images(nv, tgt_vec, metric="pearsonr", ignore_zero=True) for nv in it],
        dtype=float
    )
    if r_null.ndim != 1 or r_null.size != n_perm:
        # sanity check in case backend changes
        r_null = np.asarray(r_null, dtype=float).reshape(-1)

    # Two-tailed p with +1 smoothing
    p_spin = (np.sum(np.abs(r_null) >= abs(r_obs)) + 1) / (r_null.size + 1)

    # Round only the scalars; keep r_null full precision for downstream plots
    return round(float(r_obs), 6), round(float(p_spin), 6), r_null

def fdr(df, alpha=0.05):
    '''
    Apply Benjamini–Hochberg FDR correction to spin-test p-values.

    :param df: Pandas DataFrame with a 'p_spin' column.
    :param alpha: FDR control level (default 0.05).
    :return: DataFrame with 'p_fdr' and 'sig_fdr', rounded to 6 decimals.
    '''
    df = df.copy()
    reject, p_fdr, _, _ = multipletests(df["p_spin"].values, alpha=alpha, method="fdr_bh")
    df["p_fdr"] = np.round(p_fdr, 6)
    df["sig_fdr"] = reject
    return df.sort_values("p_fdr")

def nulls_to_corrs(nulls: np.ndarray, tgt_vec: np.ndarray) -> np.ndarray:
    """
    Convert a nulls matrix of spun maps to a 1-D vector of Pearson r
    against tgt_vec (one r per permutation). Accepts either
    (n_vertices, n_perm) or (n_perm, n_vertices).
    """
    nulls = np.asarray(nulls, float)
    tgt_vec = np.asarray(tgt_vec, float)

    if nulls.ndim == 1:
        return nulls  # already correlation vector

    if nulls.ndim != 2:
        raise ValueError(f"Unexpected nulls ndim={nulls.ndim}")

    # choose iteration over permutations robustly
    if nulls.shape[0] == tgt_vec.shape[0]:
        it = (nulls[:, i] for i in range(nulls.shape[1]))
    elif nulls.shape[1] == tgt_vec.shape[0]:
        it = (nulls[i, :] for i in range(nulls.shape[0]))
    else:
        raise ValueError(
            f"Nulls shape {nulls.shape} not compatible with target length {tgt_vec.shape[0]}"
        )

    r_null = np.array([
        compare_images(nv, tgt_vec, metric="pearsonr", ignore_zero=True)
        for nv in it
    ], dtype=float)
    return r_null
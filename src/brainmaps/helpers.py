'''
helpers.py contains functions for loading and normalizing brain map data
into left (L) and right (R) hemisphere arrays.
'''

from neuromaps.datasets import fetch_annotation, fetch_atlas, ALIAS
from neuromaps.images import load_gifti
import nibabel as nib
import numpy as np, tempfile, os
from pathlib import Path


def data1d(x):
    '''
    Load a GIFTI image or file path and return a 1D NumPy array.

    :param x: Path to a .gii file or a preloaded nibabel GiftiImage object.
    :return: A 1D NumPy array of vertex data values.
    '''
    img = x if isinstance(x, nib.gifti.gifti.GiftiImage) else nib.load(x)
    return np.asarray(img.agg_data()).squeeze()


def lr_arrays(entry, nL=32492, nR=32492, fill='nan', single='R', L_mean=None):
    '''
    Normalize an input into (L, R) hemisphere arrays.

    :param entry: Input data. Can be:
        - A tuple/list of length 2: (L, R) paths or GiftiImage objects.
        - A tuple/list of length 1: (hemi_only,) where hemi is given by `single`.
        - A single path or GiftiImage (hemi-only).
    :param nL: Number of vertices in the left hemisphere (default 32492).
    :param nR: Number of vertices in the right hemisphere (default 32492).
    :param fill: Strategy if only one hemisphere is provided.
        'nan'    → pad the missing hemisphere with NaNs (default)
        'mirror' → duplicate the same array for both hemispheres.
    :param single: Which hemisphere a single array represents: 'L' or 'R' (default 'R').
    :param L_mean: Placeholder parameter (unused).
    :return: (L_arr, R_arr) as 1D NumPy arrays.
    '''
    if isinstance(entry, (tuple, list)):
        # Case 1: both hemispheres provided
        if len(entry) == 2:
            L, R = entry
            L_arr = data1d(L) if L is not None else np.full(nL, np.nan)
            R_arr = data1d(R) if R is not None else np.full(nR, np.nan)
            return L_arr, R_arr
        
        # Case 2: single hemisphere only (not a tuple/list)
        elif len(entry) == 1:
            arr = data1d(entry[0])
            if fill == 'mirror':
                return arr.copy(), arr.copy()
            # pad the other hemisphere with NaNs
            if str(single).upper() == 'L':
                return arr, np.full(nR, np.nan)
            else:
                return np.full(nL, np.nan), arr
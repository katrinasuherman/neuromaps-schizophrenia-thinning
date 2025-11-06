"""
transforms.py contains helper functions for converting brain maps between
different surface coordinate systems and densities using the Neuromaps
transformation utilities.

Each function wraps a Neuromaps transform to resample or reproject data into
the standard fsLR-32k surface space, which is used throughout this project
for cross-map correlation and visualization.
"""

from neuromaps.transforms import fsaverage_to_fslr, fslr_to_fslr, civet_to_fslr

def to_fsLR32k_from_fsaverage10k(annotation):
    """
    Transform a brain annotation from fsaverage-10k space to fsLR-32k.

    :param annotation: Path or GiftiImage object representing the input
        annotation in fsaverage-10k space.
    :return: A tuple of transformed (L, R) hemisphere GiftiImage objects
        in fsLR-32k space.
    """
    return fsaverage_to_fslr(annotation, target_density="32k")

def to_fsLR32k_from_fsLR164k(annotation, hemi=None):
    """
    Downsample a high-resolution fsLR-164k annotation to fsLR-32k.

    :param annotation: Path or GiftiImage object representing the input
        annotation in fsLR-164k space.
    :param hemi: Optional string specifying the hemisphere ('L' or 'R')
        if the input contains only one hemisphere.
    :return: A tuple of transformed (L, R) hemisphere GiftiImage objects
        in fsLR-32k space, or a single GiftiImage if only one hemisphere
        was provided.
    """
    return fslr_to_fslr(annotation, target_density="32k", hemi=hemi)

def to_fsLR32k_from_civet41k(annotation):
    """
    Transform a CIVET-41k surface annotation to fsLR-32k space.

    :param annotation: Path or GiftiImage object representing the input
        annotation in CIVET-41k surface space.
    :return: A tuple of transformed (L, R) hemisphere GiftiImage objects
        in fsLR-32k space.
    """
    return civet_to_fslr(annotation, target_density="32k")

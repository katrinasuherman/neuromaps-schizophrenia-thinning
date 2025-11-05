from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
import matplotlib.colors as mcolors
import numpy as np
from nilearn.plotting import plot_surf
from neuromaps.datasets import ALIAS, fetch_atlas
from neuromaps.images import load_gifti
from neuromaps.transforms import _check_hemi

HEMI = dict(L='left', R='right')

def plot_surf_lateral_only(
    data, template, density, *,
    surf='inflated', hemi=None, mask_medial=True,
    colorbar='shared',          # 'none', 'each', or 'shared'
    cbar_location='right',      # 'right' or 'bottom' when colorbar='shared'
    vmin=None, vmax=None,       # shared limits. If None they are computed
    wspace=0.12, figsize=(10, 4), **kwargs
):
    '''
    Create side-by-side lateral surface plots (one per hemisphere) with optional shared colorbar.

    :param data: Surface data to plot. Accepts:
        - Iterable of (img, hemi) pairs understood by neuromaps `_check_hemi`
        - Path(s) to GIFTI file(s) or preloaded GIFTI images
        - A single hemisphere with `hemi="L"` or `"R"` or auto-detected by `_check_hemi`.
    :param template: Surface template key (e.g., 'fsLR', 'fsaverage'; not 'MNI152').
    :param density: Surface density string for the template (e.g., '32k', '164k').
    :param surf: Surface geometry to use from the atlas (default 'inflated').
    :param hemi: Hemisphere selector if `data` doesn’t encode it. One of {'L','R'} or None.
    :param mask_medial: If True, mask out non-cortical/medial wall vertices before plotting.
    :param colorbar: Colorbar mode. One of:
        - 'shared' (default): one colorbar for all panels
        - 'each': per-panel colorbars
        - 'none': no colorbar
    :param cbar_location: Location for shared colorbar; 'right' (default) or 'bottom'.
    :param vmin: Lower bound for colormap. If None, computed from finite values across hemispheres.
    :param vmax: Upper bound for colormap. If None, computed from finite values across hemispheres.
    :param wspace: Horizontal spacing between hemisphere subplots.
    :param figsize: Matplotlib figure size tuple (width, height).
    :param kwargs: Extra keyword args forwarded to `nilearn.plotting.plot_surf`
        (e.g., cmap, bg_map, alpha, threshold, etc.).

    :return: (fig, cbar)
        - fig: The Matplotlib Figure containing the hemisphere plots.
        - cbar: The colorbar object when `colorbar='shared'`, otherwise None.
    '''
    atlas = fetch_atlas(template, density, verbose=0)
    template = ALIAS.get(template, template)
    if template == 'MNI152':
        raise ValueError("Cannot plot MNI152 on the surface. Project to a surface first.")

    surf_geom = atlas[surf]
    medial = atlas['medial']

    # nilearn defaults
    opts = dict(alpha=1.0, threshold=np.spacing(1))
    opts.update(**kwargs)
    if kwargs.get('bg_map') is not None and kwargs.get('alpha') is None:
        opts['alpha'] = 'auto'

    data, hemispheres = zip(*_check_hemi(data, hemi))
    n_hemi = len(data)

    # collect data arrays to compute shared limits if needed
    arrays = []
    for img in data:
        arr = load_gifti(img).agg_data().astype('float32')
        if mask_medial:
            # we do not know hemisphere here yet, compute after
            arrays.append(arr)
        else:
            arrays.append(arr)

    if vmin is None or vmax is None:
        # robust range across hemispheres
        allvals = np.concatenate([np.ravel(a[np.isfinite(a)]) for a in arrays])
        vmin = np.nanmin(allvals) if vmin is None else vmin
        vmax = np.nanmax(allvals) if vmax is None else vmax
    opts.update(dict(vmin=vmin, vmax=vmax))

    fig, axes = plt.subplots(1, n_hemi, subplot_kw={'projection': '3d'}, figsize=figsize)
    axes = (axes,) if n_hemi == 1 else axes

    # if we will add a shared colorbar, turn per-panel colorbars off
    per_panel_cbar = (colorbar == 'each')
    if colorbar == 'shared':
        opts['colorbar'] = False
    else:
        opts['colorbar'] = per_panel_cbar

    for ax, h, img in zip(axes, hemispheres, data):
        geom = load_gifti(getattr(surf_geom, h)).agg_data()
        arr = load_gifti(img).agg_data().astype('float32')
        if mask_medial:
            med = load_gifti(getattr(medial, h)).agg_data().astype(bool)
            arr[~med] = np.nan

        ax.disable_mouse_rotation()
        plot_surf(geom, arr, hemi=HEMI[h], axes=ax, view='lateral', **opts)

    # spacing between hemispheres
    fig.subplots_adjust(wspace=wspace)

    # optional shared colorbar
    cbar = None
    if colorbar == 'shared':
        sm = ScalarMappable(norm=mcolors.Normalize(vmin=vmin, vmax=vmax),
                            cmap=opts.get('cmap', None))
        sm.set_array([])
        if cbar_location == 'right':
            cbar = fig.colorbar(sm, ax=axes, location='right', fraction=0.05, pad=0.03)
        else:  # bottom
            cbar = fig.colorbar(sm, ax=axes, location='bottom', fraction=0.05, pad=0.3)

    return fig, cbar

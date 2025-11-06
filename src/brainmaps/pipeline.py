# src/brainmaps/pipeline.py

'''
pipeline.py orchestrates the end-to-end workflow:
environment check → fetch/transform maps → correlations/spin test → FDR → figures.
All outputs are written under cfg.out_dir (default: "out/").
'''

import os
from pathlib import Path
from neuromaps.images import load_gifti
from neuromaps.datasets import fetch_annotation
import numpy as np
import pandas as pd

from .config import Config
from . import env
from . import boxplot
from . import catalog
# from .catalog import get_plot_kwargs
from . import transforms as T
from . import plotting
from . import helpers
from . import stats


def _ensure_outdir(cfg: Config) -> Path:
    '''
    Ensure the output directory exists.

    :param cfg: Config object with cfg.out_dir.
    :return: Path to the output directory.
    '''
    out = Path(cfg.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    return out

def step_env(cfg: Config):
    '''
    Optional environment checks (e.g., Connectome Workbench presence).

    :param cfg: Config object.
    :return: None
    '''
    try:
        version = env.ensure_workbench()  # prints wb_command -version
        print(f"[env] workbench OK: {version}")
    except Exception as e:
        # Not fatal for many neuromaps ops, but helpful to know.
        print(f"[env] warning: {e}")

def step_transforms(cfg: Config):
    '''
    Fetch the source map and all target maps, transforming each to fsLR-32k as needed.

    :param cfg: Config object.
    :return: dict mapping map_name -> data, plus:
             {"__source_name__": str, "__source__": (L_img, R_img) or single hemi only}
    '''

    def _fetch(spec: dict):
        # pass only the kwargs fetch_annotation accepts
        kwargs = {k: spec[k] for k in ("source", "desc", "space", "den") if k in spec}
        return fetch_annotation(**kwargs)

    out = {}

    # source
    src = catalog.get_source()       # if this already returns a fetched annotation, keep it
    out["__source_name__"] = "source_thickness"
    out["__source__"] = src

    # targets
    targets = catalog.get_targets()
    for name, spec in targets.items():
        space = spec.get("space")
        den = spec.get("den")
        hemi = spec.get("hemi")  # may be None

        if space == "fsLR" and den == "32k":
            out[name] = _fetch(spec)

        elif space == "fsaverage":
            out[name] = T.to_fsLR32k_from_fsaverage10k(_fetch(spec))

        elif space == "fsLR" and den == "164k":
            out[name] = T.to_fsLR32k_from_fsLR164k(_fetch(spec), hemi=hemi)

        elif space == "civet":
            out[name] = T.to_fsLR32k_from_civet41k(_fetch(spec))

        else:
            raise NotImplementedError(
                f"Transform not defined for {name}: space={space}, den={den}"
            )

    return out

def _vectors_from_maps(maps_dict: dict) -> tuple[str, np.ndarray, dict]:
    '''
    Build concatenated LH||RH vectors for source and targets.

    :param maps_dict: dict from step_transforms().
    :return: (src_name, src_vec, vectors) where:
             - src_name is the key for the source
             - src_vec is a 1D NumPy vector (LH||RH, NaN where missing)
             - vectors is {name: 1D vector} for each target
    '''
    src_name = maps_dict["__source_name__"]
    src_data = maps_dict["__source__"]

    # Source vector
    L, R = helpers.lr_arrays(src_data, fill="nan")
    src_vec = np.concatenate([L, R], 0)

    # Target vectors
    vectors = {}
    for name, data in maps_dict.items():
        if name.startswith("__"):
            continue
        L, R = helpers.lr_arrays(data, fill="nan")
        vec = np.concatenate([L, R], 0)
        vectors[name] = vec

    return src_name, src_vec, vectors

def step_stats(cfg: Config, maps_dict: dict) -> pd.DataFrame:
    '''
    Compute observed Pearson r (via compare_images) and spin-test p-values
    for each target against the source; write CSV to out/correlations.csv.

    :param cfg: Config object (uses n_perm and seed).
    :param maps_dict: dict produced by step_transforms().
    :return: Pandas DataFrame with columns ['map', 'r', 'p_spin'] sorted by p_spin.
    '''
    outdir = _ensure_outdir(cfg)
    null_dir = outdir / "nulls"
    null_dir.mkdir(parents=True, exist_ok=True)
    src_name, src_vec, vectors = _vectors_from_maps(maps_dict)

    rows = []
    for name, tgt_vec in vectors.items():
        # observed r and spin-test
        r_obs, p_spin, r_null = stats.spin_test(
            src_vec, tgt_vec,
            atlas="fsLR", density="32k",
            n_perm=cfg.n_perm, seed=cfg.seed
        )
        rows.append({"map": name, "r": float(r_obs), "p_spin": float(p_spin)})

        if r_null.ndim != 1:
            print(f"[stats] converting raw spins to r_null for {name} from shape {r_null.shape}")
            r_null = stats.nulls_to_corrs(r_null, tgt_vec)

        # save null correlations for the boxplot
        np.save(null_dir / f"{name}.npy", np.asarray(r_null, float).reshape(-1))

    df = pd.DataFrame(rows).sort_values("p_spin")
    (outdir / "correlations.csv").write_text(df.to_csv(index=False))
    print(f"[stats] wrote {outdir/'correlations.csv'}")
    print(f"[stats] wrote nulls to {null_dir}")
    return df

def step_fdr(cfg: Config, df_stats: pd.DataFrame) -> pd.DataFrame:
    '''
    Apply BH-FDR to the spin-test p-values and save to out/correlations_fdr.csv.

    :param cfg: Config object.
    :param df_stats: DataFrame from step_stats().
    :return: DataFrame with added 'p_fdr' and 'sig_fdr' columns, sorted by p_fdr.
    '''
    outdir = _ensure_outdir(cfg)
    df_adj = stats.fdr(df_stats, alpha=0.05)
    (outdir / "correlations_fdr.csv").write_text(df_adj.to_csv(index=False))
    print(f"[fdr] wrote {outdir/'correlations_fdr.csv'}")
    return df_adj

def _per_hemi_length(img) -> int:
    """Return number of vertices in the hemi data array for an image/path."""
    arr = load_gifti(img).agg_data()
    arr = np.asarray(arr).squeeze()
    return int(arr.shape[0])

def _infer_template_density(img) -> tuple[str, str]:
    """Infer (template, density) from per-hemi vertex count."""
    n = _per_hemi_length(img)
    # Common cases:
    if n == 32492:   return ("fsLR", "32k")
    if n == 10242:   return ("fsaverage", "10k")
    if n == 40962:   return ("civet", "41k")
    if n == 163842:  return ("fsLR", "164k")
    # Fallback: assume fsLR-32k
    return ("fsLR", "32k")

def step_viz(cfg: Config, maps_dict: dict, *, save_png=True):
    from matplotlib import pyplot as plt
    from .catalog import get_plot_kwargs

    outdir = _ensure_outdir(cfg)
    figs_dir = outdir / "figs"
    if save_png:
        figs_dir.mkdir(parents=True, exist_ok=True)

    def _extract_imgs(data):
        return [x for x in (data if isinstance(data, (tuple, list)) else [data]) if x is not None]

    def _guess_hemi_from_name(obj) -> str:
        s = str(obj).lower()
        if any(t in s for t in ("lh.", "left", ".l.", "_l.", "-l.")): return "L"
        if any(t in s for t in ("rh.", "right", ".r.", "_r.", "-r.")): return "R"
        return "R"

    def _apply_inferred_template_density(imgs, kwargs, label):
        if not imgs:
            return kwargs
        tpl_inf, den_inf = _infer_template_density(imgs[0])
        tpl = kwargs.get("template")
        den = kwargs.get("density")
        if (tpl, den) != (tpl_inf, den_inf):
            print(f"[viz] note: overriding template/density for {label} "
                  f"from ({tpl}, {den}) to ({tpl_inf}, {den_inf}) based on data.")
            kwargs["template"], kwargs["density"] = tpl_inf, den_inf
        return kwargs

    # --- source (before + transformed if present in specs) ---
    src_name = maps_dict["__source_name__"]
    src_imgs = _extract_imgs(maps_dict["__source__"])

    for plot_key, out_suffix in [(src_name, ""), (f"{src_name}_fsLR32k", "_fsLR32k")]:
        kwargs = get_plot_kwargs(plot_key)
        if not kwargs and out_suffix:  # no transformed spec; skip second pass
            continue
        title = kwargs.pop("title", plot_key)
        if len(src_imgs) == 1 and "hemi" not in kwargs:
            kwargs["hemi"] = _guess_hemi_from_name(src_imgs[0])
        kwargs = _apply_inferred_template_density(src_imgs, kwargs, plot_key)

        fig, _ = plotting.plot_surf_lateral_only(data=src_imgs, **kwargs)
        if save_png:
            fp = figs_dir / f"{src_name}{out_suffix}.png"
            fig.suptitle(title, y=0.98)
            fig.savefig(fp, bbox_inches="tight", dpi=150)
            plt.close(fig)
            print(f"[viz] wrote {fp}")

    # --- targets (before + transformed) ---
    for name, data in maps_dict.items():
        if name.startswith("__"):
            continue

        imgs = _extract_imgs(data)
        if not imgs:
            print(f"[viz] skip {name}: no images")
            continue

        for plot_key, out_suffix in [(name, ""), (f"{name}_fsLR32k", "_fsLR32k")]:
            kwargs = get_plot_kwargs(plot_key)
            if not kwargs and out_suffix:
                # Only plot transformed if a spec exists for *_fsLR32k
                continue

            title = kwargs.pop("title", plot_key)
            if len(imgs) == 1 and "hemi" not in kwargs:
                kwargs["hemi"] = _guess_hemi_from_name(imgs[0])
            kwargs = _apply_inferred_template_density(imgs, kwargs, plot_key)

            fig, _ = plotting.plot_surf_lateral_only(data=imgs, **kwargs)
            if save_png:
                fp = figs_dir / f"{name}{out_suffix}.png"
                fig.suptitle(title, y=0.98)
                fig.savefig(fp, bbox_inches="tight", dpi=150)
                plt.close(fig)
                print(f"[viz] wrote {fp}")

def step_result(cfg: Config, use_fdr: bool = True):
    """
    Generate the summary boxplot of null distributions vs. empirical correlations.
    """
    df = boxplot.load_df(out_dir=cfg.out_dir, use_fdr=use_fdr)
    nulls = boxplot.load_nulls(out_dir=cfg.out_dir, maps=df["map"])
    boxplot.boxplot_nulls_vs_empirical(
        df, nulls, outpath=Path(cfg.out_dir) / "figs" / "boxplots.png"
    )

def run_all(cfg: Config) -> pd.DataFrame:
    '''
    Convenience function that runs transforms → stats → FDR → viz.

    :param cfg: Config object.
    :return: FDR-adjusted DataFrame (sorted by p_fdr).
    '''
    step_env(cfg)
    maps_dict = step_transforms(cfg)
    df = step_stats(cfg, maps_dict)
    df_fdr = step_fdr(cfg, df)
    step_viz(cfg, maps_dict, save_png=True)
    step_result(cfg)
    return df

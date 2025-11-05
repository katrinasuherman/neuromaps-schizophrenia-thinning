from neuromaps.datasets import fetch_annotation


PLOT_DEFAULTS = {
    "template": "fsLR",
    "density": "32k",
    "colorbar": "shared",
    "cbar_location": "right",
    "wspace": 0.18,
}

PLOT_SPECS = {
    # genepc1
    "genepc1":               {"cmap": "magma", "title": "PC1 fsaverage 10k", "template": "fsaverage", "density": "10k"},
    "genepc1_fsLR32k":       {"cmap": "magma", "title": "PC1 fsLR 32k", "template": "fsLR", "density":'32k'},

    # myelin
    "myelin":                {"cmap": "nipy_spectral", "vmin": 0.98, "vmax": 1.9, "title": "T1w/T2w ratio fsLR 32k"},
    "myelinmap_fsLR32k":     {"cmap": "nipy_spectral", "vmin": 0.98, "vmax": 1.9, "title": "T1w/T2w ratio fsLR 32k"},

    # developmental expansion (R only)
    "devexp":                {"hemi": "R", "cmap": "blue_orange", "vmin": -0.61, "vmax": 0.61, "title": "Developmental expansion fsLR 164k", "template": "fsLR", "density": "164k"},
    "devexp_fsLR32k":        {"hemi": "R", "cmap": "blue_orange", "vmin": -0.6,  "vmax": 0.59, "title": "Developmental expansion fsLR 32k", "template": "fsLR", "density": "32k"},

    # evolutionary expansion (R only)
    "evoexp":                {"hemi": "R", "cmap": "blue_orange", "vmin": -2.7, "vmax": 2.3, "title": "Evolutionary expansion fsLR 164k", "density": "164k"},
    "evoexp_fsLR32k":        {"hemi": "R", "cmap": "blue_orange", "vmin": -2.7, "vmax": 2.3, "title": "Evolutionary expansion fsLR 32k"},

    # functional gradient
    "gradient_pc1_fsLR32k":  {"cmap": "jet", "vmin": -5.4, "vmax": 6.8, "title": "Functional gradient fsLR 32k"},

    # intersubject variability
    "isv":                   {"cmap": "blue_orange", "vmin": 0.53, "vmax": 0.79, "title": "Intersubject variability fsLR 164k", "template": "fsLR", "density": "164k"},
    "isv_fsLR32k":           {"cmap": "blue_orange", "vmin": 0.53, "vmax": 0.79, "title": "Intersubject variability fsLR 32k", "template": "fsLR", "density": "32k"},

    # physiology
    "cbf":                   {"cmap": "nipy_spectral", "vmin": 4600, "vmax": 7000, "title": "Cerebral blood flow fsLR 164k", "template": "fsLR", "density": "164k"},
    "cbf_fsLR32k":           {"cmap": "nipy_spectral", "vmin": 4600, "vmax": 7000, "title": "Cerebral blood flow fsLR 32k", "template": "fsLR", "density": "32k"},

    "cbv":                   {"cmap": "nipy_spectral", "vmin": 0, "vmax": 26000, "title": "Cerebral blood volume fsLR 164k", "template": "fsLR", "density": "164k"},
    "cbv_fsLR32k":           {"cmap": "nipy_spectral", "vmin": 2100, "vmax": 13000, "title": "Cerebral blood volume fsLR 32k", "template": "fsLR", "density": "32k"},

    "cmro2":                 {"cmap": "nipy_spectral", "vmin": 4000, "vmax": 7500, "title": "Oxygen metabolism fsLR 164k", "density": "164k"},
    "cmro2_fsLR32k":         {"cmap": "nipy_spectral", "vmin": 4000, "vmax": 7500, "title": "Oxygen metabolism fsLR 32k", "template": "fsLR", "density": "32k"},

    "cmrglc":                {"cmap": "nipy_spectral", "vmin": 360, "vmax": 8500, "title": "Glucose metabolism fsLR 164k", "density": "164k"},
    "cmrglc_fsLR32k":        {"cmap": "nipy_spectral", "vmin": 360, "vmax": 8500, "title": "Glucose metabolism fsLR 32k", "template": "fsLR", "density": "32k"},

    # allometric scaling?
    "scalingnih":            {"cmap": "seismic", "vmin": 0.45, "vmax": 1.6, "title": "Allometric scaling (NIH) civet 41k", "template": "civet", "density": "41k"},
    "scalingnih_fsLR32k":    {"cmap": "seismic", "vmin": 0.0,  "vmax": 1.6, "title": "Allometric scaling (NIH) fsLR 32k", "template": "fsLR", "density": "32k"},

    "scalingpnc":            {"cmap": "seismic", "vmin": 0.4, "vmax": 1.7, "title": "Allometric scaling (PNC) civet 41k", "template": "civet", "density": "41k"},
    "scalingpnc_fsLR32k":    {"cmap": "seismic", "vmin": 0.0, "vmax": 1.7, "title": "Allometric scaling (PNC) fsLR 32k", "template": "fsLR", "density": "32k"},

    # source
    "source_thickness":      {"title": "Source map"},
}

def get_plot_kwargs(name: str) -> dict:
    """
    Merge defaults with per-map overrides.
    Tries exact match, then common suffix/prefix variants.
    """
    out = dict(PLOT_DEFAULTS)
    spec = PLOT_SPECS.get(name)

    if spec is None:
        candidates = [
            name.replace("_fsLR32k", ""),        # strip suffix
            f"{name}_fsLR32k",                   # add suffix
        ]
        for cand in candidates:
            if cand in PLOT_SPECS:
                spec = PLOT_SPECS[cand]
                break

    if spec:
        out.update(spec)
    return out

def get_source():
    # hcps1200 thickness fsLR-32k
    return fetch_annotation(source="hcps1200", desc="thickness", space="fsLR", den="32k")

def get_targets():
    return dict(
        genepc1 = dict(source="abagen", desc="genepc1",  space="fsaverage", den="10k"),
        myelin  = dict(source="hcps1200", desc="myelinmap", space="fsLR", den="32k"),
        devexp  = dict(source="hill2010", desc="devexp",    space="fsLR", den="164k", hemi="R"),
        evoexp  = dict(source="hill2010", desc="evoexp",    space="fsLR", den="164k", hemi="R"),
        gradient_pc1 = dict(source="margulies2016", desc="fcgradient01", space="fsLR", den="32k"),
        isv     = dict(source="mueller2013", desc="intersubjvar", space="fsLR", den="164k"),
        cbf     = dict(source="raichle", desc="cbf",   space="fsLR", den="164k"),
        cbv     = dict(source="raichle", desc="cbv",   space="fsLR", den="164k"),
        cmro2   = dict(source="raichle", desc="cmr02", space="fsLR", den="164k"),
        cmrglc  = dict(source="raichle", desc="cmrglc",space="fsLR", den="164k"),
        scalingnih = dict(source="reardon2018", desc="scalingnih", space="civet", den="41k"),
        scalingpnc = dict(source="reardon2018", desc="scalingpnc", space="civet", den="41k"),
    )

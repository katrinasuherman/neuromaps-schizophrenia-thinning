# neuromaps-schizophrenia-thinning

This repository reproduces the analytic workflow of Figure 4a from Markello & Mišić (2021) using the neuromaps Python library.
It compares a schizophrenia‐related cortical‐thickness map to multiple structural / functional / genetic brain maps by transforming surfaces, correlating them, applying spatial spin tests, correcting via FDR, and generating both surface figures and null-distribution boxplots.

## Repository Layout
.
├─ run.py                       # command-line entrypoint
├─ configs/config.json          # global parameters
├─ requirements.txt
└─ src/brainmaps/
   ├─ config.py                 # Config class
   ├─ env.py                    # environment check for Workbench
   ├─ helpers.py, catalog.py    # helper functions
   ├─ transforms.py             # brain map conversions
   ├─ stats.py                  # spin tests + FDR
   ├─ plotting.py               # surface plotting
   ├─ pipeline.py               # orchestrates the steps
   └─ boxplot.py                # generates p_spin-based boxplot


## Output Structure (out/)
out/
├─ correlations.csv             # correlation and p_spin results
├─ correlations_fdr.csv         # correlation, p_spin, and FDR
├─ figs/
│   ├─ <map>.png                # individual surface plots
│   └─ boxplots.png             # result
└─ nulls/
    └─ <map>.npy                # 1-D null correlation arrays

## Setup
### clone and enter
git clone https://github.com/<yourname>/neuromaps-schizophrenia-thinning.git
cd neuromaps-schizophrenia-thinning

### create virtualenv
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt


On macOS, install Connectome Workbench and ensure /Applications/Workbench/bin_macosx64 is in your PATH.

## Configuration
Edit configs/config.json:
- `out_dir` – output directory
- `seed` – random seed for reproducibility
- `n_perm` – number of spin-test permutations

## Run the Pipeline
| Step                      | Command                    | Description                                         |
| ------------------------- | -------------------------- | --------------------------------------------------- |
| Environment check         | `python run.py env`        | Verify neuromaps + Workbench availability           |
| Transform maps            | `python run.py transforms` | Fetch & transform all target maps                   |
| Correlations + spin tests | `python run.py stats`      | Saves `correlations.csv` and 1-D nulls              |
| FDR correction            | `python run.py fdr`        | Adds `p_fdr` column                                 |
| Surface plots             | `python run.py viz`        | Generates per-map figures                           |
| Boxplot                   | `python run.py results`    | Creates `figs/boxplots.png` using *p_spin*          |
| Full pipeline             | `python run.py all`        | Runs env → transforms → stats → fdr → viz → results |

You can delete old outputs anytime:
rm -rf out/

### Interpretation
- Red dots in `boxplots.png` = maps where `p_spin < 0.05` (significant correlation).
- Gray boxes = null distributions from spatial permutations.
- `correlations.csv` lists Pearson r and p_spin per target map.
- `correlations_fdr.csv` additionally provides p_fdr values for FDR-controlled significance.

### Citation
Markello, R. D., & Mišić, B. (2021). Comparing spatially resolved brain maps using the neuromaps toolbox. Nature Communications, 12(1), 1–13.
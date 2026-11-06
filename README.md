# Adiabatic_CEE_solver

Code and data accompanying **“The Common Envelope Evolution Outcome. III. The Improvement of
Stellar Binding Energy with the Envelope Residual”**.

You can use python to call the RBF interpolation table, or use the formulae that we fit.

Given the initial mass Mi, the initial radius Ri, the residual envelope mass Menv and the
evolutionary stage of a donor star, this repository returns

- the binding-energy parameter of the adiabatic model, λ_ad (and its building blocks
  λ_core, F_res, F_ad), and
- the final stellar radius after the common-envelope phase, R_f.

Two independent interfaces are provided:

| Interface | Location | Method | Accuracy | Speed |
|---|---|---|---|---|
| RBF interpolation | `Adiab_CEE/RBF/` | Radial-basis-function tables on a 3-D grid | Grid data, highest fidelity | Very fast |
| Polynomial fits | `Adiab_CEE/formulae/` | Analytic fitting formulae (3-D and 2-D) | Fitting accuracy (Sect. 4 of the paper) | Instantaneous, no data files |

Notation (paper and code):

λ_ad = λ_core × F_res × F_ad  ⇔  `lam_ad = lam0 * F1 * F2`

where `F1` ≡ F_res (envelope-residual correction) and `F2` ≡ F_ad (adiabatic-expansion correction).

---

## Repository layout

```
Adiabatic_CEE_solver/
├── README.md
└── Adiab_CEE/
    ├── RBF/                       # RBF interpolation tables (recommended)
    │   ├── Adiab_CEE.py           #   <- key module
    │   ├── test.py                #   <- minimal working example
    │   ├── formula_test.py        #   cross-check against the formulae
    │   ├── Readme.txt             #   short notes (legacy)
    │   └── data/                  #   .npy grids (~500 MB in total)
    │       ├── loglam0_grid_<stage>_<core>.npy
    │       ├── F1_grid_<stage>_<core>_<scheme>.npy
    │       ├── F2_grid_<stage>_<core>_<scheme>.npy
    │       └── Rf_grid_<stage>_<core>.npy
    └── formulae/                  # polynomial fitting formulae
        ├── formula_test.py        #   <- key module
        ├── test_poly.py           #   <- minimal working example
        ├── Readme.txt             #   short notes (legacy)
        └── lambda_CE_formula/     #   .csv coefficient tables
```

`.ipynb_checkpoints/` folders are leftovers from Jupyter and can be deleted.

---

## Requirements

- Python >= 3.8
- NumPy (interpolation and file I/O)
- pandas (only for reading the `.csv` coefficient tables in `formulae/`)

```bash
pip install numpy pandas
```

---

## Using the RBF tables

`Adiab_CEE/RBF/Adiab_CEE.py` loads the `.npy` grids through the **relative** path `data/...`, so run
your script **from inside `Adiab_CEE/RBF/`** (or copy `data/` next to your script).

```python
from Adiab_CEE import RecCEE_solver, Rfad_solver

Mi, Ri, Menv = 1.4, 120.0, 1e-1     # Msun, Rsun, Msun
stage       = 'RGB'                 # 'HG' | 'RGB' | 'AGB' | 'massive'
core        = 'core1'               # 'core1' (X_H = 0.01) | 'core2' (X_H = 0.1)
use_grid    = 'abs'                 # 'abs' (Menv = eta * 1 Msun) | 'rel' (Menv = eta * Mi)

lam0, F1, F2 = RecCEE_solver(Mi, Ri, Menv, stage, core, use_grid)
lam_ad       = lam0 * F1 * F2

Rf = Rfad_solver(Mi, Ri, Menv, stage, core)   # Rsun
```

### `RecCEE_solver(Mi, Ri, Menv, stage, core, use_grid)`

| Argument | Unit / value | Meaning |
|---|---|---|
| `Mi` | Msun | initial mass of the donor before the CE |
| `Ri` | Rsun | initial radius before the CE |
| `Menv` | Msun | residual envelope mass after the CE, 1e-4 <= Menv <= 0.1 |
| `stage` | str | `'HG'`, `'RGB'`, `'AGB'` or `'massive'` |
| `core` | str | `'core1'` (X_H = 0.01) or `'core2'` (X_H = 0.1) |
| `use_grid` | str | `'abs'`: tables built with Menv = eta x 1 Msun (Sets 1, 2); `'rel'`: Menv = eta x Mi (Sets 3, 4) |

Returns `[lam0, F1, F2]` = [λ_core, F_res, F_ad]; multiply the three to obtain λ_ad. `F1` and `F2`
are floored at 1 (they are correction factors and are never below unity in the tables).

### `Rfad_solver(Mi, Ri, Menv, stage, core)`

Returns the final radius R_f in Rsun after the CE. Internally the tables are indexed by
eta = Menv/Mi with 0 <= eta <= 0.1, so **`Menv` must satisfy `Menv < 0.1 * Mi`**.

### Tabulated ranges (inputs must lie inside them)

| `stage` | Mi (Msun) | Ri (Rsun) |
|---|---|---|
| `HG` | 1.2 - 110 | 2 - 2500 |
| `RGB` | 0.2 - 21 | 0.9 - 1200 |
| `AGB` | 0.5 - 21 | 15 - 1300 |
| `massive` | 19 - 110 | 15 - 3300 |

The residual envelope mass axis is common to all stages: 1e-4 to 0.1 Msun (31 points in log space).
The interpolation is tri-linear in (log Mi, log Ri, log Menv).

### File-naming convention

```
loglam0_grid_<stage>_<core>.npy          # log10(lambda_core) on the (Mi, Ri) grid, 100 x 200
F1_grid_<stage>_<core>_<scheme>.npy      # F_res, 100 x 200 x 31
F2_grid_<stage>_<core>_<scheme>.npy      # F_ad,  100 x 200 x 31
Rf_grid_<stage>_<core>.npy               # R_f,   100 x 200 x 11  (eta axis)
```

with `<core>` in {`0.01`, `0.1`} and `<scheme>` in {`abs`, `rel`}.

---

## Using the polynomial formulae

```python
from formula_test import lambda_ce_3d_fgb, lambda_ce_fgb_eta_1e2

path = 'lambda_CE_formula'                 # folder holding the .csv coefficients

lam_3d = lambda_ce_3d_fgb(1.5, 60.0, 0.0056, path)   # (Mi, Ri, eta, path)
lam_2d = lambda_ce_fgb_eta_1e2(1.5, 60.0)            # (Mi, Ri) at eta = 1e-2
```

Available functions (all return a single λ_ad, which reproduces λ_core F_res F_ad from the RBF
tables to fitting accuracy):

| 3-D (Mi, Ri, eta) | 2-D (fixed eta) |
|---|---|
| `lambda_ce_3d_hg(Mi, Ri, eta, path)` | `lambda_ce_hg_eta_1e2(Mi, Ri)`, `lambda_ce_hg_eta_1e3(Mi, Ri)` |
| `lambda_ce_3d_fgb(Mi, Ri, eta, path)` | `lambda_ce_fgb_eta_1e2(Mi, Ri)`, `lambda_ce_fgb_eta_1e3(Mi, Ri)` |
| `lambda_ce_3d_agb(Mi, Ri, eta, path)` | `lambda_ce_agb_eta_1e2(Mi, Ri)`, `lambda_ce_agb_eta_1e3(Mi, Ri)` |
| `lambda_ce_3d_massive(Mi, Ri, eta, path)` | `lambda_ce_massive_eta_1e2(Mi, Ri)`, `lambda_ce_massive_eta_1e3(Mi, Ri)` |

Here eta is the **relative** residual envelope mass, eta = Menv/Mi (0 <= eta <= 0.1); `_1e2` and
`_1e3` denote eta = 1e-2 and eta = 1e-3. The coefficient files in `lambda_CE_formula/` must be
present; run the script from `Adiab_CEE/formulae/`.

---

## Notes and pitfalls

- **Working directory matters.** Both modules read their data with relative paths; run scripts inside
  the corresponding folder (`RBF/` or `formulae/`) or pass the correct `folder_path` for the CSV files.
- **Stay inside the grid.** `RecCEE_solver` performs a bracketing search on the log-grids; inputs
  outside the tabulated Mi, Ri or Menv ranges (see table above) are not extrapolated and will fail.
  Use the polynomial formulae for exploratory work outside the grid.
- The edge value `Menv = 0.1` is internally shifted by 1e-10 to stay inside the last cell.
- `Rfad_solver` requires `Menv < 0.1 * Mi` because it interpolates in eta over 0 to 0.1.
- **Which `core`/`use_grid` to choose:** `core` sets the core-envelope boundary definition used to
  build the grid (X_H = 0.01 or 0.1); `use_grid` selects the family of tables built with the two
  prescriptions for the residual mass (Menv = eta x 1 Msun for lower-mass stars, and Menv = eta x Mi
  for higher-mass stars; see Sect. 3 of the paper).
- The RBF tables are data: they inherit the resolution of the underlying stellar grids. Differences
  between `core1` and `core2` are <~15% for most stars (<~50% for 3-10 Msun RGB/EAGB stars), see
  Fig. 4 and Sect. 3.1 of the paper.
- `formula_test.py` (in `RBF/`) can be used to verify that λ_core F_res F_ad and the polynomial fits
  agree.

---

## Reproducing the figures/tables of the paper

- Figs. 5-7: evaluate `RecCEE_solver` (or the fits) over the (Mi, Ri) grid with the four parameter
  sets and compute the relative deviations against the raw model values.
- Table 2: the last column is the error of λ_core F_res F_ad relative to the raw λ_ad, i.e.
  |(λ_core F_res F_ad)_pred - (λ_ad)_raw| / (λ_ad)_raw — *not* the product of the three individual
  relative errors.

---

## Data provenance

The grids were computed from adiabatic mass-loss sequences of the Eggleton stellar-evolution code for
the four model sets described in the paper (two residual-mass prescriptions x two core-boundary
definitions). The RBF tables and the polynomial coefficients are two representations of the same
underlying models.

---

## Citation

The article has been submitted to ApJ. If you use this code or these tables, please cite the paper:

```bibtex
@ARTICLE{2026arXiv260924010Z,
       author = {{Zhang}, Lifu and {Ge}, Hongwei and {Hebrail}, Dylan and {Church}, Ross and
                 {Yang}, Mingkuan and {Li}, Zhenwei and {Chen}, Hailiang and {Jiang}, Dengkai and
                 {Chen}, Xuefei and {Han}, Zhanwen},
        title = "{The Common Envelope Evolution Outcome. III. the Improvement of Stellar Binding
                  Energy with the Envelope Residual}",
      journal = {arXiv e-prints},
     keywords = {Solar and Stellar Astrophysics},
         year = 2026,
        month = sep,
          eid = {arXiv:2609.24010},
        pages = {arXiv:2609.24010},
          doi = {10.48550/arXiv.2609.24010},
archivePrefix = {arXiv},
       eprint = {2609.24010},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2026arXiv260924010Z},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
Please also cite the archived version of this repository (Zenodo release):

```bibtex
@misc{<zenodo-key>,
  author    = {Zhang, Lifu and Ge, Hongwei and ...},
  title     = {{The data for The Common Envelope Evolution Outcome. III. The Improvement of
               Stellar Binding Energy with the Envelope Residual}},
  year      = 2026,
  publisher = {Zenodo},
  version   = {1.0.0},
  doi       = {10.5281/zenodo.24014183},
  url       = {https://doi.org/10.5281/zenodo.24014182}
}
```
---

## License

The code and data in this repository are released under the
**Creative Commons Zero v1.0 Universal (CC0 1.0)** public-domain dedication - see [`LICENSE`](LICENSE).
You may copy, modify and redistribute the material for any purpose without asking permission.

---

## Contact

Hongwei Ge (Yunnan Observatories, CAS) - for questions about the models and the data.
`Adiab_CEE.py` was written by Lifu Zhang (v1.0.0).


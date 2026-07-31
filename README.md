# Weak Factor Analysis via Temporally Weighted PCA

This folder reorganizes the original project into a single replication package for the JRSSB revision.

The package follows two rules:

- the historical source files and data are retained under `original/`
- the editor-facing scripts are organized by paper output under `scripts/`

This audit copy makes only the corrections needed to align the executable code with
the manuscript: symmetric WPCA uses the largest algebraic eigenvalues, the simulation
factor transition is \(0.9I_r\), the loading matrix is \(1.2U_sD_L\), and reported
projector errors are raw Frobenius norms.

The inference-figure normalization uses the theorem-consistent right-alignment
convention and rotation-compatible whitening factors. For orthonormal \(A,B\),
the returned matrix \(R\) satisfies \(AR\approx B\); the loading and factor
whiteners satisfy \(W\Sigma W^\top=I_r\).

The pre-revision numerical results are organized from the original code base:

- `simulation_est.R`
- `simulation_gamma.R`
- `simulation_inf_F.R`
- `simulation_inf_L.R`
- `fred.R`
- `functions.R`

The two revision-only additions are kept separate:

- `Figure2.R` uses the minimally corrected revision simulation script for the line plot
- `Figure7.R` uses the shared corrected runtime helpers together with a small interpretive wrapper for the new real-data figure (Figure 3 in the manuscript)

## Folder Layout

- `original/codes/`: untouched copies of the original core scripts and FRED data files
- `original/revision_figure2/`: revision Figure 2 script with the same minimal estimator, DGP, and metric corrections
- `original/revision_figure7/`: untouched copy of the revision Figure 7 script
- `R/`: shared helper modules used by the editor-facing scripts
  - `R/functions.R` is the runtime copy used by the organized scripts
  - `original/codes/functions.R` remains the untouched original reference copy
- `scripts/`: one script per paper figure/table
- `outputs/figs/`: generated figures
- `outputs/tables/`: generated table files and supporting CSV files

## Output Mapping

### Figures

- Manuscript Figure 1: `scripts/Figure1.R`
  - Generates `outputs/figs/sim_inf_L_N100_T800_hist.png`
  - Generates `outputs/figs/sim_inf_L_PCA_N100_T800_hist.png`
  - Generates `outputs/tables/Figure1_loading_error_summary.csv`
  - Organized from original `simulation_inf_L.R`

- Manuscript Figure 2: `scripts/Figure2.R`
  - Generates `outputs/figs/fig_error_lines.png`
  - Generates supporting CSV files in `outputs/tables/`
  - Runs the corrected revision script `original/revision_figure2/run_uv_v4_spd_noise.R`
  - Default official run uses `n_rep = 300` and `seed = 20260315`, matching the manuscript revision figure

- Manuscript Figure 3: `scripts/Figure7.R`
  - Generates `outputs/figs/fig_interpretive_addon.png`
  - Generates supporting CSV files in `outputs/tables/`
  - Uses the shared corrected runtime helpers, original FRED data files, and an interpretive wrapper
  - Uses 10 expanding windows: 24-month steps with minimum 240 observations for FRED-MD,
    and 8-quarter steps with minimum 80 observations for FRED-QD

- Manuscript Figure 4: `scripts/Figure3.R`
  - Generates `outputs/figs/sim_inf_L_N200_T200.png`
  - Generates `outputs/figs/sim_inf_L_N200_T200_hist.png`
  - Generates `outputs/tables/Figure3_loading_error_summary.csv`
  - Organized from original `simulation_inf_L.R`

- Manuscript Figure 5: `scripts/Figure4.R`
  - Generates `outputs/figs/sim_inf_L_N100_T800.png`
  - Generates `outputs/figs/sim_inf_L_N100_T800_hist.png`
  - Generates `outputs/tables/Figure4_loading_error_summary.csv`
  - Organized from original `simulation_inf_L.R`

- Manuscript Figure 6: `scripts/Figure5.R`
  - Generates `outputs/figs/sim_inf_N100_T100.png`
  - Generates `outputs/figs/sim_inf_N100_T100_hist.png`
  - Generates `outputs/tables/Figure5_factor_error_summary.csv`
  - Organized from original `simulation_inf_F.R`

- Manuscript Figure 7: `scripts/Figure6.R`
  - Generates `outputs/figs/sim_inf_N100_T500.png`
  - Generates `outputs/figs/sim_inf_N100_T500_hist.png`
  - Generates `outputs/tables/Figure6_factor_error_summary.csv`
  - Organized from original `simulation_inf_F.R`

### Tables

- Manuscript Table 2 (`tab:fred-md`): `scripts/Table5.R`
  - Generates `outputs/tables/Table5_fred_reconstruction.csv`
  - Organized from original `fred.R`
  - Recomputes by default; `cached` mode is retained only for comparison with the historical output

- Manuscript Table 3 (`tab:sim_est_off_diag`): `scripts/Table2.R`
  - Generates `outputs/tables/Table2_sim_est_off_diag.csv`
  - Organized from original `simulation_est.R`
  - Off-diagonal covariance setting

- Manuscript Table 4 (`tab:sim_est_diag`): `scripts/Table3.R`
  - Generates `outputs/tables/Table3_sim_est_diag.csv`
  - Organized from original `simulation_est.R`
  - Diagonal covariance setting
  - Recomputes by default; `cached` mode is retained only for comparison with the historical output

- Manuscript Table 5 (`tab:sim_cv`): `scripts/Table4.R`
  - Generates `outputs/tables/Table4_sim_cv.csv`
  - Organized from original `simulation_gamma.R`
  - Recomputes by default; `cached` mode is retained only for comparison with the historical output

Table 1 in the manuscript is a conceptual comparison table rather than a numerical replication result, so it is not generated by code here.

### Diagnostic

- `diagnostics/check_inference_normalization.R`
  - Checks the compatible-whitening identities numerically
  - Reports coordinate-wise means, medians, and standard deviations for the
    theorem-normalized loading and factor statistics
  - Optional first argument: number of replications

## How To Run

Run any single script from this folder, for example:

```r
Rscript scripts/Table2.R
Rscript scripts/Figure7.R  # manuscript Figure 3
```

To run the full paper package:

```r
Rscript scripts/run_all_paper_outputs.R
```

This default is the formal `paper` mode. It runs Figure 2 with the manuscript settings and recomputes Tables 2 to 5 using the corrected estimator.

For a faster smoke test:

```r
Rscript scripts/run_all_paper_outputs.R quick
```

The `quick` mode is only for pipeline checking. It overwrites the same files in `outputs/`, so rerun the default `paper` mode before preparing the final repository snapshot.

## Optional Arguments

Most scripts accept optional command-line arguments:

- simulation tables and inference figures: first argument `n_rep`, second argument `seed`
- `Figure2.R`: first argument `n_rep`, second argument `seed`
- `Figure7.R` (manuscript Figure 3): first argument `seed`
- `Table3.R` and `Table4.R`: default mode is `recompute`; use `cached` only to inspect the historical outputs
- `Table5.R`: default mode is `recompute`; use `cached` only to inspect the historical output

If `fig_error_lines.png` looks less smooth than the manuscript version, rerun `Figure2.R` without arguments or with `300 20260315`. Small-replicate test runs overwrite the same output file and can make the trend visibly noisier.

If no arguments are supplied, each script uses the manuscript-style default values.

## Official Run Settings

| Outputs | Replications | Seed |
|---|---:|---:|
| Manuscript Figures 1, 4, and 5 | 500 | 123 |
| Figure 2 | 300 | 20260315 |
| Manuscript Figures 6 and 7 | 500 | 1234 |
| Manuscript Figure 3 | deterministic data analysis | 20260315 |
| Tables 2 to 5 | 100 | 1234 |

## Software Notes

The paper-output scripts use `MASS`, `ggplot2`, and base R. The retained legacy
`fredqd()` helper additionally requires `readr` if it is called directly.

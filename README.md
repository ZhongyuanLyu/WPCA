# Replication Code for Weak Factor Analysis via Temporally Weighted PCA

This public repository contains the code and data needed to reproduce the numerical
results in the paper *Weak Factor Analysis via Temporally Weighted PCA*. It provides
one executable entry point for each paper figure and numerical table, together with
a run-all script and a faster pipeline check.

The main replication scripts are organized by paper output under
<code>scripts/</code>. Shared implementation functions are under <code>R/</code>,
the source scripts and data from which the package was organized are retained under
<code>original/</code>, and generated files are written to <code>outputs/</code>.

## Quick Start

To reproduce all numerical paper outputs from the repository root:

    Rscript scripts/run_all_paper_outputs.R

This default <code>paper</code> mode uses the paper settings listed below.

For a faster check that the complete pipeline runs:

    Rscript scripts/run_all_paper_outputs.R quick

The <code>quick</code> mode uses fewer simulation replications and overwrites the
same files in <code>outputs/</code>. Rerun the default <code>paper</code> mode
before using the generated files as the paper replication outputs.

## Implementation Details

The implementation used for the paper has the following conventions:

- symmetric WPCA uses the eigenvectors associated with the largest algebraic
  eigenvalues;
- the simulation factor transition is \(0.9I_r\);
- the loading matrix is \(1.2U_sD_L\);
- projector errors are reported as raw Frobenius norms; and
- the inference figures use right alignment and rotation-compatible whitening.

For orthonormal \(A,B\), the alignment matrix \(R\) returned by the helper functions
satisfies \(AR\approx B\). The loading and factor whiteners satisfy
\(W\Sigma W^\top=I_r\).

## Folder Layout

- <code>scripts/</code>: one executable script per paper figure or numerical table
- <code>R/</code>: shared runtime functions used by the replication scripts
- <code>diagnostics/</code>: numerical checks for inference normalization
- <code>original/codes/</code>: source scripts and included FRED data files
- <code>original/revision_figure2/</code>: source implementation used for manuscript Figure 2
- <code>original/revision_figure7/</code>: source implementation underlying manuscript Figure 3
- <code>outputs/figs/</code>: generated figures
- <code>outputs/tables/</code>: generated tables and supporting CSV files

The executable scripts use <code>R/functions.R</code>.
<code>original/codes/functions.R</code> is retained as the corresponding source copy.

## Output Mapping

### Figures

**Manuscript Figure 1 — <code>scripts/Figure1.R</code>**

- Generates <code>outputs/figs/sim_inf_L_N100_T800_hist.png</code>
- Generates <code>outputs/figs/sim_inf_L_PCA_N100_T800_hist.png</code>
- Generates <code>outputs/tables/Figure1_loading_error_summary.csv</code>

**Manuscript Figure 2 — <code>scripts/Figure2.R</code>**

- Generates <code>outputs/figs/fig_error_lines.png</code>
- Generates supporting CSV files in <code>outputs/tables/</code>
- Runs <code>original/revision_figure2/run_uv_v4_spd_noise.R</code>
- The paper run uses <code>n_rep = 300</code> and <code>seed = 20260315</code>

**Manuscript Figure 3 — <code>scripts/Figure7.R</code>**

- Generates <code>outputs/figs/fig_interpretive_addon.png</code>
- Generates supporting CSV files in <code>outputs/tables/</code>
- Uses the shared runtime helpers, included FRED data files, and the paper's
  data-processing workflow
- Uses 10 expanding windows: 24-observation increments with a minimum allowed
  training size of 240 observations for FRED-MD, and 8-observation increments
  with a minimum allowed training size of 80 observations for FRED-QD

**Manuscript Figure 4 — <code>scripts/Figure3.R</code>**

- Generates <code>outputs/figs/sim_inf_L_N200_T200.png</code>
- Generates <code>outputs/figs/sim_inf_L_N200_T200_hist.png</code>
- Generates <code>outputs/tables/Figure3_loading_error_summary.csv</code>

**Manuscript Figure 5 — <code>scripts/Figure4.R</code>**

- Generates <code>outputs/figs/sim_inf_L_N100_T800.png</code>
- Generates <code>outputs/figs/sim_inf_L_N100_T800_hist.png</code>
- Generates <code>outputs/tables/Figure4_loading_error_summary.csv</code>

**Manuscript Figure 6 — <code>scripts/Figure5.R</code>**

- Generates <code>outputs/figs/sim_inf_N100_T100.png</code>
- Generates <code>outputs/figs/sim_inf_N100_T100_hist.png</code>
- Generates <code>outputs/tables/Figure5_factor_error_summary.csv</code>

**Manuscript Figure 7 — <code>scripts/Figure6.R</code>**

- Generates <code>outputs/figs/sim_inf_N100_T500.png</code>
- Generates <code>outputs/figs/sim_inf_N100_T500_hist.png</code>
- Generates <code>outputs/tables/Figure6_factor_error_summary.csv</code>

### Tables

**Manuscript Table 2 — <code>scripts/Table5.R</code>**

- Generates <code>outputs/tables/Table5_fred_reconstruction.csv</code>
- Recomputes by default; <code>cached</code> mode reads the included reference output

**Manuscript Table 3 — <code>scripts/Table2.R</code>**

- Generates <code>outputs/tables/Table2_sim_est_off_diag.csv</code>
- Uses the off-diagonal covariance setting

**Manuscript Table 4 — <code>scripts/Table3.R</code>**

- Generates <code>outputs/tables/Table3_sim_est_diag.csv</code>
- Uses the diagonal covariance setting
- Recomputes by default; <code>cached</code> mode reads the included reference output

**Manuscript Table 5 — <code>scripts/Table4.R</code>**

- Generates <code>outputs/tables/Table4_sim_cv.csv</code>
- Recomputes by default; <code>cached</code> mode reads the included reference output

Table 1 is a conceptual comparison table rather than a numerical result, so it is
not generated by code.

## Running Individual Outputs

Run any single script from the repository root, for example:

    Rscript scripts/Table2.R
    Rscript scripts/Figure7.R

Most scripts accept optional command-line arguments:

- <code>Figure1.R</code> through <code>Figure6.R</code>: first argument
  <code>n_rep</code>, second argument <code>seed</code>
- <code>Figure7.R</code>: first argument <code>seed</code>
- <code>Table2.R</code>: first argument <code>n_rep</code>, second argument
  <code>seed</code>
- <code>Table3.R</code>, <code>Table4.R</code>, and <code>Table5.R</code>:
  first argument <code>recompute</code> or <code>cached</code>, followed in
  <code>recompute</code> mode by <code>n_rep</code> and <code>seed</code>. The
  default mode is <code>recompute</code>; <code>cached</code> reads the included
  reference output.

If <code>fig_error_lines.png</code> looks less smooth after a small-replication
pipeline check, rerun:

    Rscript scripts/Figure2.R 300 20260315

If no arguments are supplied, each script uses the paper run settings.

## Paper Run Settings

| Outputs | Replications | Seed |
|---|---:|---:|
| Manuscript Figures 1, 4, and 5 | 500 | 123 |
| Manuscript Figure 2 | 300 | 20260315 |
| Manuscript Figures 6 and 7 | 500 | 1234 |
| Manuscript Figure 3 | randomized CV with (J=10) over 10 expanding windows | 20260315 |
| Manuscript Tables 2 to 5 | 100 | 1234 |

## Diagnostic

Run:

    Rscript diagnostics/check_inference_normalization.R

This checks the compatible-whitening identities numerically and reports
coordinate-wise means, medians, and standard deviations for the normalized
loading and factor statistics. An optional first argument sets the number of
replications.

## Software

The replication scripts use <code>MASS</code>, <code>ggplot2</code>, and base R.
The retained <code>fredqd()</code> source helper additionally requires
<code>readr</code> if it is called directly.

## Reproducibility Checks

[REPRODUCIBILITY_CHECKS.md](REPRODUCIBILITY_CHECKS.md) records the implementation
checks, output hashes, and paper-to-output comparisons used to validate the
generated figures and tables.

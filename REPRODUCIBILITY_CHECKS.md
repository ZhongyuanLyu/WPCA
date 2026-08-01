# Reproducibility checks

The checks below use the default settings in `scripts/run_all_paper_outputs.R`. The checksums do not apply to quick mode.

## Run settings

- Figures 1, 4, and 5: 500 replications, seed 123
- Figure 2: 300 replications, seed 20260315
- Figures 6 and 7: 500 replications, seed 1234
- Figure 3: randomized CV with seed 20260315
- Tables 2--5: 100 replications, seed 1234

## Table checksums

These reference files were generated with R 4.5.3 on aarch64-apple-darwin20, MASS 7.3.65, and ggplot2 4.0.3. Other R, BLAS, or package versions may not produce byte-identical files. The line counts include the header.

| Paper item | File in `outputs/tables/` | Lines | SHA-256 |
|---|---|---:|---|
| Table 2 | `Table5_fred_reconstruction.csv` | 55 | `cfc0b1487549d7ec6cbbde3aa312ca4e72243a1a66018b15256287e47d5695da` |
| Table 3 | `Table2_sim_est_off_diag.csv` | 17 | `acebcb49ff8478ef8cc6c75c92f5dfffbff8d9c94dddafbef172b997d7b3e2fe` |
| Table 4 | `Table3_sim_est_diag.csv` | 17 | `334bf24a024a762ea31fedfbce3df0a71eb8ba0b95f84101c37f256f543d2924` |
| Table 5 | `Table4_sim_cv.csv` | 17 | `1817aff00129f4c6d29486862878487689b2b5082298d126bf5a98d455218d7d` |

On macOS or Linux, a checksum can be verified with

```bash
shasum -a 256 outputs/tables/Table5_fred_reconstruction.csv
```

## Normalization diagnostic

The largest residual in the loading and factor whitening identities was below `1.1e-14` in the reference run.

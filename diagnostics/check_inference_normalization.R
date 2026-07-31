#!/usr/bin/env Rscript

# Diagnostic only: compare the current rotation convention with the convention
# in Theorems 7--8. Run from the replication-package root.

source(file.path("R", "functions.R"))
source(file.path("R", "simulation_core.R"))

args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args)) as.integer(args[[1]]) else 100L

check_whitening_identities <- function() {
  set.seed(17)
  N <- 60L
  Times <- 80L
  r <- 3L
  signal <- generate_signal_panel(N = N, Times = Times, r = r)
  Sigma_C <- draw_sigma_c(N = N, rho = 0.6, covariance = "offdiag")
  Sigma_T <- diag(Times)
  X <- signal$M + t(MASS::mvrnorm(Times, rep(0, N), Sigma_C))
  Q <- make_weight_matrix(Times = Times, gamma = 0, symmetric = TRUE)
  U_hat <- APCA(X, gamma = 0, r = r, symmetric = TRUE)
  svd_hat <- svd(U_hat %*% t(U_hat) %*% X / sqrt(Times), nu = r, nv = r)

  rotation_L <- compute_rotation_and_covariance_L(
    svd_hat$u, svd_hat$v, signal$Fmat, signal$L,
    signal$M / sqrt(Times), signal$M %*% Q %*% t(signal$M),
    Sigma_C, Sigma_T, Q, ind = as.integer(0.5 * N)
  )
  rotation_F <- compute_rotation_and_covariance_F(
    svd_hat$v, signal$Fmat, signal$M / sqrt(Times),
    Sigma_C, Sigma_T, ind = as.integer(0.5 * Times)
  )

  identity_error_L <- max(abs(rotation_L$W_L %*% rotation_L$Sigma_L %*%
                                t(rotation_L$W_L) - diag(r)))
  inverse_error_L <- max(abs(t(rotation_L$W_L) %*% rotation_L$W_L -
                               solve(rotation_L$Sigma_L)))
  identity_error_F <- max(abs(rotation_F$W_F %*% rotation_F$Sigma_F %*%
                                t(rotation_F$W_F) - diag(r)))
  inverse_error_F <- max(abs(t(rotation_F$W_F) %*% rotation_F$W_F -
                               solve(rotation_F$Sigma_F)))

  errors <- c(identity_L = identity_error_L, inverse_L = inverse_error_L,
              identity_F = identity_error_F, inverse_F = inverse_error_F)
  print(errors)
  stopifnot(max(errors) < 1e-8)
}

summarize_coordinates <- function(x, label) {
  out <- data.frame(
    setting = label,
    coordinate = seq_len(nrow(x)),
    mean = rowMeans(x),
    sd = apply(x, 1, sd),
    median = apply(x, 1, median)
  )
  print(out, row.names = FALSE)
}

check_whitening_identities()

set.seed(123)
summarize_coordinates(
  simulate_loading_distribution(
    N = 100L, Times = 800L, rho = 0.6, gamma = 0,
    N_rep = n_rep
  ),
  "loading_N100_T800"
)

set.seed(123)
summarize_coordinates(
  simulate_loading_distribution(
    N = 100L, Times = 800L, rho = 0.6, gamma = 1,
    N_rep = n_rep
  ),
  "loading_PCA_N100_T800"
)

set.seed(1234)
summarize_coordinates(
  simulate_factor_distribution(
    N = 100L, Times = 500L, rho = 0.6, gamma = 0,
    N_rep = n_rep
  ),
  "factor_N100_T500"
)

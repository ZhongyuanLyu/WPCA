#!/usr/bin/env Rscript

# v4 simulation package (strict SPD-noise form):
#   E = Sigma_C^{1/2} Z Sigma_T^{1/2}
# with Z i.i.d. entries.
#
# Methods:
#   - PCA
#   - WPCA (lag-1)
#
# Outputs:
#   - sim_raw_u_error.csv
#   - sim_summary_u_error.csv
#   - sim_raw_v_error.csv
#   - sim_summary_v_error.csv
#   - fig_error_lines.png

args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args) >= 1) as.integer(args[1]) else 300L
seed <- if (length(args) >= 2) as.integer(args[2]) else 20260315L
if (is.na(n_rep) || n_rep <= 0) stop("n_rep must be a positive integer")
set.seed(seed)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
  output_dir <- dirname(script_path)
} else {
  output_dir <- getwd()
}

orth_basis <- function(M) {
  qr.Q(qr(M))
}

std_t <- function(n, df = 5L) {
  rt(n, df = df) / sqrt(df / (df - 2))
}

draw_iid_matrix <- function(nr, nc, dist, df = 5L) {
  if (dist == "gaussian") {
    return(matrix(rnorm(nr * nc), nrow = nr, ncol = nc))
  }
  if (dist == "t5") {
    return(matrix(std_t(nr * nc, df = df), nrow = nr, ncol = nc))
  }
  stop("Unknown distribution")
}

# simulation_est.R style metric scaling
subspace_error_simest <- function(Uhat, Utrue, r) {
  sqrt(sum((Uhat %*% t(Uhat) - Utrue %*% t(Utrue))^2) / sqrt(r))
}

VAR1_sim <- function(B0, Phi, n, sigvare, burn_in = 50L) {
  nob <- burn_in + n
  k <- nrow(B0)
  x <- matrix(0, nrow = nob, ncol = k)
  tem <- chol(sigvare)
  for (i in 2:nob) {
    x[i, ] <- B0 + Phi %*% matrix(x[i - 1, ], nrow = k) + crossprod(tem, matrix(rnorm(k), nrow = k))
  }
  x[(burn_in + 1):nob, , drop = FALSE]
}

make_q_lag1 <- function(T) {
  if (T < 2) stop("T must be >= 2")
  Q <- matrix(0, nrow = T, ncol = T)
  idx <- 1:(T - 1)
  Q[cbind(idx, idx + 1)] <- 1 / (T - 1)
  Q
}

fit_pca_u <- function(X, r) {
  svd(X, nu = r, nv = 0)$u
}

fit_wpca_lag1_u <- function(X, r) {
  Q <- make_q_lag1(ncol(X))
  svd(X %*% Q %*% t(X), nu = r, nv = 0)$u
}

fit_v_from_u <- function(X, Uhat, r) {
  svd(Uhat %*% t(Uhat) %*% X, nu = 0, nv = r)$v
}

build_sigma_t <- function(T, type, rho_t = 0) {
  if (type == "identity") {
    return(diag(1, T))
  }
  if (type == "ar1") {
    idx <- 1:T
    SigmaT <- outer(idx, idx, function(i, j) rho_t^abs(i - j))
    return(SigmaT)
  }
  stop("Unknown sigma_t type")
}

make_rest_basis <- function(core, n_rest) {
  N <- nrow(core)
  if (n_rest <= 0) return(matrix(0, nrow = N, ncol = 0))
  for (try_id in 1:20) {
    Rraw <- matrix(rnorm(N * n_rest), nrow = N, ncol = n_rest)
    Rraw <- Rraw - core %*% (t(core) %*% Rraw)
    qrR <- qr(Rraw)
    if (qrR$rank >= n_rest) {
      return(qr.Q(qrR)[, 1:n_rest, drop = FALSE])
    }
  }
  stop("Failed to build orthogonal complement basis")
}

build_sigma_c <- function(L, cfg) {
  N <- nrow(L)
  r <- cfg$r
  q <- cfg$q_spike
  n_rest <- N - r - q
  if (n_rest < 0) stop("q_spike too large")

  U_sig <- orth_basis(L)[, 1:r, drop = FALSE]
  if (q > 0) {
    raw <- cbind(U_sig, matrix(rnorm(N * q), nrow = N, ncol = q))
    full <- orth_basis(raw)
    B_perp <- full[, (r + 1):(r + q), drop = FALSE]
  } else {
    B_perp <- matrix(0, nrow = N, ncol = 0)
  }

  core <- cbind(U_sig, B_perp)
  Rrest <- make_rest_basis(core, n_rest)
  U_full <- cbind(core, Rrest)

  eig_vals <- c(
    rep(cfg$eig_signal, r),
    rep(cfg$eig_spike, q),
    rep(cfg$eig_base, n_rest)
  )
  if (any(eig_vals <= 0)) stop("Sigma_C eigenvalues must be positive")

  SigmaC <- U_full %*% diag(eig_vals, nrow = N, ncol = N) %*% t(U_full)
  SigmaC <- (SigmaC + t(SigmaC)) / 2
  SigmaC
}

draw_noise_spd <- function(L, T, cfg) {
  N <- nrow(L)
  SigmaC <- build_sigma_c(L, cfg)
  SigmaT <- build_sigma_t(T, cfg$sigma_t_type, cfg$rho_t)

  # Strict form: E = Sigma_C^{1/2} Z Sigma_T^{1/2}
  Z <- draw_iid_matrix(N, T, cfg$z_dist, df = cfg$z_df)
  C_half <- t(chol(SigmaC))
  T_half <- chol(SigmaT)
  E <- C_half %*% Z %*% T_half
  E
}

simulate_panel <- function(cfg, T) {
  N <- cfg$N
  r <- cfg$r

  # simulation_est.R-aligned factor DGP
  B0 <- matrix(0, nrow = r)
  diagPhivec <- rep(cfg$phi_level, r)
  O1 <- svd(matrix(rnorm(r * r), nrow = r))$u
  O2 <- svd(matrix(rnorm(r * r), nrow = r))$u
  Phi <- O1 %*% diag(diagPhivec, nrow = r) %*% t(O2)
  sigvare <- diag(1, r, r) - Phi %*% t(Phi)
  sigvare <- (sigvare + t(sigvare)) / 2
  min_ev <- min(eigen(sigvare, symmetric = TRUE, only.values = TRUE)$values)
  if (min_ev <= 1e-8) {
    sigvare <- sigvare + diag(abs(min_ev) + 1e-6, r)
  }
  Fmat <- VAR1_sim(B0, Phi, T, sigvare, burn_in = 50L)

  # simulation_est.R-aligned loading DGP
  svd_L <- svd(matrix(rnorm(N * r), nrow = N))
  L <- svd_L$u %*% diag(svd_L$d, nrow = r, ncol = r)

  M <- cfg$signal_scale * (L %*% t(Fmat))
  E <- draw_noise_spd(L, T, cfg)
  X <- M + E

  svM <- svd(M, nu = r, nv = r)
  list(
    X = X,
    U_true = svM$u[, 1:r, drop = FALSE],
    V_true = svM$v[, 1:r, drop = FALSE]
  )
}

run_single_rep <- function(cfg, T) {
  dat <- simulate_panel(cfg, T)
  U_pca <- fit_pca_u(dat$X, cfg$r)
  U_wp <- fit_wpca_lag1_u(dat$X, cfg$r)
  V_pca <- fit_v_from_u(dat$X, U_pca, cfg$r)
  V_wp <- fit_v_from_u(dat$X, U_wp, cfg$r)

  rbind(
    data.frame(
      method = "PCA",
      u_error = subspace_error_simest(U_pca, dat$U_true, cfg$r),
      v_error = subspace_error_simest(V_pca, dat$V_true, cfg$r),
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = "WPCAlag1",
      u_error = subspace_error_simest(U_wp, dat$U_true, cfg$r),
      v_error = subspace_error_simest(V_wp, dat$V_true, cfg$r),
      stringsAsFactors = FALSE
    )
  )
}

summarize_metric <- function(raw_df, metric_name) {
  out <- aggregate(raw_df[[metric_name]],
                   by = list(scenario = raw_df$scenario, T = raw_df$T, method = raw_df$method),
                   FUN = function(x) c(mean = mean(x), sd = sd(x), se = sd(x) / sqrt(length(x))))
  stats <- if (is.matrix(out$x)) out$x else do.call(rbind, out$x)
  out$mean <- stats[, "mean"]
  out$sd <- stats[, "sd"]
  out$se <- stats[, "se"]
  out$x <- NULL
  out$metric <- metric_name
  out[, c("scenario", "T", "method", "metric", "mean", "sd", "se")]
}

plot_error_lines <- function(summary_df, out_file) {
  scenarios <- c("gaussian_noise", "student_t_noise")
  metrics <- c("u_error", "v_error")
  method_order <- c("PCA", "WPCAlag1")
  method_cols <- c("PCA" = "#223A5E", "WPCAlag1" = "#C44536")
  method_pch <- c("PCA" = 16, "WPCAlag1" = 15)
  method_lab <- c("PCA" = "PCA", "WPCAlag1" = "WPCA (lag-1)")

  png(filename = out_file, width = 1900, height = 1200, res = 200, bg = "white")
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)

  par(
    mfrow = c(2, 2),
    mar = c(4.6, 4.8, 2.2, 1.2),
    oma = c(0.8, 0.6, 0.4, 0.4),
    bg = "white"
  )

  for (metric in metrics) {
    for (sc in scenarios) {
      dat <- summary_df[summary_df$scenario == sc & summary_df$metric == metric, ]
      dat <- dat[order(dat$T, dat$method), ]
      Ts <- sort(unique(dat$T))

      y_rng <- range(dat$mean, na.rm = TRUE)
      y_span <- max(1e-8, y_rng[2] - y_rng[1])
      ylim <- c(max(0, y_rng[1] - 0.06 * y_span), y_rng[2] + 0.06 * y_span)

      panel_title <- if (sc == "gaussian_noise") "Gaussian Noise" else "Student-t Noise"
      ylab <- if (metric == "u_error") "U error" else "V error"

      plot(
        Ts, Ts,
        type = "n",
        xlab = "T",
        ylab = ylab,
        main = panel_title,
        ylim = ylim,
        xaxt = "n",
        bty = "l",
        cex.main = 1.0,
        font.main = 2
      )
      axis(1, at = Ts)

      for (m in method_order) {
        dm <- dat[dat$method == m, ]
        dm <- dm[order(dm$T), ]
        lines(dm$T, dm$mean, col = method_cols[[m]], lwd = 2.6)
        points(dm$T, dm$mean, pch = method_pch[[m]], col = method_cols[[m]], cex = 1.1)
      }

      legend(
        "right",
        legend = c(method_lab[["PCA"]], method_lab[["WPCAlag1"]]),
        col = method_cols[method_order],
        lwd = 2.6,
        pch = method_pch[method_order],
        pt.cex = 1.0,
        bty = "n",
        cex = 0.9
      )
    }
  }
}

# Tuning-ready defaults (at least 5 T points per line)
scenarios <- list(
  gaussian_noise = list(
    description = "Gaussian i.i.d. Z with SPD Sigma_C and Sigma_T = I",
    N = 100,
    r = 3,
    T_grid = c(120, 180, 260, 380, 560),
    phi_level = 0.90,
    signal_scale = 1.2,
    q_spike = 3,
    eig_signal = 2.5,
    eig_spike = 220,
    eig_base = 10,
    z_dist = "gaussian",
    z_df = 5L,
    sigma_t_type = "identity",
    rho_t = 0
  ),
  student_t_noise = list(
    description = "Student-t i.i.d. Z with SPD Sigma_C and AR(1)-Toeplitz Sigma_T",
    N = 100,
    r = 3,
    T_grid = c(150, 220, 320, 480, 720),
    phi_level = 0.90,
    signal_scale = 1.2,
    q_spike = 3,
    eig_signal = 2.5,
    eig_spike = 220,
    eig_base = 10,
    z_dist = "t5",
    z_df = 5L,
    sigma_t_type = "ar1",
    rho_t = 0.18
  )
)

cat("Start simulation (strict SPD noise: E = Sigma_C^{1/2} Z Sigma_T^{1/2})\n")
cat("n_rep =", n_rep, ", seed =", seed, "\n")

all_rows <- list()
idx <- 1L
for (scenario_name in names(scenarios)) {
  cfg <- scenarios[[scenario_name]]
  cat("\nScenario:", scenario_name, "-", cfg$description, "\n")
  for (T in cfg$T_grid) {
    cat("  T =", T, "\n")
    for (rep_id in seq_len(n_rep)) {
      rep_res <- run_single_rep(cfg, T)
      rep_res$scenario <- scenario_name
      rep_res$T <- T
      rep_res$rep <- rep_id
      all_rows[[idx]] <- rep_res
      idx <- idx + 1L
      if (rep_id %% 20 == 0) cat("    completed", rep_id, "replicates\n")
    }
  }
}

raw_results <- do.call(rbind, all_rows)
raw_results <- raw_results[, c("scenario", "T", "rep", "method", "u_error", "v_error")]
raw_results <- raw_results[order(raw_results$scenario, raw_results$T, raw_results$rep, raw_results$method), ]

summary_u <- summarize_metric(raw_results, "u_error")
summary_v <- summarize_metric(raw_results, "v_error")
summary_all <- rbind(summary_u, summary_v)
summary_all <- summary_all[order(summary_all$scenario, summary_all$metric, summary_all$T, summary_all$method), ]

write.csv(raw_results[, c("scenario", "T", "rep", "method", "u_error")],
          file.path(output_dir, "sim_raw_u_error.csv"), row.names = FALSE)
write.csv(summary_u[, c("scenario", "T", "method", "metric", "mean", "sd", "se")],
          file.path(output_dir, "sim_summary_u_error.csv"), row.names = FALSE)
write.csv(raw_results[, c("scenario", "T", "rep", "method", "v_error")],
          file.path(output_dir, "sim_raw_v_error.csv"), row.names = FALSE)
write.csv(summary_v[, c("scenario", "T", "method", "metric", "mean", "sd", "se")],
          file.path(output_dir, "sim_summary_v_error.csv"), row.names = FALSE)

plot_error_lines(summary_all, file.path(output_dir, "fig_error_lines.png"))

cat("\nSimulation finished. Saved files:\n")
cat("  - sim_raw_u_error.csv\n")
cat("  - sim_summary_u_error.csv\n")
cat("  - sim_raw_v_error.csv\n")
cat("  - sim_summary_v_error.csv\n")
cat("  - fig_error_lines.png\n\n")

cat("Mean errors by scenario/T/method:\n")
print(summary_all[, c("scenario", "T", "method", "metric", "mean", "se")])

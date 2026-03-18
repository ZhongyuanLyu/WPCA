#!/usr/bin/env Rscript

# Interpretive real-data add-on:
# 1) rolling selected gamma path
# 2) first-factor alignment with key macro series
# 3) correlation summary for AdaWPCA vs PCA
#
# For interpretability, we use transformed FRED series here.

args <- commandArgs(trailingOnly = TRUE)
seed <- if (length(args) >= 1) as.integer(args[1]) else 20260315L
set.seed(seed)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
  output_dir <- dirname(script_path)
} else {
  output_dir <- getwd()
}

base_dir <- "C:/Users/poetlyu/Nutstore/1/Nutstore/Research/Weighted-PCA/codes-for-codex/Realdata_FRED"

H_mat <- function(X) {
  diag(X) <- 0
  X
}

heteroPCA <- function(R, K, T0 = 10) {
  M <- H_mat(R %*% t(R))
  M_no_diag <- M
  if (T0 > 0) {
    for (t in 1:T0) {
      svd_res <- svd(M, nu = K, nv = K)
      if (K == 1) {
        M_bar <- svd_res$u %*% t(svd_res$v) * svd_res$d
      } else {
        M_bar <- svd_res$u %*% diag(svd_res$d[1:K], nrow = K, ncol = K) %*% t(svd_res$v)
      }
      M <- M_no_diag + diag(diag(M_bar))
    }
  }
  svd(M, nu = K, nv = K)$u
}

APCA <- function(X, gamma, r, symmetric = TRUE) {
  Times <- ncol(X)
  Q <- diag(gamma, nrow = Times)
  if (symmetric) {
    diag(Q[-nrow(Q), -1]) <- (1 - gamma)
    Q[lower.tri(Q)] <- t(Q)[lower.tri(Q)]
    Uhat <- svd(X %*% Q %*% t(X), nu = r)$u
  } else {
    diag(Q[-nrow(Q), -1]) <- (1 - gamma)
    Uhat <- svd(X %*% Q %*% t(X), nu = r)$u
  }
  Uhat
}

CV_APCA <- function(X, r, p_star = 0.8, grid_len = 10, J = 10) {
  N <- nrow(X)
  Times <- ncol(X)
  gamma_list <- seq(0, 1, length.out = grid_len)
  err_mat <- matrix(0, nrow = J, ncol = grid_len)
  for (jj in 1:J) {
    mask <- matrix(runif(N * Times) < p_star, nrow = N, ncol = Times)
    Y_train <- X * mask
    Y_test <- X * (1 - mask)
    for (j in 1:length(gamma_list)) {
      Ugammahat_j <- APCA(Y_train, gamma_list[j], r)
      err_mat[jj, j] <- sum((Y_test - (Ugammahat_j %*% t(Ugammahat_j) %*% Y_train / p_star) * (1 - mask))^2) / sum(1 - mask)
    }
  }
  err_vec <- colMeans(err_mat)
  best_gamma <- gamma_list[which.min(err_vec)]
  U_best <- APCA(X, best_gamma, r)
  list(U = U_best, gamma = best_gamma, err = err_vec)
}

transxf <- function(x, tcode) {
  n <- length(x)
  small <- 1e-06
  y <- rep(NA_real_, n)
  y1 <- rep(NA_real_, n)
  if (tcode == 1) {
    y <- x
  } else if (tcode == 2) {
    y[2:n] <- x[2:n] - x[1:(n - 1)]
  } else if (tcode == 3) {
    y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)]
  } else if (tcode == 4) {
    if (min(x, na.rm = TRUE) > small) y <- log(x)
  } else if (tcode == 5) {
    if (min(x, na.rm = TRUE) > small) {
      x <- log(x)
      y[2:n] <- x[2:n] - x[1:(n - 1)]
    }
  } else if (tcode == 6) {
    if (min(x, na.rm = TRUE) > small) {
      x <- log(x)
      y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)]
    }
  } else if (tcode == 7) {
    y1[2:n] <- (x[2:n] - x[1:(n - 1)]) / x[1:(n - 1)]
    y[3:n] <- y1[3:n] - y1[2:(n - 1)]
  }
  y
}

load_fred_balanced <- function(file, transform = TRUE) {
  raw <- utils::read.csv(file, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)
  header <- as.character(unlist(raw[1, ]))
  keep_nonempty <- !(is.na(header) | header == "")
  raw <- raw[, keep_nonempty, drop = FALSE]
  header <- header[keep_nonempty]

  data_rows <- grepl("^[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}$", raw[[1]])
  first_data <- which(data_rows)[1]
  dat <- raw[data_rows, , drop = FALSE]
  colnames(dat) <- header

  for (j in 2:ncol(dat)) {
    dat[[j]] <- suppressWarnings(as.numeric(dat[[j]]))
  }

  if (transform) {
    trow <- first_data - 1L
    tcodes <- suppressWarnings(as.integer(unlist(raw[trow, -1, drop = TRUE])))
    for (j in 2:ncol(dat)) {
      dat[[j]] <- transxf(dat[[j]], tcodes[j - 1])
    }
  }

  col_na_prop <- apply(is.na(dat[, -1, drop = FALSE]), 2, mean)
  keep_cols <- c(TRUE, col_na_prop < 0.05)
  data_select <- dat[, keep_cols, drop = FALSE]
  data_bal <- stats::na.omit(data_select)

  dates <- as.Date(data_bal[[1]], format = "%m/%d/%Y")
  rownames(data_bal) <- data_bal[[1]]
  X_bal <- data_bal[, -1, drop = FALSE]
  X <- t(as.matrix(X_bal))
  rownames(X) <- colnames(X_bal)

  list(X = X, dates = dates, var_names = rownames(X))
}

row_standardize <- function(X) {
  row_means <- rowMeans(X)
  row_sds <- apply(X, 1, stats::sd)
  row_sds[!is.finite(row_sds) | row_sds < 1e-8] <- 1
  Xs <- (X - row_means) / row_sds
  list(X = Xs, means = row_means, sds = row_sds)
}

zscore <- function(x) {
  s <- stats::sd(x)
  if (!is.finite(s) || s < 1e-8) return(x - mean(x))
  as.numeric((x - mean(x)) / s)
}

get_first_factor <- function(X, method = c("AdaWPCA", "PCA")) {
  method <- match.arg(method)
  if (method == "AdaWPCA") {
    res <- CV_APCA(X, r = 1L, p_star = 0.8, grid_len = 10, J = 10)
    Uhat <- res$U
    gamma <- res$gamma
  } else {
    Uhat <- svd(X, nu = 1, nv = 0)$u
    gamma <- 1
  }
  score <- as.numeric(t(Uhat) %*% X)
  list(score = score, gamma = gamma, U = Uhat)
}

select_gamma_path <- function(X, dates, horizon_step, n_windows, min_train) {
  Times <- ncol(X)
  last_end <- Times
  first_end <- max(min_train, last_end - horizon_step * (n_windows - 1L))
  train_ends <- seq(from = first_end, to = last_end, by = horizon_step)
  out <- data.frame(date = dates[train_ends], train_end = train_ends, gamma = NA_real_)
  for (i in seq_along(train_ends)) {
    Xi <- X[, 1:train_ends[i], drop = FALSE]
    out$gamma[i] <- CV_APCA(Xi, r = 1L, p_star = 0.8, grid_len = 10, J = 10)$gamma
  }
  out
}

align_sign <- function(factor, target) {
  if (stats::cor(factor, target, use = "complete.obs") < 0) {
    return(-factor)
  }
  factor
}

make_corr_table <- function(factor_df, benchmarks) {
  rows <- list()
  idx <- 1L
  for (m in c("AdaWPCA", "PCA")) {
    f <- factor_df[[m]]
    for (bname in names(benchmarks)) {
      y <- benchmarks[[bname]]
      rows[[idx]] <- data.frame(
        method = m,
        variable = bname,
        corr0 = stats::cor(f, y, use = "complete.obs"),
        corr_lead1 = stats::cor(f[-length(f)], y[-1], use = "complete.obs"),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

plot_interpretive <- function(md_gamma, qd_gamma, md_dates, md_series, qd_dates, qd_series, out_file) {
  clip_for_plot <- function(x) {
    qs <- stats::quantile(x, probs = c(0.01, 0.99), na.rm = TRUE)
    pmin(pmax(x, qs[1]), qs[2])
  }
  smooth_for_plot <- function(x, k) {
    as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
  }

  png(filename = out_file, width = 2000, height = 1300, res = 180, bg = "white")
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)

  par(mfrow = c(2, 2), mar = c(4.2, 4.8, 2.2, 1.2), bg = "white")

  legend_labels <- c(
    "AdaWPCA factor",
    "Payroll employment",
    "Inverse unemployment rate"
  )

  plot(md_gamma$date, md_gamma$gamma, type = "b", pch = 16, col = "#C44536", lwd = 2.2,
       ylim = c(-0.02, 1.02), xlab = "", ylab = expression(hat(gamma)), main = "FRED-MD: Rolling Selected Gamma", bty = "l")
  abline(h = 1, lty = 2, col = "#666666")

  plot(qd_gamma$date, qd_gamma$gamma, type = "b", pch = 16, col = "#C44536", lwd = 2.2,
       ylim = c(-0.02, 1.02), xlab = "", ylab = expression(hat(gamma)), main = "FRED-QD: Rolling Selected Gamma", bty = "l")
  abline(h = 1, lty = 2, col = "#666666")

  md_factor_plot <- clip_for_plot(smooth_for_plot(md_series$factor, 12))
  md_activity_plot <- clip_for_plot(smooth_for_plot(md_series$activity, 12))
  md_slack_plot <- clip_for_plot(smooth_for_plot(md_series$slack, 12))
  ylim_md <- range(c(md_factor_plot, md_activity_plot, md_slack_plot), na.rm = TRUE)
  plot(md_dates, md_factor_plot, type = "l", lwd = 2.4, col = "#C44536", ylim = ylim_md,
       xlab = "", ylab = "Standardized index (z-score)", main = "FRED-MD: Factor and Macro Alignment", bty = "l")
  lines(md_dates, md_activity_plot, lwd = 2.0, col = "#223A5E")
  lines(md_dates, md_slack_plot, lwd = 2.0, col = "#5C7C5A")
  legend("topleft", legend = legend_labels, col = c("#C44536", "#223A5E", "#5C7C5A"),
         lwd = c(2.4, 2.0, 2.0), bty = "n", cex = 0.8, bg = "white")

  qd_factor_plot <- clip_for_plot(smooth_for_plot(qd_series$factor, 4))
  qd_activity_plot <- clip_for_plot(smooth_for_plot(qd_series$activity, 4))
  qd_slack_plot <- clip_for_plot(smooth_for_plot(qd_series$slack, 4))
  ylim_qd <- range(c(qd_factor_plot, qd_activity_plot, qd_slack_plot), na.rm = TRUE)
  plot(qd_dates, qd_factor_plot, type = "l", lwd = 2.4, col = "#C44536", ylim = ylim_qd,
       xlab = "", ylab = "Standardized index (z-score)", main = "FRED-QD: Factor and Macro Alignment", bty = "l")
  lines(qd_dates, qd_activity_plot, lwd = 2.0, col = "#223A5E")
  lines(qd_dates, qd_slack_plot, lwd = 2.0, col = "#5C7C5A")
  legend("bottomleft", legend = legend_labels, col = c("#C44536", "#223A5E", "#5C7C5A"),
         lwd = c(2.4, 2.0, 2.0), bty = "n", cex = 0.8, bg = "white")
}

md <- load_fred_balanced(file.path(base_dir, "FRED-MD-2025-02.csv"), transform = TRUE)
qd <- load_fred_balanced(file.path(base_dir, "FRED-QD-2025-02.csv"), transform = TRUE)

md_std <- row_standardize(md$X)$X
qd_std <- row_standardize(qd$X)$X

md_factor_ada <- get_first_factor(md_std, "AdaWPCA")
md_factor_pca <- get_first_factor(md_std, "PCA")
qd_factor_ada <- get_first_factor(qd_std, "AdaWPCA")
qd_factor_pca <- get_first_factor(qd_std, "PCA")

md_activity <- zscore(md_std["PAYEMS", ])
md_output <- zscore(md_std["INDPRO", ])
md_unrate <- zscore(-md_std["UNRATE", ])
qd_activity <- zscore(qd_std["PAYEMS", ])
qd_output <- zscore(qd_std["GDPC1", ])
qd_unrate <- zscore(-qd_std["UNRATE", ])

md_factor_ada_score <- align_sign(zscore(md_factor_ada$score), md_activity)
md_factor_pca_score <- align_sign(zscore(md_factor_pca$score), md_activity)
qd_factor_ada_score <- align_sign(zscore(qd_factor_ada$score), qd_activity)
qd_factor_pca_score <- align_sign(zscore(qd_factor_pca$score), qd_activity)

md_gamma_path <- select_gamma_path(md_std, md$dates, horizon_step = 24L, n_windows = 10L, min_train = 240L)
qd_gamma_path <- select_gamma_path(qd_std, qd$dates, horizon_step = 8L, n_windows = 10L, min_train = 80L)

md_corr <- make_corr_table(
  factor_df = list(AdaWPCA = md_factor_ada_score, PCA = md_factor_pca_score),
  benchmarks = list(PAYEMS = md_activity, INDPRO = md_output, minus_UNRATE = md_unrate)
)
md_corr$dataset <- "FRED-MD"

qd_corr <- make_corr_table(
  factor_df = list(AdaWPCA = qd_factor_ada_score, PCA = qd_factor_pca_score),
  benchmarks = list(PAYEMS = qd_activity, GDPC1 = qd_output, minus_UNRATE = qd_unrate)
)
qd_corr$dataset <- "FRED-QD"

corr_table <- rbind(md_corr, qd_corr)
corr_table <- corr_table[, c("dataset", "method", "variable", "corr0", "corr_lead1")]

gamma_summary <- data.frame(
  dataset = c("FRED-MD", "FRED-QD"),
  gamma_mean = c(mean(md_gamma_path$gamma), mean(qd_gamma_path$gamma)),
  gamma_min = c(min(md_gamma_path$gamma), min(qd_gamma_path$gamma)),
  gamma_max = c(max(md_gamma_path$gamma), max(qd_gamma_path$gamma)),
  share_below_one = c(mean(md_gamma_path$gamma < 0.999), mean(qd_gamma_path$gamma < 0.999))
)

plot_interpretive(
  md_gamma = md_gamma_path,
  qd_gamma = qd_gamma_path,
  md_dates = md$dates,
  md_series = list(factor = md_factor_ada_score, activity = md_activity, slack = md_unrate),
  qd_dates = qd$dates,
  qd_series = list(factor = qd_factor_ada_score, activity = qd_activity, slack = qd_unrate),
  out_file = file.path(output_dir, "fig_interpretive_addon.png")
)

write.csv(md_gamma_path, file.path(output_dir, "md_gamma_path.csv"), row.names = FALSE)
write.csv(qd_gamma_path, file.path(output_dir, "qd_gamma_path.csv"), row.names = FALSE)
write.csv(gamma_summary, file.path(output_dir, "gamma_summary.csv"), row.names = FALSE)
write.csv(corr_table, file.path(output_dir, "factor_macro_correlations.csv"), row.names = FALSE)

cat("Saved files:\n")
cat("  - fig_interpretive_addon.png\n")
cat("  - md_gamma_path.csv\n")
cat("  - qd_gamma_path.csv\n")
cat("  - gamma_summary.csv\n")
cat("  - factor_macro_correlations.csv\n\n")

cat("Gamma summary:\n")
print(gamma_summary)
cat("\nFactor/macro correlations:\n")
print(corr_table)

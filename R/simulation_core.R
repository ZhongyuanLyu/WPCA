generate_signal_panel <- function(N, Times, r = 3L, phi_level = 0.9) {
  B0 <- matrix(0, nrow = r)
  diagPhivec <- seq(phi_level, phi_level, length.out = r)
  U_phi <- svd(matrix(stats::rnorm(r * r), nrow = r))$u
  V_phi <- svd(matrix(stats::rnorm(r * r), nrow = r))$u
  Phi <- U_phi %*% diag(diagPhivec) %*% t(V_phi)
  sigvare <- diag(1, r, r) - Phi %*% t(Phi)
  Fmat <- VAR1.sim(B0, Phi, Times, sigvare)
  Fmat <- stats::ts(Fmat)

  svd_L <- svd(matrix(stats::rnorm(N * r), nrow = N))
  L <- svd_L$u %*% diag(svd_L$d)
  M <- L %*% t(Fmat)

  list(
    Fmat = Fmat,
    L = L,
    M = M,
    U_true = svd(M, nu = r)$u,
    V_true = svd(M, nv = r)$v
  )
}

draw_sigma_c <- function(N, rho = 0.6, covariance = c("offdiag", "diagonal", "isotropic")) {
  covariance <- match.arg(covariance)
  if (covariance == "isotropic") {
    return(10 * diag(N))
  }
  if (covariance == "diagonal") {
    return(diag(stats::runif(N, 1, 20)))
  }
  diag(stats::runif(N, 1, 20)) + rho * (rep(1, N) %*% t(rep(1, N)) - diag(N))
}

make_weight_matrix <- function(Times, gamma, symmetric = TRUE) {
  Q <- diag(gamma, nrow = Times)
  diag(Q[-nrow(Q), -1]) <- 1 - gamma
  if (symmetric) {
    Q[lower.tri(Q)] <- t(Q)[lower.tri(Q)]
  }
  Q
}

simulate_estimation_errors <- function(N, Times, rho, covariance, N_rep = 100L, r = 3L) {
  err_res <- matrix(NA_real_, nrow = N_rep, ncol = 6)
  for (ii in seq_len(N_rep)) {
    signal <- generate_signal_panel(N = N, Times = Times, r = r)
    sigeC <- draw_sigma_c(N = N, rho = rho, covariance = covariance)
    E <- t(MASS::mvrnorm(Times, rep(0, N), sigeC))
    X <- signal$M + E

    U_PCA <- svd(X, nu = r, nv = r)$u
    U_HPCA <- heteroPCA(X, r, 10)
    U_APCA <- CV_APCA(X, r, p_star = 0.8, J = 10)$U

    err_res[ii, 1] <- sqrt(sum((U_PCA %*% t(U_PCA) - signal$U_true %*% t(signal$U_true))^2) / sqrt(r))
    err_res[ii, 2] <- sqrt(sum((U_HPCA %*% t(U_HPCA) - signal$U_true %*% t(signal$U_true))^2) / sqrt(r))
    err_res[ii, 3] <- sqrt(sum((U_APCA %*% t(U_APCA) - signal$U_true %*% t(signal$U_true))^2) / sqrt(r))

    V_PCA <- svd(U_PCA %*% t(U_PCA) %*% X, nu = r, nv = r)$v
    V_HPCA <- svd(U_HPCA %*% t(U_HPCA) %*% X, nu = r, nv = r)$v
    V_APCA <- svd(U_APCA %*% t(U_APCA) %*% X, nu = r, nv = r)$v

    err_res[ii, 4] <- sqrt(sum((V_PCA %*% t(V_PCA) - signal$V_true %*% t(signal$V_true))^2) / sqrt(r))
    err_res[ii, 5] <- sqrt(sum((V_HPCA %*% t(V_HPCA) - signal$V_true %*% t(signal$V_true))^2) / sqrt(r))
    err_res[ii, 6] <- sqrt(sum((V_APCA %*% t(V_APCA) - signal$V_true %*% t(signal$V_true))^2) / sqrt(r))
  }
  err_res
}

build_estimation_table <- function(N_values, T_values, rho, covariance, N_rep = 100L) {
  rows <- list()
  idx <- 1L
  for (N in N_values) {
    for (Times in T_values) {
      err_res <- simulate_estimation_errors(N = N, Times = Times, rho = rho, covariance = covariance, N_rep = N_rep)
      rows[[idx]] <- data.frame(
        Panel = "U",
        N = N,
        T = Times,
        AdaWPCA = format_mean_sd(err_res[, 3]),
        PCA = format_mean_sd(err_res[, 1]),
        HeteroPCA = format_mean_sd(err_res[, 2]),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        Panel = "V",
        N = N,
        T = Times,
        AdaWPCA = format_mean_sd(err_res[, 6]),
        PCA = format_mean_sd(err_res[, 4]),
        HeteroPCA = format_mean_sd(err_res[, 5]),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

simulate_cv_selection <- function(N, Times, rho, covariance, N_rep = 100L, r = 3L) {
  err_res <- matrix(NA_real_, nrow = N_rep, ncol = 13)
  gamma_list <- seq(0, 1, length.out = 10)

  for (ii in seq_len(N_rep)) {
    signal <- generate_signal_panel(N = N, Times = Times, r = r)
    sigeC <- draw_sigma_c(N = N, rho = rho, covariance = covariance)
    E <- t(MASS::mvrnorm(Times, rep(0, N), sigeC))
    X <- signal$M + E

    U_PCA <- svd(X, nu = r, nv = r)$u
    U_HPCA <- heteroPCA(X, r, 10)
    APCA_res <- CV_APCA(X, r, p_star = 0.8, J = 10)
    U_APCA <- APCA_res$U

    true_err <- rep(NA_real_, length(gamma_list))
    for (j in seq_along(gamma_list)) {
      Ugammahat_j <- APCA(X, gamma_list[j], r)
      true_err[j] <- sqrt(sum((Ugammahat_j %*% t(Ugammahat_j) - signal$U_true %*% t(signal$U_true))^2) / sqrt(r))
    }

    ind_true <- order(true_err)[1:5]
    ind_bad <- order(true_err)[8:10]

    err_res[ii, 1] <- sqrt(sum((U_PCA %*% t(U_PCA) - signal$U_true %*% t(signal$U_true))^2) / sqrt(r))
    err_res[ii, 2] <- sqrt(sum((U_HPCA %*% t(U_HPCA) - signal$U_true %*% t(signal$U_true))^2) / sqrt(r))
    err_res[ii, 3] <- sqrt(sum((U_APCA %*% t(U_APCA) - signal$U_true %*% t(signal$U_true))^2) / sqrt(r))
    err_res[ii, 4] <- true_err[ind_true[1]]
    err_res[ii, 5] <- APCA_res$gamma
    err_res[ii, 6:10] <- gamma_list[ind_true]
    err_res[ii, 11:13] <- gamma_list[ind_bad]
  }

  err_res
}

build_cv_table <- function(N_values, T_values, rho_values, covariance_values, setting_labels, N_rep = 100L) {
  rows <- list()
  idx <- 1L
  for (kk in seq_along(setting_labels)) {
    for (N in N_values) {
      for (Times in T_values) {
        mat <- simulate_cv_selection(
          N = N,
          Times = Times,
          rho = rho_values[[kk]],
          covariance = covariance_values[[kk]],
          N_rep = N_rep
        )
        top3 <- mean(apply(mat, 1, function(x) x[5] %in% x[6:8]))
        nonbottom3 <- mean(1 - apply(mat, 1, function(x) x[5] %in% x[11:13]))
        rows[[idx]] <- data.frame(
          Scenario = setting_labels[[kk]],
          N = N,
          T = Times,
          Top3 = round(top3, 2),
          NonBottom3 = round(nonbottom3, 2),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  do.call(rbind, rows)
}

simulate_loading_distribution <- function(N, Times, rho = 0.6, gamma = 0, N_rep = 500L, r = 3L) {
  signal <- generate_signal_panel(N = N, Times = Times, r = r)
  sigeC <- draw_sigma_c(N = N, rho = rho, covariance = "offdiag")
  Q <- make_weight_matrix(Times = Times, gamma = gamma, symmetric = TRUE)
  MQMt <- signal$M %*% Q %*% t(signal$M)
  ind <- as.integer(0.5 * N)

  LResnorm <- matrix(NA_real_, nrow = r, ncol = N_rep)
  for (ii in seq_len(N_rep)) {
    E <- t(MASS::mvrnorm(Times, rep(0, N), sigeC))
    X <- signal$M + E
    U_WPCA <- APCA(X, gamma, r, symmetric = TRUE)
    svd_nUUtX <- svd(U_WPCA %*% t(U_WPCA) %*% X / sqrt(Times), nu = r, nv = r)
    rotation_cov <- compute_rotation_and_covariance_L(
      svd_nUUtX$u,
      svd_nUUtX$v,
      signal$Fmat,
      signal$L,
      signal$M / sqrt(Times),
      MQMt,
      sigeC,
      diag(Times),
      Q,
      ind
    )
    LRes <- svd_nUUtX$u %*% diag(svd_nUUtX$d[1:r]) - signal$L %*% rotation_cov$R_LQ
    LResnorm[, ii] <- LRes[ind, ] %*% solve(expm::sqrtm(rotation_cov$Sigma_L))
  }

  LResnorm
}

simulate_factor_distribution <- function(N, Times, rho = 0.6, gamma = 0, N_rep = 500L, r = 3L) {
  signal <- generate_signal_panel(N = N, Times = Times, r = r)
  sigeC <- draw_sigma_c(N = N, rho = rho, covariance = "offdiag")
  ind <- as.integer(0.5 * Times)

  FResnorm <- matrix(NA_real_, nrow = r, ncol = N_rep)
  for (ii in seq_len(N_rep)) {
    E <- t(MASS::mvrnorm(Times, rep(0, N), sigeC))
    X <- signal$M + E
    U_WPCA <- APCA(X, gamma, r)
    V_WPCA <- svd(U_WPCA %*% t(U_WPCA) %*% X, nu = r, nv = r)$v
    rotation_cov <- compute_rotation_and_covariance_F(
      V_WPCA,
      signal$Fmat,
      signal$M / sqrt(Times),
      sigeC,
      diag(Times),
      ind
    )
    FRes <- sqrt(Times) * V_WPCA - signal$Fmat %*% rotation_cov$R_FQ
    FResnorm[, ii] <- solve(expm::sqrtm(rotation_cov$Sigma_F)) %*% FRes[ind, ]
  }

  FResnorm
}

plot_qq <- function(x, title_text, y_label) {
  qq_line_layer <- if (utils::packageVersion("ggplot2") >= "3.4.0") {
    ggplot2::stat_qq_line(color = "steelblue", linewidth = 1)
  } else {
    ggplot2::stat_qq_line(color = "steelblue", size = 1)
  }

  ggplot2::ggplot(data.frame(x = x), ggplot2::aes(sample = x)) +
    ggplot2::stat_qq(size = 2) +
    qq_line_layer +
    ggplot2::labs(title = title_text, x = "Theoretical Quantiles", y = y_label) +
    ggplot2::theme_minimal(base_size = 30) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.background = ggplot2::element_blank()
    )
}

plot_hist_standardized <- function(x, title_text, x_label) {
  hist_mapping <- if (exists("after_stat", where = asNamespace("ggplot2"), inherits = FALSE)) {
    ggplot2::aes(y = ggplot2::after_stat(density))
  } else {
    ggplot2::aes(y = ..density..)
  }

  stat_function_layer <- if (utils::packageVersion("ggplot2") >= "3.4.0") {
    ggplot2::stat_function(fun = stats::dnorm, args = list(mean = 0, sd = 1), color = "black", linewidth = 1)
  } else {
    ggplot2::stat_function(fun = stats::dnorm, args = list(mean = 0, sd = 1), color = "black", size = 1)
  }

  ggplot2::ggplot(data.frame(EstimationError = x), ggplot2::aes(x = EstimationError)) +
    ggplot2::geom_histogram(mapping = hist_mapping, bins = 12, fill = "royalblue3", color = "white") +
    stat_function_layer +
    ggplot2::labs(title = title_text, x = x_label, y = "Frequency") +
    ggplot2::theme_minimal(base_size = 30) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.background = ggplot2::element_blank()
    )
}

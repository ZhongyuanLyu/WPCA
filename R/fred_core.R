transxf_fred <- function(x, tcode) {
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
    if (min(x, na.rm = TRUE) > small) {
      y <- log(x)
    }
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

load_fred_csv_panel <- function(file, transform = FALSE) {
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
      dat[[j]] <- transxf_fred(dat[[j]], tcodes[j - 1])
    }
  }

  col_na_prop <- apply(is.na(dat[, -1, drop = FALSE]), 2, mean)
  keep_cols <- c(TRUE, col_na_prop < 0.05)
  data_select <- dat[, keep_cols, drop = FALSE]
  data_bal <- stats::na.omit(data_select)
  rownames(data_bal) <- data_bal[[1]]
  X_bal <- data_bal[, -1, drop = FALSE]
  X <- t(as.matrix(X_bal))
  rownames(X) <- colnames(X_bal)

  list(
    data = data_bal,
    X = X,
    dates = as.Date(data_bal[[1]], format = "%m/%d/%Y")
  )
}

prepare_fred_panel <- function(data_dir, dataset = c("FRED-MD", "FRED-QD"), transform = FALSE) {
  dataset <- match.arg(dataset)
  file_name <- if (dataset == "FRED-MD") "FRED-MD-2025-02.csv" else "FRED-QD-2025-02.csv"
  load_fred_csv_panel(file = file.path(data_dir, file_name), transform = transform)
}

one_fred_sim <- function(X, train_ratio, r) {
  N <- nrow(X)
  Times <- ncol(X)
  err <- rep(NA_real_, 3)

  train_mat <- matrix(stats::runif(N * Times) < train_ratio, nrow = N, ncol = Times)
  X_train <- X * train_mat
  X_test <- X * (1 - train_mat)

  svd_res <- svd(X_train, nu = r, nv = r)
  if (r > 1) {
    X_PCA <- svd_res$u %*% diag(svd_res$d[1:r]) %*% t(svd_res$v) / train_ratio
  } else {
    X_PCA <- svd_res$u %*% t(svd_res$v) * svd_res$d[r] / train_ratio
  }

  U_heteroPCA <- heteroPCA(X_train, r, T0 = 10)
  res_APCA <- CV_APCA(X_train, r, p_star = 0.8, J = 10)

  err[1] <- sum((X_test - (res_APCA$U %*% t(res_APCA$U) %*% X_train) / train_ratio * (1 - train_mat))^2) / sum(X_test^2)
  err[2] <- sum((X_test - X_PCA * (1 - train_mat))^2) / sum(X_test^2)
  err[3] <- sum((X_test - (U_heteroPCA %*% t(U_heteroPCA) %*% X_train) / train_ratio * (1 - train_mat))^2) / sum(X_test^2)
  err
}

build_fred_reconstruction_table <- function(data_dir, ranks = 1:3, train_ratios = c(0.7, 0.8, 0.9), n_rep = 100L) {
  rows <- list()
  idx <- 1L
  for (dataset in c("FRED-MD", "FRED-QD")) {
    panel <- prepare_fred_panel(data_dir = data_dir, dataset = dataset, transform = FALSE)
    X <- panel$X
    for (r in ranks) {
      for (train_ratio in train_ratios) {
        errs <- replicate(n_rep, one_fred_sim(X, train_ratio = train_ratio, r = r))
        rows[[idx]] <- data.frame(
          Dataset = dataset,
          Rank = r,
          q_tr = sprintf("%.1f", train_ratio),
          Method = c("AdaWPCA", "PCA", "HeteroPCA"),
          MeanError = round(rowMeans(errs), 3),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  do.call(rbind, rows)
}

row_standardize <- function(X) {
  row_means <- rowMeans(X)
  row_sds <- apply(X, 1, stats::sd)
  row_sds[!is.finite(row_sds) | row_sds < 1e-8] <- 1
  list(
    X = (X - row_means) / row_sds,
    means = row_means,
    sds = row_sds
  )
}

zscore <- function(x) {
  s <- stats::sd(x)
  if (!is.finite(s) || s < 1e-8) {
    return(as.numeric(x - mean(x)))
  }
  as.numeric((x - mean(x)) / s)
}

align_sign <- function(factor, target) {
  if (stats::cor(factor, target, use = "complete.obs") < 0) {
    return(-factor)
  }
  factor
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
  list(
    score = as.numeric(t(Uhat) %*% X),
    gamma = gamma,
    U = Uhat
  )
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

make_corr_table <- function(factor_df, benchmarks) {
  rows <- list()
  idx <- 1L
  for (method in c("AdaWPCA", "PCA")) {
    f <- factor_df[[method]]
    for (var_name in names(benchmarks)) {
      y <- benchmarks[[var_name]]
      rows[[idx]] <- data.frame(
        method = method,
        variable = var_name,
        corr0 = stats::cor(f, y, use = "complete.obs"),
        corr_lead1 = stats::cor(f[-length(f)], y[-1], use = "complete.obs"),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

plot_interpretive_figure <- function(md_gamma, qd_gamma, md_dates, md_series, qd_dates, qd_series, out_file) {
  clip_for_plot <- function(x) {
    qs <- stats::quantile(x, probs = c(0.01, 0.99), na.rm = TRUE)
    pmin(pmax(x, qs[1]), qs[2])
  }

  smooth_for_plot <- function(x, k) {
    as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
  }

  grDevices::png(filename = out_file, width = 2000, height = 1300, res = 180, bg = "white")
  op <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(op)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.8, 2.2, 1.2), bg = "white")

  legend_labels <- c("AdaWPCA factor", "Payroll employment", "Inverse unemployment rate")

  graphics::plot(md_gamma$date, md_gamma$gamma, type = "b", pch = 16, col = "#C44536", lwd = 2.2,
                 ylim = c(-0.02, 1.02), xlab = "", ylab = expression(hat(gamma)),
                 main = "FRED-MD: Rolling Selected Gamma", bty = "l")
  graphics::abline(h = 1, lty = 2, col = "#666666")

  graphics::plot(qd_gamma$date, qd_gamma$gamma, type = "b", pch = 16, col = "#C44536", lwd = 2.2,
                 ylim = c(-0.02, 1.02), xlab = "", ylab = expression(hat(gamma)),
                 main = "FRED-QD: Rolling Selected Gamma", bty = "l")
  graphics::abline(h = 1, lty = 2, col = "#666666")

  md_factor_plot <- clip_for_plot(smooth_for_plot(md_series$factor, 12))
  md_activity_plot <- clip_for_plot(smooth_for_plot(md_series$activity, 12))
  md_slack_plot <- clip_for_plot(smooth_for_plot(md_series$slack, 12))
  ylim_md <- range(c(md_factor_plot, md_activity_plot, md_slack_plot), na.rm = TRUE)
  graphics::plot(md_dates, md_factor_plot, type = "l", lwd = 2.4, col = "#C44536", ylim = ylim_md,
                 xlab = "", ylab = "Standardized index (z-score)", main = "FRED-MD: Factor and Macro Alignment", bty = "l")
  graphics::lines(md_dates, md_activity_plot, lwd = 2, col = "#223A5E")
  graphics::lines(md_dates, md_slack_plot, lwd = 2, col = "#5C7C5A")
  graphics::legend("topleft", legend = legend_labels, col = c("#C44536", "#223A5E", "#5C7C5A"),
                   lwd = c(2.4, 2, 2), bty = "n", cex = 0.8, bg = "white")

  qd_factor_plot <- clip_for_plot(smooth_for_plot(qd_series$factor, 4))
  qd_activity_plot <- clip_for_plot(smooth_for_plot(qd_series$activity, 4))
  qd_slack_plot <- clip_for_plot(smooth_for_plot(qd_series$slack, 4))
  ylim_qd <- range(c(qd_factor_plot, qd_activity_plot, qd_slack_plot), na.rm = TRUE)
  graphics::plot(qd_dates, qd_factor_plot, type = "l", lwd = 2.4, col = "#C44536", ylim = ylim_qd,
                 xlab = "", ylab = "Standardized index (z-score)", main = "FRED-QD: Factor and Macro Alignment", bty = "l")
  graphics::lines(qd_dates, qd_activity_plot, lwd = 2, col = "#223A5E")
  graphics::lines(qd_dates, qd_slack_plot, lwd = 2, col = "#5C7C5A")
  graphics::legend("bottomleft", legend = legend_labels, col = c("#C44536", "#223A5E", "#5C7C5A"),
                   lwd = c(2.4, 2, 2), bty = "n", cex = 0.8, bg = "white")
}

build_interpretive_outputs <- function(data_dir) {
  md <- prepare_fred_panel(data_dir = data_dir, dataset = "FRED-MD", transform = TRUE)
  qd <- prepare_fred_panel(data_dir = data_dir, dataset = "FRED-QD", transform = TRUE)

  md_scaled <- row_standardize(md$X)
  qd_scaled <- row_standardize(qd$X)

  md_ada <- get_first_factor(md_scaled$X, method = "AdaWPCA")
  md_pca <- get_first_factor(md_scaled$X, method = "PCA")
  qd_ada <- get_first_factor(qd_scaled$X, method = "AdaWPCA")
  qd_pca <- get_first_factor(qd_scaled$X, method = "PCA")

  extract_series <- function(obj, code) {
    idx <- match(code, rownames(obj$X))
    if (is.na(idx)) {
      stop("Series not found in balanced panel: ", code)
    }
    as.numeric(obj$X[idx, ])
  }

  md_activity <- zscore(extract_series(md_scaled, "PAYEMS"))
  md_output <- zscore(extract_series(md_scaled, "INDPRO"))
  md_slack <- zscore(-extract_series(md_scaled, "UNRATE"))
  qd_activity <- zscore(extract_series(qd_scaled, "PAYEMS"))
  qd_output <- zscore(extract_series(qd_scaled, "GDPC1"))
  qd_slack <- zscore(-extract_series(qd_scaled, "UNRATE"))

  md_factor <- align_sign(zscore(md_ada$score), md_activity)
  md_factor_pca <- align_sign(zscore(md_pca$score), md_activity)
  qd_factor <- align_sign(zscore(qd_ada$score), qd_activity)
  qd_factor_pca <- align_sign(zscore(qd_pca$score), qd_activity)

  list(
    md_gamma = select_gamma_path(md_scaled$X, md$dates, horizon_step = 12L, n_windows = 8L, min_train = 120L),
    qd_gamma = select_gamma_path(qd_scaled$X, qd$dates, horizon_step = 4L, n_windows = 8L, min_train = 60L),
    md_dates = md$dates,
    qd_dates = qd$dates,
    md_series = list(factor = md_factor, activity = md_activity, slack = md_slack),
    qd_series = list(factor = qd_factor, activity = qd_activity, slack = qd_slack),
    corr_table = rbind(
      transform(make_corr_table(list(AdaWPCA = md_factor, PCA = md_factor_pca),
                                list(PAYEMS = md_activity, INDPRO = md_output, minus_UNRATE = md_slack)),
                dataset = "FRED-MD"),
      transform(make_corr_table(list(AdaWPCA = qd_factor, PCA = qd_factor_pca),
                                list(PAYEMS = qd_activity, GDPC1 = qd_output, minus_UNRATE = qd_slack)),
                dataset = "FRED-QD")
    )
  )
}

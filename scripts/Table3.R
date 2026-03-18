get_script_path <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
  }
  frame_files <- vapply(sys.frames(), function(env) {
    if (!is.null(env$ofile)) env$ofile else NA_character_
  }, FUN.VALUE = character(1))
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files)) {
    return(normalizePath(frame_files[[length(frame_files)]], winslash = "/", mustWork = TRUE))
  }
  NA_character_
}

locate_repo_root <- function() {
  script_path <- get_script_path()
  if (!is.na(script_path)) {
    script_dir <- dirname(script_path)
    candidates <- c(
      normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE),
      normalizePath(script_dir, winslash = "/", mustWork = FALSE)
    )
  } else {
    cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    candidates <- c(
      cwd,
      normalizePath(file.path(cwd, ".."), winslash = "/", mustWork = FALSE)
    )
  }
  for (cand in unique(candidates)) {
    if (file.exists(file.path(cand, "R", "helpers.R")) && dir.exists(file.path(cand, "scripts"))) {
      return(cand)
    }
  }
  stop("Cannot locate the replication package root. Run the script from the package root, the scripts directory, or via Rscript.")
}

repo_root <- locate_repo_root()

source(file.path(repo_root, "R", "helpers.R"))
source(file.path(repo_root, "R", "functions.R"))
source(file.path(repo_root, "R", "simulation_core.R"))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) && args[1] %in% c("cached", "recompute")) args[1] else "cached"

if (mode == "cached") {
  cache_dir <- file.path(repo_root, "original", "codes")
  load_res_list <- function(file) {
    env <- new.env(parent = emptyenv())
    load(file, envir = env)
    get("res_list", envir = env)
  }
  valid_res_list <- function(x, expected_len) {
    is.list(x) &&
      length(x) == expected_len &&
      all(vapply(x, function(item) is.matrix(item) && is.numeric(item), logical(1)))
  }

  files <- c(`100` = "sim_est_random_N100_diag.Rdata", `200` = "sim_est_random_N200_diag.Rdata")
  Ts <- c(250L, 500L, 750L, 1000L)
  rows <- list()
  idx <- 1L
  for (N in names(files)) {
    res_list <- load_res_list(file.path(cache_dir, files[[N]]))
    if (!valid_res_list(res_list, length(Ts))) {
      set.seed(1234L + as.integer(N))
      res_list <- lapply(Ts, function(T_now) {
        simulate_estimation_errors(
          N = as.integer(N),
          Times = T_now,
          rho = 0.6,
          covariance = "diagonal",
          N_rep = 20L
        )
      })
    }
    for (i in seq_along(Ts)) {
      mat <- res_list[[i]]
      rows[[idx]] <- data.frame(
        Panel = "U",
        N = as.integer(N),
        T = Ts[i],
        AdaWPCA = format_mean_sd(mat[, 3]),
        PCA = format_mean_sd(mat[, 1]),
        HeteroPCA = format_mean_sd(mat[, 2]),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        Panel = "V",
        N = as.integer(N),
        T = Ts[i],
        AdaWPCA = format_mean_sd(mat[, 6]),
        PCA = format_mean_sd(mat[, 4]),
        HeteroPCA = format_mean_sd(mat[, 5]),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  table3 <- do.call(rbind, rows)
} else {
  n_rep <- get_int_arg(args, 2, 100L)
  seed <- get_int_arg(args, 3, 1234L)
  set.seed(seed)
  table3 <- build_estimation_table(
    N_values = c(100L, 200L),
    T_values = c(250L, 500L, 750L, 1000L),
    rho = 0.6,
    covariance = "diagonal",
    N_rep = n_rep
  )
}

out_file <- write_table_csv(table3, repo_root, "Table3_sim_est_diag.csv")
report_saved_files(out_file)

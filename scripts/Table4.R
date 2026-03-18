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

  Ts <- c(250L, 500L, 750L, 1000L)
  files <- list(
    "Diagonal Sigma_C" = c(`100` = "sim_cv_random_N100_diag20.Rdata", `200` = "sim_cv_random_N200_diag20.Rdata"),
    "Off-diagonal Sigma_C (rho = 0.5)" = c(`100` = "sim_cv_random_N100_rho0_5.Rdata", `200` = "sim_cv_random_N200_rho0_5.Rdata")
  )

  rows <- list()
  idx <- 1L
  for (setting_name in names(files)) {
    for (N in names(files[[setting_name]])) {
      res_list <- load_res_list(file.path(cache_dir, files[[setting_name]][[N]]))
      for (i in seq_along(Ts)) {
        mat <- res_list[[i]]
        top3 <- mean(apply(mat, 1, function(r) r[5] %in% r[6:8]))
        nonbottom3 <- mean(1 - apply(mat, 1, function(r) r[5] %in% r[11:13]))
        rows[[idx]] <- data.frame(
          Scenario = setting_name,
          N = as.integer(N),
          T = Ts[i],
          Top3 = round(top3, 2),
          NonBottom3 = round(nonbottom3, 2),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  table4 <- do.call(rbind, rows)
} else {
  n_rep <- get_int_arg(args, 2, 100L)
  seed <- get_int_arg(args, 3, 1234L)
  set.seed(seed)
  table4 <- build_cv_table(
    N_values = c(100L, 200L),
    T_values = c(250L, 500L, 750L, 1000L),
    rho_values = list(0, 0.5),
    covariance_values = list("diagonal", "offdiag"),
    setting_labels = list("Diagonal Sigma_C", "Off-diagonal Sigma_C (rho = 0.5)"),
    N_rep = n_rep
  )
}

out_file <- write_table_csv(table4, repo_root, "Table4_sim_cv.csv")
report_saved_files(out_file)

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
n_rep <- get_int_arg(args, 1, 500L)
seed <- get_int_arg(args, 2, 1234L)

set.seed(seed)
dist_mat <- simulate_factor_distribution(N = 100L, Times = 500L, rho = 0.6, gamma = 0, N_rep = n_rep)

p_qq <- plot_qq(dist_mat[1, ], "Normal Q-Q Plot", "Sample Quantiles")
p_hist <- plot_hist_standardized(dist_mat[1, ], "Histogram of Factor Errors", "Estimation Error")

fig_dir <- get_output_dir(repo_root, "figs")
file_qq <- file.path(fig_dir, "sim_inf_N100_T500.png")
file_hist <- file.path(fig_dir, "sim_inf_N100_T500_hist.png")

ggplot2::ggsave(file_qq, p_qq, width = 6, height = 6, dpi = 300)
ggplot2::ggsave(file_hist, p_hist, width = 6, height = 6, dpi = 300)

report_saved_files(c(file_qq, file_hist))

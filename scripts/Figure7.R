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
source(file.path(repo_root, "R", "fred_core.R"))

args <- commandArgs(trailingOnly = TRUE)
seed <- get_int_arg(args, 1, 20260315L)
set.seed(seed)

data_dir <- file.path(repo_root, "original", "codes")
res <- build_interpretive_outputs(data_dir = data_dir)

gamma_summary <- data.frame(
  dataset = c("FRED-MD", "FRED-QD"),
  gamma_mean = c(mean(res$md_gamma$gamma), mean(res$qd_gamma$gamma)),
  gamma_min = c(min(res$md_gamma$gamma), min(res$qd_gamma$gamma)),
  gamma_max = c(max(res$md_gamma$gamma), max(res$qd_gamma$gamma)),
  share_below_one = c(mean(res$md_gamma$gamma < 0.999), mean(res$qd_gamma$gamma < 0.999))
)

fig_dir <- get_output_dir(repo_root, "figs")
table_dir <- get_output_dir(repo_root, "tables")

fig_file <- file.path(fig_dir, "fig_interpretive_addon.png")
plot_interpretive_figure(
  md_gamma = res$md_gamma,
  qd_gamma = res$qd_gamma,
  md_dates = res$md_dates,
  md_series = res$md_series,
  qd_dates = res$qd_dates,
  qd_series = res$qd_series,
  out_file = fig_file
)

md_gamma_file <- file.path(table_dir, "Figure7_FRED_MD_gamma_path.csv")
qd_gamma_file <- file.path(table_dir, "Figure7_FRED_QD_gamma_path.csv")
gamma_summary_file <- file.path(table_dir, "Figure7_gamma_summary.csv")
corr_file <- file.path(table_dir, "Figure7_factor_macro_correlations.csv")

utils::write.csv(res$md_gamma, md_gamma_file, row.names = FALSE)
utils::write.csv(res$qd_gamma, qd_gamma_file, row.names = FALSE)
utils::write.csv(gamma_summary, gamma_summary_file, row.names = FALSE)
utils::write.csv(res$corr_table, corr_file, row.names = FALSE)

report_saved_files(c(fig_file, md_gamma_file, qd_gamma_file, gamma_summary_file, corr_file))

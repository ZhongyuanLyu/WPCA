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

args <- commandArgs(trailingOnly = TRUE)
n_rep <- get_int_arg(args, 1, 300L)
seed <- get_int_arg(args, 2, 20260315L)

source_script <- file.path(repo_root, "original", "revision_figure2", "run_uv_v4_spd_noise.R")
rscript_path <- get_rscript_path()
work_dir <- ensure_dir(file.path(repo_root, "outputs", "tmp", "figure2_work"))
work_script <- file.path(work_dir, basename(source_script))

file.copy(source_script, work_script, overwrite = TRUE)

cmd_out <- system2(
  rscript_path,
  c(normalizePath(work_script, winslash = "/"), as.character(n_rep), as.character(seed)),
  stdout = TRUE,
  stderr = TRUE
)
status <- attr(cmd_out, "status")
if (!is.null(status) && status != 0) {
  cat(cmd_out, sep = "\n")
  stop("Figure2.R failed while running the preserved revision script.")
}

fig_dir <- get_output_dir(repo_root, "figs")
table_dir <- get_output_dir(repo_root, "tables")

copy_map <- c(
  fig_error_lines.png = file.path(fig_dir, "fig_error_lines.png"),
  sim_raw_u_error.csv = file.path(table_dir, "Figure2_sim_raw_u_error.csv"),
  sim_summary_u_error.csv = file.path(table_dir, "Figure2_sim_summary_u_error.csv"),
  sim_raw_v_error.csv = file.path(table_dir, "Figure2_sim_raw_v_error.csv"),
  sim_summary_v_error.csv = file.path(table_dir, "Figure2_sim_summary_v_error.csv")
)

for (source_name in names(copy_map)) {
  source_file <- file.path(work_dir, source_name)
  if (!file.exists(source_file)) {
    stop("Expected output not found: ", source_file)
  }
  file.copy(source_file, copy_map[[source_name]], overwrite = TRUE)
}

report_saved_files(unname(copy_map))

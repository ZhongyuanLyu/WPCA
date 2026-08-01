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
  stop("Cannot locate the repository root.")
}

repo_root <- locate_repo_root()
setwd(repo_root)

source(file.path(repo_root, "R", "helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1]] else "paper"

rscript_path <- get_rscript_path()
script_runs <- switch(
  mode,
  paper = list(
    list(name = "Figure1.R", args = character()),
    list(name = "Figure2.R", args = c("300", "20260315")),
    list(name = "Figure3.R", args = character()),
    list(name = "Figure4.R", args = character()),
    list(name = "Figure5.R", args = character()),
    list(name = "Figure6.R", args = character()),
    list(name = "Figure7.R", args = c("20260315")),
    list(name = "Table2.R", args = c("100", "1234")),
    list(name = "Table3.R", args = c("100", "1234")),
    list(name = "Table4.R", args = c("recompute", "100", "1234")),
    list(name = "Table5.R", args = c("recompute", "100", "1234"))
  ),
  quick = list(
    list(name = "Figure1.R", args = c("75", "123")),
    list(name = "Figure2.R", args = c("20", "20260315")),
    list(name = "Figure3.R", args = c("75", "123")),
    list(name = "Figure4.R", args = c("75", "123")),
    list(name = "Figure5.R", args = c("75", "1234")),
    list(name = "Figure6.R", args = c("75", "1234")),
    list(name = "Figure7.R", args = c("20260315")),
    list(name = "Table2.R", args = c("2", "1234")),
    list(name = "Table3.R", args = c("2", "1234")),
    list(name = "Table4.R", args = c("recompute", "2", "1234")),
    list(name = "Table5.R", args = c("recompute", "2", "1234"))
  ),
  stop("Unknown mode: ", mode, ". Use 'paper' or 'quick'.")
)

cat("run_all_paper_outputs.R mode:", mode, "\n")
if (mode == "quick") {
  cat("quick mode uses fewer replications and overwrites the same output files. Rerun paper mode to regenerate the outputs with the full simulation settings.\n")
}

for (script_run in script_runs) {
  script_name <- script_run$name
  script_file <- file.path("scripts", script_name)
  run_args <- c(script_file, script_run$args)
  if (length(script_run$args)) {
    cat("Running ", script_name, " with args: ", paste(script_run$args, collapse = " "), "...\n", sep = "")
  } else {
    cat("Running ", script_name, "...\n", sep = "")
  }
  out <- system2(rscript_path, run_args, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    cat(out, sep = "\n")
    stop("run_all_paper_outputs.R failed at ", script_name)
  }
}

cat("All paper output scripts completed in", mode, "mode.\n")

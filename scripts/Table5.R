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
mode <- if (length(args) && args[1] %in% c("cached", "recompute")) args[1] else "cached"

if (mode == "cached") {
  cache_dir <- file.path(repo_root, "original", "codes")
  datasets <- c("FRED-MD" = "", "FRED-QD" = "qd_")
  rows <- list()
  idx <- 1L
  for (dataset_name in names(datasets)) {
    prefix <- datasets[[dataset_name]]
    for (rank_id in 1:3) {
      for (ratio_key in c("0_7", "0_8", "0_9")) {
        env <- new.env(parent = emptyenv())
        load(file.path(cache_dir, paste0("error_", prefix, "r", rank_id, "_tr_", ratio_key, ".RData")), envir = env)
        errs <- get("errs", envir = env)
        rows[[idx]] <- data.frame(
          Dataset = dataset_name,
          Rank = rank_id,
          q_tr = sub("_", ".", ratio_key),
          Method = c("AdaWPCA", "PCA", "HeteroPCA"),
          MeanError = round(rowMeans(errs), 3),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  table5 <- do.call(rbind, rows)
} else {
  n_rep <- get_int_arg(args, 2, 100L)
  seed <- get_int_arg(args, 3, 1234L)
  set.seed(seed)
  table5 <- build_fred_reconstruction_table(
    data_dir = file.path(repo_root, "original", "codes"),
    ranks = 1:3,
    train_ratios = c(0.7, 0.8, 0.9),
    n_rep = n_rep
  )
}

out_file <- write_table_csv(table5, repo_root, "Table5_fred_reconstruction.csv")
report_saved_files(out_file)

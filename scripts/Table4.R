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

source(file.path(repo_root, "R", "helpers.R"))
source(file.path(repo_root, "R", "functions.R"))
source(file.path(repo_root, "R", "simulation_core.R"))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) && args[1] %in% c("cached", "recompute")) args[1] else "recompute"

if (mode == "cached") {
  cache_dir <- file.path(repo_root, "data", "cache")
  manifest_file <- file.path(cache_dir, "table_cache_manifest.csv")
  if (!file.exists(manifest_file)) {
    stop("Verified table-cache manifest not found: ", manifest_file)
  }

  manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
  manifest_columns <- c("Table", "CacheVersion", "ImplementationVersion", "File", "MD5", "Rows", "Columns")
  if (!identical(names(manifest), manifest_columns)) {
    stop("Invalid table-cache manifest schema. Expected: ", paste(manifest_columns, collapse = ", "))
  }

  entry <- manifest[manifest$Table == "Table4", , drop = FALSE]
  expected_cache_version <- "1"
  expected_implementation <- "table4-paper-cache"
  expected_columns <- c("Scenario", "N", "T", "Top3", "NonBottom3")
  if (nrow(entry) != 1L) {
    stop("The table-cache manifest must contain exactly one Table4 entry.")
  }
  if (!identical(as.character(entry$CacheVersion), expected_cache_version) ||
      !identical(entry$ImplementationVersion, expected_implementation)) {
    stop("Table4 cache implementation version mismatch; run Table4.R in recompute mode.")
  }
  if (is.na(entry$MD5) || !nzchar(entry$MD5) || identical(entry$MD5, "PENDING")) {
    stop("Table4 verified cache is not populated; run Table4.R in recompute mode.")
  }
  if (!identical(basename(entry$File), entry$File)) {
    stop("Invalid Table4 cache filename in manifest.")
  }

  cache_file <- file.path(cache_dir, entry$File)
  if (!file.exists(cache_file)) {
    stop("Verified Table4 cache not found: ", cache_file,
         ". Run `Rscript scripts/Table4.R recompute 100 1234`.")
  }
  actual_md5 <- unname(tools::md5sum(cache_file))
  if (is.na(actual_md5) || !identical(tolower(actual_md5), tolower(entry$MD5))) {
    stop("Table4 cache MD5 mismatch; refusing to use an unexpected or modified cache.")
  }

  table4 <- utils::read.csv(cache_file, stringsAsFactors = FALSE, check.names = FALSE)
  manifest_schema <- strsplit(entry$Columns, "\\|", fixed = FALSE)[[1]]
  if (!identical(manifest_schema, expected_columns) || !identical(names(table4), expected_columns)) {
    stop("Table4 cache schema mismatch.")
  }
  if (nrow(table4) != as.integer(entry$Rows) || nrow(table4) != 16L) {
    stop("Table4 cache row-count mismatch.")
  }

  expected_scenarios <- c("Diagonal Sigma_C", "Off-diagonal Sigma_C (rho = 0.5)")
  expected_keys <- as.vector(outer(
    expected_scenarios,
    as.vector(outer(c(100L, 200L), c(250L, 500L, 750L, 1000L), paste, sep = "|")),
    paste,
    sep = "|"
  ))
  actual_keys <- paste(table4$Scenario, table4$N, table4$T, sep = "|")
  if (anyDuplicated(actual_keys) || !setequal(actual_keys, expected_keys)) {
    stop("Table4 cache does not contain the expected scenario/N/T combinations.")
  }
  if (!is.numeric(table4$Top3) || !is.numeric(table4$NonBottom3) ||
      any(!is.finite(table4$Top3)) || any(!is.finite(table4$NonBottom3)) ||
      any(table4$Top3 < 0 | table4$Top3 > 1) ||
      any(table4$NonBottom3 < 0 | table4$NonBottom3 > 1)) {
    stop("Table4 cache contains invalid probability values.")
  }
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

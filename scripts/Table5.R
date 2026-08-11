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
data_dir <- file.path(repo_root, "data")

source(file.path(repo_root, "R", "helpers.R"))
source(file.path(repo_root, "R", "functions.R"))
source(file.path(repo_root, "R", "fred_core.R"))

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

  entry <- manifest[manifest$Table == "Table5", , drop = FALSE]
  expected_cache_version <- "1"
  expected_implementation <- "table5-paper-cache"
  expected_columns <- c("Dataset", "Rank", "q_tr", "Method", "MeanError")
  if (nrow(entry) != 1L) {
    stop("The table-cache manifest must contain exactly one Table5 entry.")
  }
  if (!identical(as.character(entry$CacheVersion), expected_cache_version) ||
      !identical(entry$ImplementationVersion, expected_implementation)) {
    stop("Table5 cache implementation version mismatch; run Table5.R in recompute mode.")
  }
  if (is.na(entry$MD5) || !nzchar(entry$MD5) || identical(entry$MD5, "PENDING")) {
    stop("Table5 verified cache is not populated; ",
         "run `Rscript scripts/Table5.R recompute 100 1234`.")
  }
  if (!identical(basename(entry$File), entry$File)) {
    stop("Invalid Table5 cache filename in manifest.")
  }

  cache_file <- file.path(cache_dir, entry$File)
  if (!file.exists(cache_file)) {
    stop("Verified Table5 cache not found: ", cache_file,
         ". Run `Rscript scripts/Table5.R recompute 100 1234`.")
  }
  actual_md5 <- unname(tools::md5sum(cache_file))
  if (is.na(actual_md5) || !identical(tolower(actual_md5), tolower(entry$MD5))) {
    stop("Table5 cache MD5 mismatch; refusing to use an unexpected or modified cache.")
  }

  table5 <- utils::read.csv(cache_file, stringsAsFactors = FALSE, check.names = FALSE)
  manifest_schema <- strsplit(entry$Columns, "\\|", fixed = FALSE)[[1]]
  if (!identical(manifest_schema, expected_columns) || !identical(names(table5), expected_columns)) {
    stop("Table5 cache schema mismatch.")
  }
  if (nrow(table5) != as.integer(entry$Rows) || nrow(table5) != 54L) {
    stop("Table5 cache row-count mismatch.")
  }

  table5$q_tr <- suppressWarnings(as.numeric(table5$q_tr))
  expected_keys <- as.vector(outer(
    as.vector(outer(c("FRED-MD", "FRED-QD"), 1:3, paste, sep = "|")),
    as.vector(outer(c(0.7, 0.8, 0.9), c("AdaWPCA", "PCA", "HeteroPCA"), paste, sep = "|")),
    paste,
    sep = "|"
  ))
  actual_keys <- paste(table5$Dataset, table5$Rank, table5$q_tr, table5$Method, sep = "|")
  if (anyDuplicated(actual_keys) || !setequal(actual_keys, expected_keys)) {
    stop("Table5 cache does not contain the expected dataset/rank/q_tr/method combinations.")
  }
  if (!is.numeric(table5$MeanError) || any(!is.finite(table5$MeanError)) ||
      any(table5$MeanError < 0)) {
    stop("Table5 cache contains invalid reconstruction errors.")
  }
  table5$q_tr <- sprintf("%.1f", table5$q_tr)
} else {
  n_rep <- get_int_arg(args, 2, 100L)
  seed <- get_int_arg(args, 3, 1234L)
  set.seed(seed)
  table5 <- build_fred_reconstruction_table(
    data_dir = data_dir,
    ranks = 1:3,
    train_ratios = c(0.7, 0.8, 0.9),
    n_rep = n_rep
  )
}

rank_diagnostic <- build_fred_rank_diagnostic(data_dir)
if (!identical(rank_diagnostic$Dataset, c("FRED-MD", "FRED-QD")) ||
    !identical(rank_diagnostic$EstimatedRank, c(1L, 1L))) {
  stop("FRED rank diagnostic differs from the documented reference result.")
}

out_file <- write_table_csv(table5, repo_root, "Table5_fred_reconstruction.csv")
rank_file <- write_table_csv(
  rank_diagnostic, repo_root, "Table5_fred_rank_diagnostic.csv"
)
report_saved_files(c(out_file, rank_file))

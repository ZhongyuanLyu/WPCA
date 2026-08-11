get_int_arg <- function(args, index, default) {
  if (length(args) < index) {
    return(as.integer(default))
  }
  value <- suppressWarnings(as.integer(args[[index]]))
  if (is.na(value)) {
    return(as.integer(default))
  }
  value
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

get_output_dir <- function(repo_root, kind = c("figs", "tables")) {
  kind <- match.arg(kind)
  ensure_dir(file.path(repo_root, "outputs", kind))
}

write_table_csv <- function(x, repo_root, filename) {
  out_file <- file.path(get_output_dir(repo_root, "tables"), filename)
  utils::write.csv(x, out_file, row.names = FALSE)
  out_file
}

report_saved_files <- function(paths) {
  cat("Saved files:\n")
  for (path in paths) {
    cat(" -", normalizePath(path, winslash = "/", mustWork = FALSE), "\n")
  }
}

get_rscript_path <- function() {
  candidates <- c(
    file.path(R.home("bin"), "Rscript"),
    file.path(R.home("bin"), "Rscript.exe")
  )
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) {
    stop("Unable to locate Rscript in R.home('bin').")
  }
  normalizePath(existing[[1]], winslash = "/", mustWork = TRUE)
}

format_mean_sd <- function(x, digits = 3L) {
  fmt <- paste0("%.", digits, "f")
  paste0(sprintf(fmt, mean(x)), "(", sprintf(fmt, stats::sd(x)), ")")
}

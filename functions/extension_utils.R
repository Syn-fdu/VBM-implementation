# ============================================================
# Version 10.2 active extension utilities
# ============================================================

source_project_files <- function(project_dir, files, env = parent.frame()) {
  for (relative_path in files) {
    path <- file.path(project_dir, relative_path)
    if (!file.exists(path)) {
      stop(paste0("Required source file missing: ", path))
    }
    source(path, local = env, chdir = FALSE)
  }
  invisible(TRUE)
}

clear_output_files <- function(...) {
  dirs <- unlist(list(...), use.names = FALSE)
  for (d in dirs) {
    if (!dir.exists(d)) next
    files <- list.files(
      d,
      pattern = "\\.(csv|png|svg|pdf|jpg|jpeg|webp|txt|json)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (length(files) > 0) unlink(files, force = TRUE)
  }
  invisible(TRUE)
}

prepare_extension_dirs <- function(project_dir, extension_id) {
  base_dir <- file.path(project_dir, "output", "extension_results", extension_id)
  dirs <- list(
    base = base_dir,
    tables = file.path(base_dir, "tables"),
    figures = file.path(base_dir, "figures"),
    enhanced = file.path(base_dir, "figures_enhanced")
  )
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

run_python_plot <- function(project_dir, script_relative_path, input_dir, output_dir, label = "enhanced plot") {
  script_path <- file.path(project_dir, script_relative_path)

  if (!file.exists(script_path)) {
    cat("Skipping ", label, " because script is missing: ", script_path, "\n", sep = "")
    return(invisible(FALSE))
  }

  python_bin <- Sys.which("python")
  if (!nzchar(python_bin)) python_bin <- Sys.which("python3")

  if (!nzchar(python_bin)) {
    cat("Skipping ", label, " because neither python nor python3 was found.\n", sep = "")
    return(invisible(FALSE))
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\nLaunching ", label, "...\n", sep = "")
  cat("Python: ", python_bin, "\n", sep = "")
  cat("Script: ", script_path, "\n", sep = "")

  result <- tryCatch(
    system2(
      command = python_bin,
      args = c(
        normalizePath(script_path, winslash = "/", mustWork = FALSE),
        "--input",
        normalizePath(input_dir, winslash = "/", mustWork = FALSE),
        "--output",
        normalizePath(output_dir, winslash = "/", mustWork = FALSE)
      ),
      stdout = TRUE,
      stderr = TRUE
    ),
    warning = function(w) {
      cat(label, " warning: ", conditionMessage(w), "\n", sep = "")
      return(NULL)
    },
    error = function(e) {
      cat(label, " error: ", conditionMessage(e), "\n", sep = "")
      return(NULL)
    }
  )

  if (!is.null(result) && length(result) > 0) {
    cat(paste(result, collapse = "\n"), "\n", sep = "")
  }

  status <- attr(result, "status")
  if (!is.null(status) && status != 0) {
    cat(label, " exited with status ", status, ".\n", sep = "")
    return(invisible(FALSE))
  }

  cat(label, " saved under: ", output_dir, "\n", sep = "")
  invisible(TRUE)
}

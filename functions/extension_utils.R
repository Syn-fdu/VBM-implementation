# ============================================================
# Version 10 extension utilities
# ============================================================
# Shared helpers for the four refactored extension entry points.
# These helpers are deliberately dependency-light so every extension can be
# run either from main.R or directly from its own folder.

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

v10_project_root <- function(anchor_file = NULL) {
  is_root <- function(path) {
    dir.exists(file.path(path, "functions")) &&
      file.exists(file.path(path, "functions", "config.R")) &&
      file.exists(file.path(path, "main.R"))
  }

  candidates <- character(0)

  if (!is.null(anchor_file) && length(anchor_file) > 0 && nzchar(anchor_file[1])) {
    candidates <- c(candidates, dirname(normalizePath(anchor_file[1], winslash = "/", mustWork = FALSE)))
  }

  for (i in seq_len(sys.nframe())) {
    candidate <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(candidate) && length(candidate) > 0 && nzchar(candidate[1])) {
      candidates <- c(candidates, dirname(normalizePath(candidate[1], winslash = "/", mustWork = FALSE)))
    }
  }

  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    candidates <- c(candidates, dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)))
  }

  candidates <- c(candidates, normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  candidates <- unique(candidates[nzchar(candidates)])

  for (start in candidates) {
    current <- start
    for (i in seq_len(12)) {
      if (is_root(current)) return(normalizePath(current, winslash = "/", mustWork = FALSE))
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }

  for (start in candidates) {
    if (!dir.exists(start)) next
    child_dirs <- list.dirs(start, recursive = TRUE, full.names = TRUE)
    child_dirs <- child_dirs[grepl("final_project|version10|version_10|v10", basename(child_dirs), ignore.case = TRUE)]
    for (child in unique(child_dirs)) {
      if (is_root(child)) return(normalizePath(child, winslash = "/", mustWork = FALSE))
    }
  }

  stop("Cannot determine project root. Run from the project folder or source a script by full path.")
}

v10_source_files <- function(project_dir, files, env = parent.frame()) {
  for (relative_path in files) {
    path <- file.path(project_dir, relative_path)
    if (!file.exists(path)) {
      stop(paste0("Required source file missing: ", path))
    }
    source(path, local = env, chdir = FALSE)
  }
  invisible(TRUE)
}

v10_clear_output_files <- function(...) {
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

v10_prepare_extension_dirs <- function(project_dir, extension_id) {
  base_dir <- file.path(project_dir, "output", "version10_extension_results", extension_id)
  dirs <- list(
    base = base_dir,
    tables = file.path(base_dir, "tables"),
    figures = file.path(base_dir, "figures"),
    enhanced = file.path(base_dir, "figures_enhanced")
  )
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

v10_flag <- function(config, name, default = TRUE) {
  opt_name <- paste0("v10_", name)
  opt <- getOption(opt_name, NULL)
  if (!is.null(opt)) return(isTRUE(opt))
  if (!is.null(config) && !is.null(config[[name]])) return(isTRUE(config[[name]]))
  isTRUE(default)
}

v10_run_python_plot <- function(project_dir, script_relative_path, input_dir, output_dir, label = "enhanced plot") {
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

v10_write_extension_index <- function(project_dir, rows) {
  out_dir <- file.path(project_dir, "output", "version10_extension_results")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, "version10_extension_index.csv")
  utils::write.csv(do.call(rbind, rows), out_path, row.names = FALSE)
  cat("Version10 extension index saved to: ", out_path, "\n", sep = "")
  invisible(out_path)
}

v10_run_extension_script <- function(project_dir, relative_path, label) {
  script_path <- file.path(project_dir, relative_path)
  if (!file.exists(script_path)) stop(paste0("Extension script not found: ", script_path))
  cat("\n========== Running Version10 extension: ", label, " ==========\n", sep = "")
  ext_env <- new.env(parent = globalenv())
  ext_env$.v10_project_dir_from_runner <- project_dir
  source(script_path, local = ext_env, chdir = TRUE)
  invisible(TRUE)
}

run_version10_extensions <- function(project_dir = getwd(), config = NULL) {
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

  if (!v10_flag(config, "run_version10_extensions", TRUE)) {
    cat("Version10 extensions are disabled by config$run_version10_extensions.\n")
    return(invisible(FALSE))
  }

  extension_plan <- list(
    list(
      id = "01_hidden_strength_vbm_msm",
      label = "Hidden-strength VBM/MSM comparison",
      flag = "run_extension_01_hidden_strength_vbm_msm",
      path = "extension_01_hidden_strength_vbm_msm.R"
    ),
    list(
      id = "02_vbm_ps_misspecification",
      label = "VBM under PS model misspecification",
      flag = "run_extension_02_vbm_ps_misspecification",
      path = "extension_02_vbm_ps_misspecification.R"
    ),
    list(
      id = "03_good_overlap_instability",
      label = "Good-overlap denominator instability",
      flag = "run_extension_03_good_overlap_instability",
      path = "extension_03_good_overlap_instability.R"
    ),
    list(
      id = "04_rgm_vbm_msm",
      label = "RGM model and VBM/MSM comparison",
      flag = "run_extension_04_rgm_vbm_msm",
      path = "extension_04_rgm_vbm_msm.R"
    )
  )

  rows <- list()
  for (x in extension_plan) {
    enabled <- v10_flag(config, x$flag, TRUE)
    rows[[length(rows) + 1L]] <- data.frame(
      extension_id = x$id,
      extension_label = x$label,
      enabled_by_default = TRUE,
      enabled_this_run = enabled,
      entry_point = x$path,
      output_dir = file.path("output", "version10_extension_results", x$id),
      stringsAsFactors = FALSE
    )

    if (enabled) {
      v10_run_extension_script(project_dir = project_dir, relative_path = x$path, label = x$label)
    } else {
      cat("Skipping Version10 extension ", x$id, " because ", x$flag, " is FALSE.\n", sep = "")
    }
  }

  v10_write_extension_index(project_dir, rows)
  invisible(TRUE)
}

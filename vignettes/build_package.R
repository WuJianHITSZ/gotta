find_package_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (file.exists(file.path(current, "DESCRIPTION"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find DESCRIPTION by searching upward from: ", start)
    }
    current <- parent
  }
}

build_package <- function(
  pkg_root = find_package_root(getwd()),
  build_vignettes = FALSE,
  dest_dir = "C:/Users/jian/Documents/R"
) {
  pkg_root <- normalizePath(pkg_root, winslash = "/", mustWork = TRUE)
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("The 'devtools' package is required. Install it with install.packages('devtools').")
  }

  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  tarball <- devtools::build(
    path = pkg_root,
    binary = FALSE,
    vignettes = isTRUE(build_vignettes),
    quiet = FALSE
  )

  target <- file.path(dest_dir, basename(tarball))
  if (normalizePath(dirname(tarball), winslash = "/", mustWork = TRUE) !=
      normalizePath(dest_dir, winslash = "/", mustWork = TRUE)) {
    file.copy(tarball, target, overwrite = TRUE)
    tarball <- target
  }

  invisible(tarball)
}

if (identical(environment(), globalenv())) {
  build_package(
    pkg_root = find_package_root(getwd())
  )
}

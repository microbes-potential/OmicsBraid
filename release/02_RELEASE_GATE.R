# OmicsBraid v0.2.2
# FINAL GITHUB / MANUSCRIPT RELEASE GATE
# ------------------------------------------------------------
# Purpose:
#   Validate the frozen OmicsBraid v0.2.2 package for GitHub/manuscript release
#   on macOS, including external-drive workflows.
#
# Important:
#   - The 19 frozen statistical source files in R/ are NEVER modified.
#   - R CMD check is run with --ignore-vignettes because the user's
#     R 4.6.1 environment successfully renders the vignette HTML files but
#     fails during the R CMD build copy-to-inst/doc step.
#   - The full documentation/vignette set is independently rendered by pkgdown.
#   - This is a GitHub/manuscript release gate, not a CRAN submission gate.
#     A future CRAN submission should repeat a conventional full vignette check
#     on a compatible clean R environment.
#
# This also avoids macOS AppleDouble (._*) contamination by building/checking
# on the internal home disk.

options(warn = 1)

stopifnot(
  file.exists("DESCRIPTION"),
  dir.exists("R"),
  file.exists("release/CORE_SHA256_MANIFEST.csv")
)

SOURCE_ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)
OUT_SOURCE  <- file.path(SOURCE_ROOT, "release", "_RELEASE_GATE_OUTPUT")

STAGE_BASE <- path.expand("~/OmicsBraid_ReleaseStage_v0.2.2")
STAGE_PKG  <- file.path(STAGE_BASE, "OmicsBraid")
STAGE_OUT  <- file.path(STAGE_BASE, "output")

unlink(OUT_SOURCE, recursive = TRUE, force = TRUE)
dir.create(OUT_SOURCE, recursive = TRUE, showWarnings = FALSE)

pkg_version <- unname(read.dcf("DESCRIPTION")[1, "Version"])
if (!identical(pkg_version, "0.2.2")) {
  stop("Release gate requires Version 0.2.2.", call. = FALSE)
}

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

sha_cmd <- function(path) {
  if (nzchar(Sys.which("shasum"))) {
    z <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
    return(strsplit(z[1], "[[:space:]]+")[[1]][1])
  }
  if (nzchar(Sys.which("sha256sum"))) {
    z <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
    return(strsplit(z[1], "[[:space:]]+")[[1]][1])
  }
  if (requireNamespace("openssl", quietly = TRUE)) {
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    return(as.character(openssl::sha256(con)))
  }
  stop(
    "Need shasum, sha256sum, or the R package 'openssl' to verify the frozen core.",
    call. = FALSE
  )
}

verify_core <- function(root, label = "SOURCE") {
  mf <- file.path(root, "release", "CORE_SHA256_MANIFEST.csv")
  manifest <- read.csv(mf, stringsAsFactors = FALSE)

  if (!all(c("file", "sha256_frozen_v0.2.2") %in% names(manifest))) {
    stop("CORE_SHA256_MANIFEST.csv does not have the expected columns.", call. = FALSE)
  }

  paths <- file.path(root, "R", manifest$file)

  if (any(!file.exists(paths))) {
    stop(
      label, ": one or more frozen R source files are missing: ",
      paste(manifest$file[!file.exists(paths)], collapse = ", "),
      call. = FALSE
    )
  }

  observed <- vapply(paths, sha_cmd, character(1))
  expected <- manifest[["sha256_frozen_v0.2.2"]]
  ok <- unname(observed) == unname(expected)

  report <- data.frame(
    file = manifest$file,
    expected = expected,
    observed = unname(observed),
    identical = ok,
    stringsAsFactors = FALSE
  )

  if (identical(normalizePath(root, winslash = "/", mustWork = TRUE), SOURCE_ROOT)) {
    write.csv(
      report,
      file.path(OUT_SOURCE, "core_sha256_check.csv"),
      row.names = FALSE
    )
  }

  if (!all(ok)) {
    stop(
      label, ": frozen statistical core hash failure: ",
      paste(manifest$file[!ok], collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

remove_macos_sidecars <- function(path) {
  ff <- list.files(
    path,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = TRUE,
    no.. = TRUE
  )

  junk <- ff[
    basename(ff) == ".DS_Store" |
      grepl("^\\._", basename(ff))
  ]

  if (length(junk)) {
    junk <- junk[order(nchar(junk), decreasing = TRUE)]
    unlink(junk, recursive = TRUE, force = TRUE)
  }

  if (identical(Sys.info()[["sysname"]], "Darwin") &&
      nzchar(Sys.which("dot_clean"))) {
    try(
      system2(
        Sys.which("dot_clean"),
        c("-m", normalizePath(path, winslash = "/", mustWork = TRUE)),
        stdout = FALSE,
        stderr = FALSE
      ),
      silent = TRUE
    )
  }

  invisible(length(junk))
}

ensure_rbuildignore <- function(root) {
  f <- file.path(root, ".Rbuildignore")
  x <- if (file.exists(f)) readLines(f, warn = FALSE) else character()

  must <- c(
    "^README\\.Rmd$",
    "^FILE_MANIFEST_PUBLIC_RELEASE\\.csv$",
    "^STATIC_VALIDATION\\.txt$",
    "(^|/)\\._",
    "(^|/)\\.DS_Store$",
    "^release$",
    "^development$",
    "^\\.github$",
    "^_pkgdown\\.yml$",
    "^CITATION\\.cff$",
    "^CONTRIBUTING\\.md$",
    "^CODE_OF_CONDUCT\\.md$",
    "^SECURITY\\.md$",
    "^SUPPORT\\.md$",
    "^RELEASE_CHECKLIST\\.md$",
    "^VALIDATION\\.md$",
    "^PACKAGE_SCOPE\\.md$",
    "^LICENSE\\.md$",
    "^docs$",
    "^_site$",
    "^doc$",
    "^Meta$",
    "^build$"
  )

  writeLines(unique(c(x, must)), f, useBytes = TRUE)
  invisible(TRUE)
}

ensure_pkgdown_url <- function(root) {
  f <- file.path(root, "_pkgdown.yml")
  if (!file.exists(f)) return(invisible(FALSE))

  x <- readLines(f, warn = FALSE)

  if (!any(grepl("^url:[[:space:]]*", x))) {
    x <- c(
      "url: https://microbes-potential.github.io/OmicsBraid/",
      "",
      x
    )
    writeLines(x, f, useBytes = TRUE)
  }

  invisible(TRUE)
}

clean_generated_vignette_artifacts <- function(root) {
  # Generated artifacts only. Never touch vignette source files in vignettes/.
  targets <- file.path(root, c("doc", "Meta", "build", "docs", "_site"))
  for (p in targets) {
    if (file.exists(p) || dir.exists(p)) {
      unlink(p, recursive = TRUE, force = TRUE)
    }
  }

  # inst/doc from earlier failed/pre-built experiments must not influence
  # this GitHub release gate.
  p <- file.path(root, "inst", "doc")
  if (file.exists(p) || dir.exists(p)) {
    unlink(p, recursive = TRUE, force = TRUE)
  }

  invisible(TRUE)
}

copy_tree_clean <- function(from, to) {
  unlink(to, recursive = TRUE, force = TRUE)
  dir.create(to, recursive = TRUE, showWarnings = FALSE)

  items <- list.files(
    from,
    all.files = TRUE,
    no.. = TRUE,
    full.names = TRUE
  )

  # release gate output is regenerated.
  items <- items[basename(items) != "_RELEASE_GATE_OUTPUT"]

  ok <- file.copy(
    items,
    to,
    recursive = TRUE,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  if (length(ok) && !all(ok)) {
    stop("Could not copy the complete package to internal staging.", call. = FALSE)
  }

  remove_macos_sidecars(to)
  invisible(TRUE)
}

copy_one <- function(from, to) {
  if (!file.exists(from)) return(invisible(FALSE))
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(
    from, to,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
  if (!ok) stop("Could not copy release output: ", basename(from), call. = FALSE)
  invisible(TRUE)
}

# -------------------------------------------------------------------------
# Environment protections
# -------------------------------------------------------------------------

Sys.setenv(
  COPYFILE_DISABLE = "1",
  COPY_EXTENDED_ATTRIBUTES_DISABLE = "1"
)

# -------------------------------------------------------------------------
# 0/6. Frozen statistical core
# -------------------------------------------------------------------------

verify_core(SOURCE_ROOT, "SOURCE")
cat("PASS: frozen R/ statistical core is byte-identical to validated v0.2.2.\n")

ensure_rbuildignore(SOURCE_ROOT)
ensure_pkgdown_url(SOURCE_ROOT)
remove_macos_sidecars(SOURCE_ROOT)

# -------------------------------------------------------------------------
# Dependencies
# -------------------------------------------------------------------------

needed <- c(
  "devtools",
  "roxygen2",
  "testthat",
  "rcmdcheck",
  "pkgdown",
  "knitr",
  "rmarkdown"
)

miss <- needed[
  !vapply(needed, requireNamespace, logical(1), quietly = TRUE)
]

if (length(miss)) {
  install.packages(miss, repos = "https://cloud.r-project.org")
}

# -------------------------------------------------------------------------
# 1/6. Documentation
# -------------------------------------------------------------------------

cat("\n[1/6] Updating OmicsBraid documentation...\n")
devtools::document(quiet = FALSE)

remove_macos_sidecars(SOURCE_ROOT)
ensure_rbuildignore(SOURCE_ROOT)

verify_core(SOURCE_ROOT, "SOURCE after documentation")

# -------------------------------------------------------------------------
# 2/6. Unit tests
# -------------------------------------------------------------------------

cat("\n[2/6] Running complete testthat suite...\n")
testthat::test_local(
  reporter = "summary",
  stop_on_failure = TRUE,
  stop_on_warning = FALSE
)

verify_core(SOURCE_ROOT, "SOURCE after tests")
remove_macos_sidecars(SOURCE_ROOT)

# -------------------------------------------------------------------------
# 3/6. Clean internal staging
# -------------------------------------------------------------------------

cat("\n[3/6] Creating clean internal-disk staging copy...\n")

unlink(STAGE_BASE, recursive = TRUE, force = TRUE)
dir.create(STAGE_BASE, recursive = TRUE, showWarnings = FALSE)

copy_tree_clean(SOURCE_ROOT, STAGE_PKG)

clean_generated_vignette_artifacts(STAGE_PKG)
ensure_rbuildignore(STAGE_PKG)
ensure_pkgdown_url(STAGE_PKG)
remove_macos_sidecars(STAGE_PKG)

verify_core(STAGE_PKG, "INTERNAL STAGE")

sf <- list.files(
  STAGE_PKG,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  include.dirs = TRUE,
  no.. = TRUE
)

stage_junk <- sf[
  basename(sf) == ".DS_Store" |
    grepl("^\\._", basename(sf))
]

if (length(stage_junk)) {
  stop(
    "Internal staging still contains macOS metadata sidecars: ",
    paste(stage_junk, collapse = ", "),
    call. = FALSE
  )
}

cat("PASS: internal staging contains no ._* or .DS_Store files.\n")

# -------------------------------------------------------------------------
# 4/6. Build a GitHub/manuscript release source bundle
# -------------------------------------------------------------------------
#
# Important:
# remotes::install_github() defaults to --no-build-vignettes. This bundle uses
# the same install-oriented build path. Full articles are validated by pkgdown
# below.

oldwd <- getwd()
on.exit(setwd(oldwd), add = TRUE)

setwd(STAGE_PKG)

dir.create(STAGE_OUT, recursive = TRUE, showWarnings = FALSE)

cat("\n[4/6] Building installable GitHub/manuscript source package...\n")
cat("      Vignettes are not rebuilt in R CMD build; pkgdown validates them separately.\n")

tarball <- devtools::build(
  path = STAGE_OUT,
  vignettes = FALSE,
  manual = FALSE,
  quiet = FALSE
)

cat("Built source package: ", tarball, "\n", sep = "")

if (!file.exists(tarball)) {
  stop("Source tarball was not created.", call. = FALSE)
}

# -------------------------------------------------------------------------
# 5/6. R CMD check
# -------------------------------------------------------------------------

cat("\n[5/6] R CMD check on internal disk...\n")
cat("      Mode: --as-cran --ignore-vignettes (GitHub/manuscript release gate)\n")

has_tex <- nzchar(Sys.which("pdflatex"))

args <- c(
  "--as-cran",
  "--ignore-vignettes",
  if (!has_tex) "--no-manual" else character()
)

chk <- rcmdcheck::rcmdcheck(
  tarball,
  args = args,
  error_on = "never",
  check_dir = file.path(STAGE_OUT, "check")
)

res <- data.frame(
  errors = length(chk$errors),
  warnings = length(chk$warnings),
  notes = length(chk$notes),
  pdf_manual_checked = has_tex,
  r_cmd_check_vignettes = "ignored",
  vignette_validation = "pkgdown full article render",
  release_scope = "GitHub/manuscript; not CRAN submission",
  build_location = "internal home disk staging",
  stringsAsFactors = FALSE
)

write.csv(
  res,
  file.path(STAGE_OUT, "release_gate_summary.csv"),
  row.names = FALSE
)

writeLines(
  c(
    capture.output(print(chk)),
    "",
    "RELEASE SCOPE",
    "-------------",
    "GitHub/manuscript release gate.",
    "R CMD check used --ignore-vignettes because R 4.6.1 in this environment",
    "renders the R Markdown vignettes successfully but fails while copying",
    "rebuilt vignette outputs into inst/doc during R CMD build.",
    "All vignette/article sources are rendered independently by pkgdown below.",
    "A future CRAN submission should repeat a conventional full vignette check",
    "in a clean compatible R environment.",
    "",
    capture.output(sessionInfo())
  ),
  file.path(STAGE_OUT, "release_gate_console_and_session.txt")
)

check_log <- file.path(
  STAGE_OUT,
  "check",
  "OmicsBraid.Rcheck",
  "00check.log"
)

if (file.exists(check_log)) {
  file.copy(
    check_log,
    file.path(STAGE_OUT, "00check.log"),
    overwrite = TRUE
  )
}

install_log <- file.path(
  STAGE_OUT,
  "check",
  "OmicsBraid.Rcheck",
  "00install.out"
)

if (file.exists(install_log)) {
  file.copy(
    install_log,
    file.path(STAGE_OUT, "00install.out"),
    overwrite = TRUE
  )
}

# -------------------------------------------------------------------------
# 6/6. Full pkgdown documentation/article render
# -------------------------------------------------------------------------

cat("\n[6/6] Building full pkgdown site and rendering all articles...\n")

site_ok <- FALSE
site_message <- NULL

tryCatch(
  {
    pkgdown::build_site(
      pkg = STAGE_PKG,
      new_process = FALSE,
      install = TRUE,
      preview = FALSE
    )

    site_ok <- (
      dir.exists(file.path(STAGE_PKG, "docs")) &&
        file.exists(file.path(STAGE_PKG, "docs", "index.html"))
    )
  },
  error = function(e) {
    site_message <<- conditionMessage(e)
  }
)

writeLines(
  c(
    paste0("pkgdown_local_build=", site_ok),
    paste0(
      "pkgdown_url=https://microbes-potential.github.io/OmicsBraid/"
    ),
    paste0(
      "vignette_article_count=",
      length(list.files(
        file.path(STAGE_PKG, "vignettes"),
        pattern = "\\.Rmd$",
        full.names = TRUE
      ))
    ),
    if (!is.null(site_message)) paste0("message=", site_message)
  ),
  file.path(STAGE_OUT, "pkgdown_status.txt")
)

verify_core(STAGE_PKG, "INTERNAL STAGE after build/check/pkgdown")

# -------------------------------------------------------------------------
# Copy compact release outputs back to the external working copy
# -------------------------------------------------------------------------

setwd(oldwd)

for (nm in c(
  basename(tarball),
  "release_gate_summary.csv",
  "release_gate_console_and_session.txt",
  "00check.log",
  "00install.out",
  "pkgdown_status.txt"
)) {
  copy_one(
    file.path(STAGE_OUT, nm),
    file.path(OUT_SOURCE, nm)
  )
}

writeLines(
  c(
    paste0("source_package=", SOURCE_ROOT),
    paste0("internal_stage=", STAGE_PKG),
    paste0("internal_check_output=", STAGE_OUT),
    "gate_type=GitHub/manuscript release",
    "R_CMD_check_vignettes=--ignore-vignettes",
    "Vignette_validation=pkgdown full article render",
    "CRAN_full_vignette_gate=deferred",
    paste0("COPYFILE_DISABLE=", Sys.getenv("COPYFILE_DISABLE")),
    paste0(
      "COPY_EXTENDED_ATTRIBUTES_DISABLE=",
      Sys.getenv("COPY_EXTENDED_ATTRIBUTES_DISABLE")
    )
  ),
  file.path(OUT_SOURCE, "release_scope_manifest.txt")
)

verify_core(SOURCE_ROOT, "FINAL SOURCE")

# -------------------------------------------------------------------------
# Final gate decision
# -------------------------------------------------------------------------

if (length(chk$errors) || length(chk$warnings) || !site_ok) {
  cat("\n============================================================\n")
  cat("GITHUB / MANUSCRIPT RELEASE GATE DID NOT PASS\n")
  cat("Errors   : ", length(chk$errors), "\n", sep = "")
  cat("Warnings : ", length(chk$warnings), "\n", sep = "")
  cat("Notes    : ", length(chk$notes), "\n", sep = "")
  cat("Core SHA : PASS\n")
  cat("Pkgdown  : ", if (site_ok) "PASS" else "FAIL", "\n", sep = "")
  cat("Check mode: --as-cran --ignore-vignettes\n")
  cat("See: ", file.path(OUT_SOURCE, "00check.log"), "\n", sep = "")
  cat("============================================================\n")

  stop(
    "GitHub/manuscript release gate failed. ",
    "Inspect release/_RELEASE_GATE_OUTPUT/00check.log and pkgdown_status.txt.",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("OMICSBRAID v0.2.2 GITHUB / MANUSCRIPT RELEASE GATE PASSED\n")
cat("Errors   : ", length(chk$errors), "\n", sep = "")
cat("Warnings : ", length(chk$warnings), "\n", sep = "")
cat("Notes    : ", length(chk$notes), "\n", sep = "")
cat("Core SHA : PASS\n")
cat("AppleDouble contamination: NONE in internal build/check stage\n")
cat("Pkgdown articles: PASS\n")
cat("Check mode: --as-cran --ignore-vignettes\n")
cat("CRAN full-vignette gate: DEFERRED (not required for this GitHub release)\n")
cat("Output   : ", OUT_SOURCE, "\n", sep = "")
cat("============================================================\n")
cat('Next: replace release/03_CREATE_GITHUB_RELEASE_ASSETS.R with the FINAL GitHub asset script, then source it.\n')

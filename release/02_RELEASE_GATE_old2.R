# OmicsBraid v0.2.2 publication release gate -- EXTERNAL DRIVE HOTFIX1
#
# Why this hotfix exists:
# macOS can create AppleDouble "._*" metadata files on non-APFS/external volumes.
# Those files can be created AFTER an initial cleanup while roxygen/vignettes/build
# write new files, causing R CMD check warnings.
#
# This gate:
#   1. verifies the frozen v0.2.2 R/ core in the user's package directory;
#   2. regenerates docs/tests there;
#   3. removes AppleDouble metadata again;
#   4. stages a clean copy on the INTERNAL HOME DISK;
#   5. builds/checks/pkgdown from that clean internal staging copy;
#   6. copies only final release reports/tarball back to release/_RELEASE_GATE_OUTPUT.
#
# It does NOT modify any R/*.R statistical source file.

options(warn = 1)

stopifnot(
  file.exists("DESCRIPTION"),
  dir.exists("R"),
  file.exists("release/CORE_SHA256_MANIFEST.csv")
)

SOURCE_ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)
OUT_SOURCE  <- file.path(SOURCE_ROOT, "release", "_RELEASE_GATE_OUTPUT")

# Internal-disk staging. On a normal macOS setup, ~ is on the internal APFS disk.
STAGE_BASE <- path.expand("~/OmicsBraid_ReleaseStage_v0.2.2")
STAGE_PKG  <- file.path(STAGE_BASE, "OmicsBraid")
STAGE_OUT  <- file.path(STAGE_BASE, "output")

unlink(OUT_SOURCE, recursive = TRUE, force = TRUE)
dir.create(OUT_SOURCE, recursive = TRUE, showWarnings = FALSE)

pkg_version <- unname(read.dcf("DESCRIPTION")[1, "Version"])
if (!identical(pkg_version, "0.2.2"))
  stop("Release gate requires Version 0.2.2.", call. = FALSE)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

sha_cmd <- function(path) {
  if (Sys.which("shasum") != "") {
    z <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
    return(strsplit(z[1], "[[:space:]]+")[[1]][1])
  }
  if (Sys.which("sha256sum") != "") {
    z <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
    return(strsplit(z[1], "[[:space:]]+")[[1]][1])
  }
  if (requireNamespace("openssl", quietly = TRUE)) {
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    return(as.character(openssl::sha256(con)))
  }
  stop("Need shasum, sha256sum, or R package 'openssl' to verify frozen core.",
       call. = FALSE)
}

remove_macos_sidecars <- function(path = ".") {
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
    # deepest paths first
    junk <- junk[order(nchar(junk), decreasing = TRUE)]
    unlink(junk, recursive = TRUE, force = TRUE)
  }

  # dot_clean is useful on macOS external volumes, but is not required.
  if (Sys.info()[["sysname"]] == "Darwin" &&
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

ensure_rbuildignore <- function(path = ".") {
  f <- file.path(path, ".Rbuildignore")
  x <- if (file.exists(f)) readLines(f, warn = FALSE) else character()

  must <- c(
    "^FILE_MANIFEST_PUBLIC_RELEASE\\.csv$",
    "^STATIC_VALIDATION\\.txt$",
    "(^|/)\\._",
    "(^|/)\\.DS_Store$",
    "^release$",
    "^development$",
    "^\\.github$",
    "^_pkgdown\\.yml$",
    "^CITATION\\.cff$",
    "^docs$",
    "^_site$"
  )

  add <- must[!must %in% x]
  if (length(add)) {
    x <- c(x, add)
    writeLines(unique(x), f, useBytes = TRUE)
  }
  invisible(f)
}

verify_core <- function(root, label) {
  manifest <- read.csv(
    file.path(root, "release", "CORE_SHA256_MANIFEST.csv"),
    stringsAsFactors = FALSE
  )

  paths <- file.path(root, "R", manifest$file)
  if (any(!file.exists(paths))) {
    stop(label, ": one or more frozen R files are missing.", call. = FALSE)
  }

  observed <- vapply(paths, sha_cmd, character(1))
  expected <- manifest$sha256_frozen_v0.2.2
  ok_each <- unname(observed) == unname(expected)

  report <- data.frame(
    file = manifest$file,
    expected = expected,
    observed = unname(observed),
    identical = ok_each,
    stringsAsFactors = FALSE
  )

  if (identical(root, SOURCE_ROOT)) {
    write.csv(
      report,
      file.path(OUT_SOURCE, "core_sha256_check.csv"),
      row.names = FALSE
    )
  }

  if (!all(ok_each)) {
    bad <- manifest$file[!ok_each]
    stop(
      label,
      ": FROZEN CORE HASH FAILURE: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
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

  # The output directory is regenerated later; don't stage it.
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
    stop("Could not copy the complete package into internal staging.",
         call. = FALSE)
  }

  remove_macos_sidecars(to)
  invisible(to)
}

copy_one <- function(from, to) {
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(from, to, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE))
    stop("Could not copy release output: ", basename(from), call. = FALSE)
}

# -------------------------------------------------------------------------
# 1) Frozen statistical-core hash gate on the ORIGINAL package.
# -------------------------------------------------------------------------

verify_core(SOURCE_ROOT, "SOURCE")
cat("PASS: frozen R/ statistical core is byte-identical to validated v0.2.2.\n")

# -------------------------------------------------------------------------
# 2) Clean external-drive metadata and ensure build exclusions.
# -------------------------------------------------------------------------

ensure_rbuildignore(SOURCE_ROOT)
remove_macos_sidecars(SOURCE_ROOT)

# Prevent macOS copy/archive helpers from materializing resource-fork sidecars.
Sys.setenv(
  COPYFILE_DISABLE = "1",
  COPY_EXTENDED_ATTRIBUTES_DISABLE = "1"
)

# -------------------------------------------------------------------------
# 3) Dependencies.
# -------------------------------------------------------------------------

needed <- c(
  "devtools", "roxygen2", "testthat",
  "rcmdcheck", "pkgdown", "knitr", "rmarkdown"
)
miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  install.packages(miss, repos = "https://cloud.r-project.org")
}

# -------------------------------------------------------------------------
# 4) Regenerate documentation in the user's actual package folder.
# -------------------------------------------------------------------------

cat("\n[1/6] Updating OmicsBraid documentation...\n")
devtools::document(quiet = FALSE)

# Roxygen writes new files; clean again AFTER documentation.
remove_macos_sidecars(SOURCE_ROOT)
ensure_rbuildignore(SOURCE_ROOT)

# -------------------------------------------------------------------------
# 5) Unit tests in the user's actual package folder.
# -------------------------------------------------------------------------

cat("\n[2/6] Running testthat suite...\n")
testthat::test_local(
  reporter = "summary",
  stop_on_failure = TRUE,
  stop_on_warning = FALSE
)

# Tests/build helpers can also trigger Finder metadata on an external drive.
remove_macos_sidecars(SOURCE_ROOT)

# Frozen core must still be identical after documentation/tests.
verify_core(SOURCE_ROOT, "SOURCE after document/tests")

# -------------------------------------------------------------------------
# 6) Stage a CLEAN package copy on the internal disk.
# -------------------------------------------------------------------------

cat("\n[3/6] Creating clean internal-disk staging copy...\n")
unlink(STAGE_BASE, recursive = TRUE, force = TRUE)
dir.create(STAGE_BASE, recursive = TRUE, showWarnings = FALSE)

copy_tree_clean(SOURCE_ROOT, STAGE_PKG)
ensure_rbuildignore(STAGE_PKG)
remove_macos_sidecars(STAGE_PKG)
verify_core(STAGE_PKG, "INTERNAL STAGE")

# Confirm no AppleDouble files are present in stage.
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
# 7) Build + R CMD check entirely on internal disk.
# -------------------------------------------------------------------------

oldwd <- getwd()
on.exit(setwd(oldwd), add = TRUE)
setwd(STAGE_PKG)

dir.create(STAGE_OUT, recursive = TRUE, showWarnings = FALSE)

cat("\n[4/6] R CMD build on internal disk...\n")
tarball <- devtools::build(
  path = STAGE_OUT,
  vignettes = TRUE,
  manual = FALSE,
  quiet = FALSE
)
cat("Built source package: ", tarball, "\n", sep = "")

cat("\n[5/6] R CMD check on internal disk...\n")
has_tex <- nzchar(Sys.which("pdflatex"))
args <- c("--as-cran", if (!has_tex) "--no-manual" else character())

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
    capture.output(sessionInfo())
  ),
  file.path(STAGE_OUT, "release_gate_console_and_session.txt")
)

# Save 00check.log explicitly.
check_log <- file.path(
  STAGE_OUT, "check", "OmicsBraid.Rcheck", "00check.log"
)
if (file.exists(check_log)) {
  file.copy(
    check_log,
    file.path(STAGE_OUT, "00check.log"),
    overwrite = TRUE
  )
}

# -------------------------------------------------------------------------
# 8) Local pkgdown build in internal staging.
# -------------------------------------------------------------------------

cat("\n[6/6] Building local pkgdown preview on internal disk...\n")
site_ok <- FALSE
site_message <- NULL

tryCatch({
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
}, error = function(e) {
  site_message <<- conditionMessage(e)
})

writeLines(
  c(
    paste0("pkgdown_local_build=", site_ok),
    if (!is.null(site_message)) paste0("message=", site_message)
  ),
  file.path(STAGE_OUT, "pkgdown_status.txt")
)

# Recheck stage core after build/check/site.
verify_core(STAGE_PKG, "INTERNAL STAGE after build/check")

# -------------------------------------------------------------------------
# 9) Copy compact final release outputs back to external package directory.
# -------------------------------------------------------------------------

setwd(oldwd)

for (nm in c(
  basename(tarball),
  "release_gate_summary.csv",
  "release_gate_console_and_session.txt",
  "00check.log",
  "pkgdown_status.txt"
)) {
  src <- file.path(STAGE_OUT, nm)
  if (file.exists(src)) {
    copy_one(src, file.path(OUT_SOURCE, nm))
  }
}

# Also record the staging location for reproducibility/debugging.
writeLines(
  c(
    paste0("source_package=", SOURCE_ROOT),
    paste0("internal_stage=", STAGE_PKG),
    paste0("internal_check_output=", STAGE_OUT),
    paste0("COPYFILE_DISABLE=", Sys.getenv("COPYFILE_DISABLE")),
    paste0(
      "COPY_EXTENDED_ATTRIBUTES_DISABLE=",
      Sys.getenv("COPY_EXTENDED_ATTRIBUTES_DISABLE")
    )
  ),
  file.path(OUT_SOURCE, "external_drive_hotfix_manifest.txt")
)

# Final original-source core check.
verify_core(SOURCE_ROOT, "FINAL SOURCE")

# -------------------------------------------------------------------------
# 10) Gate decision.
# -------------------------------------------------------------------------

if (length(chk$errors) || length(chk$warnings)) {
  cat("\n============================================================\n")
  cat("RELEASE GATE DID NOT PASS\n")
  cat("Errors   : ", length(chk$errors), "\n", sep = "")
  cat("Warnings : ", length(chk$warnings), "\n", sep = "")
  cat("Notes    : ", length(chk$notes), "\n", sep = "")
  cat("Core SHA : PASS\n")
  cat("Check was performed on INTERNAL disk staging.\n")
  cat("See: ", file.path(OUT_SOURCE, "00check.log"), "\n", sep = "")
  cat("============================================================\n")

  stop(
    "Release gate failed because R CMD check still has errors or warnings. ",
    "The external-drive AppleDouble problem has been isolated, so inspect ",
    "release/_RELEASE_GATE_OUTPUT/00check.log for any remaining real package issue.",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("OMICSBRAID v0.2.2 RELEASE GATE PASSED\n")
cat("Errors   : ", length(chk$errors), "\n", sep = "")
cat("Warnings : ", length(chk$warnings), "\n", sep = "")
cat("Notes    : ", length(chk$notes), "\n", sep = "")
cat("Core SHA : PASS\n")
cat("AppleDouble contamination: NONE in internal build/check stage\n")
cat("Pkgdown  : ", if (site_ok) "PASS" else "not built; inspect status", "\n", sep = "")
cat("Output   : ", OUT_SOURCE, "\n", sep = "")
cat("============================================================\n")
cat('Next: source("release/03_CREATE_GITHUB_RELEASE_ASSETS.R")\n')

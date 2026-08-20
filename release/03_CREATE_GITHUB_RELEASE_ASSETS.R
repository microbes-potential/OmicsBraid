# OmicsBraid v0.2.2
# FINAL GITHUB RELEASE ASSET CREATOR -- ZIP HOTFIX2 INTERNAL-DISK
# ------------------------------------------------------------
# Run ONLY after release/02_RELEASE_GATE.R has passed.
#
# This version fixes the macOS/base-utils ZIP failure by using the
# cross-platform R package `zip` with an explicit staging root.

stopifnot(
  file.exists("DESCRIPTION"),
  dir.exists("R"),
  file.exists("release/CORE_SHA256_MANIFEST.csv")
)

ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)

gate <- file.path(
  "release",
  "_RELEASE_GATE_OUTPUT",
  "release_gate_summary.csv"
)

if (!file.exists(gate)) {
  stop("Run release/02_RELEASE_GATE.R first.", call. = FALSE)
}

g <- read.csv(gate, stringsAsFactors = FALSE)

if (g$errors[1] != 0 || g$warnings[1] != 0) {
  stop(
    "Release gate did not pass with 0 errors and 0 warnings.",
    call. = FALSE
  )
}

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

sha_cmd <- function(path) {
  if (nzchar(Sys.which("shasum"))) {
    z <- system2(
      "shasum",
      c("-a", "256", path),
      stdout = TRUE,
      stderr = TRUE
    )
    return(strsplit(z[1], "[[:space:]]+")[[1]][1])
  }

  if (nzchar(Sys.which("sha256sum"))) {
    z <- system2(
      "sha256sum",
      path,
      stdout = TRUE,
      stderr = TRUE
    )
    return(strsplit(z[1], "[[:space:]]+")[[1]][1])
  }

  if (requireNamespace("openssl", quietly = TRUE)) {
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    return(as.character(openssl::sha256(con)))
  }

  stop("Need shasum, sha256sum, or openssl.", call. = FALSE)
}

remove_sidecars <- function(path) {
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

  invisible(TRUE)
}

# -------------------------------------------------------------------------
# 0. Verify frozen statistical core
# -------------------------------------------------------------------------

manifest <- read.csv(
  file.path("release", "CORE_SHA256_MANIFEST.csv"),
  stringsAsFactors = FALSE
)

r_paths <- file.path("R", manifest$file)

if (any(!file.exists(r_paths))) {
  stop("One or more frozen R/ source files are missing.", call. = FALSE)
}

obs <- vapply(r_paths, sha_cmd, character(1))
ok <- unname(obs) == unname(manifest[["sha256_frozen_v0.2.2"]])

if (!all(ok)) {
  stop(
    "Frozen core hash failure before creating release assets: ",
    paste(manifest$file[!ok], collapse = ", "),
    call. = FALSE
  )
}

cat("PASS: frozen R/ statistical core verified before asset creation.\n")

Sys.setenv(
  COPYFILE_DISABLE = "1",
  COPY_EXTENDED_ATTRIBUTES_DISABLE = "1"
)

remove_sidecars(ROOT)

# -------------------------------------------------------------------------
# Dependencies
# -------------------------------------------------------------------------

needed <- c("devtools", "zip")
miss <- needed[
  !vapply(needed, requireNamespace, logical(1), quietly = TRUE)
]

if (length(miss)) {
  install.packages(miss, repos = "https://cloud.r-project.org")
}

# -------------------------------------------------------------------------
# Output directory
# -------------------------------------------------------------------------

out <- file.path("release", "_RELEASE_ASSETS")
unlink(out, recursive = TRUE, force = TRUE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 1/3. Canonical GitHub/manuscript source tarball
# -------------------------------------------------------------------------

cat("\n[1/3] Building GitHub/manuscript source tarball...\n")

tarball <- devtools::build(
  path = out,
  vignettes = FALSE,
  manual = FALSE,
  quiet = FALSE
)

if (!file.exists(tarball)) {
  stop("Source tarball was not created.", call. = FALSE)
}

tarball <- normalizePath(tarball, winslash = "/", mustWork = TRUE)

cat("Created: ", tarball, "\n", sep = "")

# -------------------------------------------------------------------------
# 2/3. Clean GitHub repository ZIP
# -------------------------------------------------------------------------

cat("\n[2/3] Building clean GitHub repository ZIP...\n")

staging <- tempfile("OmicsBraid_repo_")
dir.create(staging, recursive = TRUE, showWarnings = FALSE)

repo_root <- file.path(staging, "OmicsBraid")
dir.create(repo_root, recursive = TRUE, showWarnings = FALSE)

files <- list.files(
  ROOT,
  recursive = TRUE,
  all.files = TRUE,
  full.names = TRUE,
  include.dirs = FALSE,
  no.. = TRUE
)

exclude_rx <- paste(
  c(
    "(^|/)\\.git(/|$)",
    "(^|/)\\.Rproj\\.user(/|$)",
    "(^|/)docs(/|$)",
    "(^|/)doc(/|$)",
    "(^|/)Meta(/|$)",
    "(^|/)build(/|$)",
    "(^|/)inst/doc(/|$)",
    "(^|/)release/_RELEASE_",
    "(^|/)\\.DS_Store$",
    "(^|/)\\._",
    "(^|/)_CHECK_OUTPUT(/|$)",
    "(^|/)_SMOKE_TEST_OUTPUT(/|$)"
  ),
  collapse = "|"
)

# Convert to forward slashes for matching.
files_norm <- gsub("\\\\", "/", files)
keep <- !grepl(exclude_rx, files_norm)
files <- files[keep]

if (!length(files)) {
  stop("No files remained for the GitHub repository ZIP.", call. = FALSE)
}

for (f in files) {
  rel <- substring(
    normalizePath(f, winslash = "/", mustWork = TRUE),
    nchar(ROOT) + 2L
  )

  dst <- file.path(repo_root, rel)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)

  ok_copy <- file.copy(
    f,
    dst,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  if (!ok_copy) {
    stop("Failed to stage repository file: ", rel, call. = FALSE)
  }
}

remove_sidecars(staging)

repo_files <- list.files(
  repo_root,
  recursive = TRUE,
  all.files = TRUE,
  full.names = TRUE,
  include.dirs = FALSE,
  no.. = TRUE
)

if (!length(repo_files)) {
  stop(
    "Repository staging folder exists but contains no files.",
    call. = FALSE
  )
}

cat(
  "Staged repository files: ",
  length(repo_files),
  "\n",
  sep = ""
)

# IMPORTANT: create the ZIP on the INTERNAL Mac temp disk first.
# Some external filesystems mounted under /Volumes can reject direct rzip writes.
zipfile <- normalizePath(
  file.path(out, "OmicsBraid_GitHub_v0.2.2.zip"),
  winslash = "/",
  mustWork = FALSE
)

internal_zip_dir <- file.path(
  path.expand("~"),
  "OmicsBraid_ReleaseStage_v0.2.2",
  "zip_output"
)
dir.create(
  internal_zip_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

internal_zip <- file.path(
  internal_zip_dir,
  "OmicsBraid_GitHub_v0.2.2.zip"
)

if (file.exists(internal_zip)) {
  unlink(internal_zip, force = TRUE)
}

zip::zipr(
  zipfile = internal_zip,
  files = "OmicsBraid",
  root = staging,
  recurse = TRUE,
  include_directories = TRUE
)

if (!file.exists(internal_zip)) {
  stop(
    "GitHub repository ZIP was not created on the internal disk.",
    call. = FALSE
  )
}

if (file.info(internal_zip)$size <= 0) {
  stop(
    "GitHub repository ZIP was created on the internal disk but is empty.",
    call. = FALSE
  )
}

# Audit the ZIP BEFORE copying it back to the external drive.
zip_listing <- zip::zip_list(internal_zip)

if (!nrow(zip_listing)) {
  stop("Internal GitHub repository ZIP contains no entries.", call. = FALSE)
}

if (!any(grepl("(^|/)DESCRIPTION$", zip_listing$filename))) {
  stop(
    "GitHub ZIP audit failed: DESCRIPTION is missing.",
    call. = FALSE
  )
}

if (!any(grepl("(^|/)R/workflow\\.R$", zip_listing$filename))) {
  stop(
    "GitHub ZIP audit failed: frozen R/ source files are missing.",
    call. = FALSE
  )
}

cat(
  "PASS: internal GitHub ZIP audit found ",
  nrow(zip_listing),
  " entries.\n",
  sep = ""
)

# Only after the audit passes, copy the completed ZIP to the external
# release-assets directory.
dir.create(
  dirname(zipfile),
  recursive = TRUE,
  showWarnings = FALSE
)

copied <- file.copy(
  internal_zip,
  zipfile,
  overwrite = TRUE,
  copy.mode = TRUE,
  copy.date = TRUE
)

if (!copied || !file.exists(zipfile)) {
  stop(
    "ZIP was valid on the internal disk but could not be copied to ",
    "release/_RELEASE_ASSETS.",
    call. = FALSE
  )
}

if (file.info(zipfile)$size != file.info(internal_zip)$size) {
  stop(
    "Copied GitHub ZIP size does not match the verified internal ZIP.",
    call. = FALSE
  )
}

cat("Created: ", zipfile, "\n", sep = "")

cat("PASS: GitHub ZIP copied from verified internal build.\n")

# -------------------------------------------------------------------------
# 3/3. Checksums + release notes
# -------------------------------------------------------------------------

cat("\n[3/3] Writing checksums and release notes...\n")

assets <- c(tarball, zipfile)

asset_manifest <- data.frame(
  file = basename(assets),
  sha256 = vapply(assets, sha_cmd, character(1)),
  stringsAsFactors = FALSE
)

write.csv(
  asset_manifest,
  file.path(out, "SHA256SUMS.csv"),
  row.names = FALSE
)

writeLines(
  paste(asset_manifest$sha256, asset_manifest$file),
  file.path(out, "SHA256SUMS.txt")
)

writeLines(
  c(
    "OmicsBraid v0.2.2 - frozen manuscript validation release",
    "",
    "Release scope",
    "-------------",
    "This is the frozen GitHub/manuscript software release used for the",
    "OmicsBraid methods study.",
    "",
    "The 19-file statistical R/ core is SHA256-verified against validated",
    "v0.2.2 before release assets are created.",
    "",
    "Validation represented by this frozen release includes pure simulation,",
    "robust calibration, CI sensitivity, CPTAC real-background semi-synthetic",
    "validation, comparator/ablation benchmarking, CPTAC-GBM application,",
    "independent CPTAC-LUAD validation, and ranked pathway analyses.",
    "",
    "Vignettes / documentation",
    "--------------------------",
    "The GitHub repository retains all vignette source files. Rendered",
    "documentation is published through pkgdown. The source tarball is built",
    "without rebuilding vignettes, matching the GitHub/manuscript release",
    "workflow validated by release/02_RELEASE_GATE.R.",
    "",
    "Repository:",
    "https://github.com/microbes-potential/OmicsBraid",
    "",
    "Documentation:",
    "https://microbes-potential.github.io/OmicsBraid/"
  ),
  file.path(out, "GITHUB_RELEASE_NOTES_v0.2.2.md")
)

writeLines(
  c(
    "OmicsBraid v0.2.2 release asset manifest",
    paste0("Generated: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "Frozen statistical core: PASS",
    "GitHub repository ZIP audit: PASS",
    "Source tarball vignette mode: no rebuild",
    "GitHub repository includes vignette source: YES",
    "Rendered documentation: pkgdown",
    "CRAN full-vignette gate: deferred"
  ),
  file.path(out, "RELEASE_ASSET_SCOPE.txt")
)

# Save ZIP listing for auditability.
write.csv(
  zip_listing,
  file.path(out, "GITHUB_ZIP_CONTENTS.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("OMICSBRAID v0.2.2 RELEASE ASSETS CREATED\n")
cat("Frozen core: PASS\n")
cat("GitHub ZIP audit: PASS\n")
cat("Source tarball: ", basename(tarball), "\n", sep = "")
cat("Repository ZIP: ", basename(zipfile), "\n", sep = "")
cat("Output folder: ", normalizePath(out), "\n", sep = "")
cat("============================================================\n")

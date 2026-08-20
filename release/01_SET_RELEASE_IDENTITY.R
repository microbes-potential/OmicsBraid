# OmicsBraid v0.2.2 release metadata helper
# This copy is pre-configured for the public GitHub owner/maintainer.
# Source from the OmicsBraid package root if you ever need to refresh metadata.

GITHUB_USERNAME <- "microbes-potential"
GIVEN_NAME      <- "Adeel"
FAMILY_NAME     <- "Farooq"
EMAIL           <- "jhwanj9@gmail.com"
ORCID           <- ""
AFFILIATION     <- "University of Guelph"

# ---- do not edit below unless you know why ----
stopifnot(file.exists("DESCRIPTION"), dir.exists("R"))

replace_or_add_dcf <- function(lines, key, value) {
  start <- grep(paste0("^", gsub("/", "\\/", key), ":"), lines)
  if (length(start)) {
    i <- start[1]
    j <- i + 1L
    while (j <= length(lines) && grepl("^[[:space:]]", lines[j])) j <- j + 1L
    lines <- lines[-seq.int(i, j - 1L)]
  }
  c(lines, paste0(key, ": ", value))
}

desc <- readLines("DESCRIPTION", warn = FALSE)
i <- grep("^Authors@R:", desc)
if (length(i)) {
  j <- i[1] + 1L
  while (j <= length(desc) && grepl("^[[:space:]]", desc[j])) j <- j + 1L
  desc <- desc[-seq.int(i[1], j - 1L)]
}
author_line <- paste0(
  'Authors@R: person("', GIVEN_NAME, '", "', FAMILY_NAME,
  '", email = "', EMAIL, '", role = c("aut", "cre")',
  if (nzchar(ORCID)) paste0(', comment = c(ORCID = "', ORCID, '")') else '',
  ')'
)
desc <- c(desc, author_line)
desc <- replace_or_add_dcf(desc, "URL",
  paste0("https://", GITHUB_USERNAME, ".github.io/OmicsBraid/, https://github.com/", GITHUB_USERNAME, "/OmicsBraid"))
desc <- replace_or_add_dcf(desc, "BugReports",
  paste0("https://github.com/", GITHUB_USERNAME, "/OmicsBraid/issues"))
writeLines(desc, "DESCRIPTION", useBytes = TRUE)

cff <- c(
  "cff-version: 1.2.0",
  'message: "If you use OmicsBraid, please cite the software release and the accompanying methods article when available."',
  'title: "OmicsBraid: Covariance-Aware Inference of Cross-Omic Effect Trajectories"',
  "type: software",
  "version: 0.2.2",
  "date-released: 2026-08-20",
  "authors:",
  paste0('  - given-names: "', GIVEN_NAME, '"'),
  paste0('    family-names: "', FAMILY_NAME, '"'),
  paste0('    affiliation: "', AFFILIATION, '"'),
  paste0('    email: "', EMAIL, '"'),
  if (nzchar(ORCID)) paste0('    orcid: "https://orcid.org/', sub('^https://orcid.org/', '', ORCID), '"') else NULL,
  paste0('repository-code: "https://github.com/', GITHUB_USERNAME, '/OmicsBraid"'),
  paste0('url: "https://', GITHUB_USERNAME, '.github.io/OmicsBraid/"'),
  "license: MIT",
  "keywords:", "  - multi-omics", "  - effect size", "  - covariance",
  "  - heterogeneity", "  - trajectory", "  - bioinformatics", "  - R"
)
writeLines(cff, "CITATION.cff", useBytes = TRUE)

orc_note <- if (nzchar(ORCID)) paste0("; ORCID: ", ORCID) else ""
cit <- c(
  'citHeader("To cite OmicsBraid in publications, please use:")', "",
  "bibentry(", '  bibtype = "Manual",',
  '  title = "OmicsBraid: Covariance-Aware Inference of Cross-Omic Effect Trajectories",',
  paste0('  author = c(person("', GIVEN_NAME, '", "', FAMILY_NAME, '")),'),
  '  year = "2026",',
  paste0('  note = "R package version 0.2.2; frozen manuscript analysis release', orc_note, '",'),
  paste0('  url = "https://github.com/', GITHUB_USERNAME, '/OmicsBraid",'),
  "  textVersion = paste0(",
  paste0('    "', FAMILY_NAME, ', ', GIVEN_NAME, ' (2026). OmicsBraid: Covariance-Aware ",'),
  '    "Inference of Cross-Omic Effect Trajectories. R package version 0.2.2."',
  "  )", ")"
)
writeLines(cit, "inst/CITATION", useBytes = TRUE)

cat("PASS: public release identity refreshed for ", GITHUB_USERNAME, ".\n", sep = "")
cat('Next: source("release/02_RELEASE_GATE.R")\n')

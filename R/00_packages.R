# 00_packages.R
# Install (if missing) and load all packages used in Project 1:
# "Statistical Analysis of SPT-N Values and Shear Strength Parameters in Dhaka"

set.seed(123)

# Use a CRAN mirror if none is set (e.g. when run via Rscript)
if (is.null(getOption("repos")) || getOption("repos")[["CRAN"]] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

packages <- c(
  "tidyverse", "janitor", "psych", "moments", "car", "lmtest",
  "corrplot", "GGally", "caret", "rsample", "boot", "broom",
  "knitr", "flextable", "officer", "sandwich", "patchwork", "ggrepel"
)

# Install any packages that are not already installed
missing <- packages[!packages %in% rownames(installed.packages())]
if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, dependencies = TRUE)
}

# Load all packages quietly
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

# Make sure the output folders exist
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables",  recursive = TRUE, showWarnings = FALSE)

# Plain-English summary
loaded <- packages[packages %in% (.packages())]
cat("\nSummary: ", length(loaded), " of ", length(packages),
    " packages are loaded and ready.\n", sep = "")
if (length(loaded) < length(packages)) {
  cat("Not loaded: ", paste(setdiff(packages, loaded), collapse = ", "), "\n")
}

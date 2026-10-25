# 99_run_all.R
# Run the whole analysis (scripts 01-10) in order, each in a fresh environment,
# then write outputs/sessionInfo.txt.
# Run from the project folder:  Rscript R/99_run_all.R

scripts <- c(
  "R/01_clean.R",
  "R/02_descriptive.R",
  "R/03_normality.R",
  "R/04_correlation.R",
  "R/05_simple_models.R",
  "R/06_multiple_regression.R",
  "R/07_diagnostics.R",
  "R/08_validation.R",
  "R/09_published_comparison.R",
  "R/10_study_area_map.R"
)

stopifnot(all(file.exists(scripts)))

timings <- numeric(0)
for (s in scripts) {
  cat("\n\n==================== ", s, " ====================\n", sep = "")
  t0 <- Sys.time()
  # A fresh environment per script so no object leaks from one script to the next
  source(s, local = new.env(), echo = FALSE)
  timings[s] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

# ---- Session info ----------------------------------------------------------
writeLines(c(paste("Run completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
             "",
             capture.output(sessionInfo())),
           "outputs/sessionInfo.txt")

# ---- Plain-English summary -------------------------------------------------
cat("\n\n==================== Run summary ====================\n")
print(data.frame(script = names(timings), seconds = round(timings, 1)), row.names = FALSE)
cat(sprintf("\nAll %d scripts ran without error in %.0f seconds.\n",
            length(scripts), sum(timings)))
cat("Tables:", length(list.files("outputs/tables", pattern = "\\.csv$")),
    "| Figures:", length(list.files("outputs/figures", pattern = "\\.png$")), "\n")
cat("Saved outputs/sessionInfo.txt\n")

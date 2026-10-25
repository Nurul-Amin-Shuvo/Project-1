# 01_clean.R
# Read the raw Dhaka SPT-N / qu data, clean it, run basic quality checks,
# and save the clean data (data/clean.rds) and Table 4.1 (data summary).

source("R/00_packages.R")
set.seed(123)

raw_file <- "data/Dhaka qu, spt N, depth.csv"
expected_n <- 134

# ---- Read ------------------------------------------------------------------
# Read everything as text first so non-numeric entries can be detected
raw <- readr::read_csv(raw_file, col_types = cols(.default = col_character()),
                       name_repair = "unique_quiet", show_col_types = FALSE)

cat("Raw file:", nrow(raw), "rows x", ncol(raw), "columns\n")

# ---- Drop empty columns and rows -------------------------------------------
# Treat blank strings as missing so janitor can spot empty rows/columns
raw <- raw %>%
  mutate(across(everything(), ~ na_if(trimws(.x), ""))) %>%
  janitor::remove_empty(which = c("rows", "cols"), quiet = TRUE)

cat("After dropping empty rows/columns:", nrow(raw), "rows x", ncol(raw),
    "columns\n")

# ---- Rename ----------------------------------------------------------------
# Match columns by name pattern rather than position, so a column shift fails loudly
find_col <- function(pattern) {
  hit <- grep(pattern, names(raw), ignore.case = TRUE)
  if (length(hit) != 1) {
    stop("Expected one column matching '", pattern, "', found ", length(hit),
         ". Columns are: ", paste(names(raw), collapse = " | "))
  }
  names(raw)[hit]
}

clean <- raw %>%
  select(depth_m = all_of(find_col("depth")),
         N       = all_of(find_col("spt")),
         qu_kPa  = all_of(find_col("qu|compressive")))

# ---- Convert to numeric ----------------------------------------------------
non_numeric <- clean %>%
  mutate(row = row_number()) %>%
  pivot_longer(-row, names_to = "variable", values_to = "value") %>%
  filter(!is.na(value), is.na(suppressWarnings(as.numeric(value))))

if (nrow(non_numeric) > 0) {
  cat("\nWARNING: non-numeric entries found (set to NA):\n")
  print(non_numeric)
}

clean <- clean %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

# ---- Quality checks --------------------------------------------------------
n_rows <- nrow(clean)
missing <- colSums(is.na(clean))
n_dup <- sum(duplicated(clean))

cat("\n--- Data checks ---\n")
cat("Rows:", n_rows, "\n")
if (n_rows != expected_n) {
  cat("WARNING: expected", expected_n, "rows but found", n_rows,
      "- check the raw file.\n")
}
cat("Missing values per column:\n")
print(missing)
cat("Fully duplicated rows:", n_dup, "\n")
if (n_dup > 0) {
  cat("Duplicated rows (kept in the data; check whether they are genuine repeats):\n")
  print(clean[duplicated(clean) | duplicated(clean, fromLast = TRUE), ])
}

# ---- Save clean data -------------------------------------------------------
saveRDS(clean, "data/clean.rds")

# ---- Table 4.1: data summary -----------------------------------------------
table4_1 <- clean %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "value") %>%
  mutate(Variable = factor(Variable,
                           levels = c("depth_m", "N", "qu_kPa"),
                           labels = c("Depth (m)", "SPT-N (blows/300 mm)",
                                      "qu (kPa)"))) %>%
  group_by(Variable) %>%
  summarise(n       = sum(!is.na(value)),
            Missing = sum(is.na(value)),
            Mean    = mean(value, na.rm = TRUE),
            SD      = sd(value, na.rm = TRUE),
            Min     = min(value, na.rm = TRUE),
            Median  = median(value, na.rm = TRUE),
            Max     = max(value, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(across(Mean:Max, ~ round(.x, 2)))

write_csv(table4_1, "outputs/tables/Table4_1_data_summary.csv")

# ---- Plain-English summary -------------------------------------------------
cat("\n--- Summary ---\n")
print(table4_1)
cat("\nThe clean dataset has ", n_rows, " samples (",
    ifelse(n_rows == expected_n, "as expected", paste("expected", expected_n)),
    "), ", sum(missing), " missing values and ", n_dup,
    " duplicated rows.\n", sep = "")
with(table4_1, cat(sprintf("%s ranges from %s to %s (mean %s).\n",
                           Variable, Min, Max, Mean), sep = ""))
cat("Saved data/clean.rds and outputs/tables/Table4_1_data_summary.csv\n")

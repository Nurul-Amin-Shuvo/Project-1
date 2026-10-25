# 03_normality.R
# Shapiro-Wilk normality tests (Table 4.3), Q-Q plots (Fig 4.3) and
# outlier screening by the 1.5 IQR rule (Table 4.3b).
# All 134 samples are used, including the 4 duplicated rows.

source("R/00_packages.R")
set.seed(123)

clean <- readRDS("data/clean.rds") %>%
  mutate(sample_id = row_number(), .before = 1)

theme_set(theme_bw(base_size = 11, base_family = "Times New Roman"))
save_fig <- function(plot, file, width, height) {
  ggsave(file.path("outputs/figures", file), plot, width = width,
         height = height, units = "in", dpi = 300, device = ragg::agg_png)
}

# Variables to test (natural log for the transformed ones)
test_vars <- clean %>%
  transmute(`Depth (m)`               = depth_m,
            `SPT-N (blows/300 mm)`    = N,
            `qu (kPa)`                = qu_kPa,
            `ln SPT-N`                = log(N),
            `ln qu`                   = log(qu_kPa))

long <- test_vars %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "value") %>%
  mutate(Variable = factor(Variable, levels = names(test_vars)))

# ---- Table 4.3: Shapiro-Wilk tests -----------------------------------------
table4_3 <- long %>%
  group_by(Variable) %>%
  summarise(n = n(),
            W = unname(shapiro.test(value)$statistic),
            p_value = shapiro.test(value)$p.value,
            .groups = "drop") %>%
  mutate(normal_at_5pct = ifelse(p_value > 0.05, "Yes", "No"),
         W = round(W, 4),
         p_value = signif(p_value, 3))

write_csv(table4_3, "outputs/tables/Table4_3_normality.csv")

# ---- Fig 4.3: Q-Q plots ----------------------------------------------------
fig4_3 <- ggplot(long, aes(sample = value)) +
  stat_qq_line(colour = "black", linewidth = 0.6) +
  stat_qq(shape = 21, fill = "grey60", colour = "black", size = 1.6) +
  facet_wrap(~ Variable, scales = "free_y", nrow = 2) +
  labs(x = "Theoretical quantiles (standard normal)",
       y = "Sample quantiles")

save_fig(fig4_3, "Fig4_3_qq_plots.png", width = 8, height = 5.5)

# ---- Outliers: 1.5 IQR rule ------------------------------------------------
iqr_flag <- function(x) {
  q <- quantile(x, c(0.25, 0.75))
  fence <- 1.5 * diff(q)
  x < q[1] - fence | x > q[2] + fence
}

outliers <- clean %>%
  mutate(out_depth = iqr_flag(depth_m),
         out_N     = iqr_flag(N),
         out_qu    = iqr_flag(qu_kPa)) %>%
  filter(out_depth | out_N | out_qu) %>%
  mutate(outlier_in = pmap_chr(list(out_depth, out_N, out_qu), function(d, n, q) {
    paste(c("Depth", "N", "qu")[c(d, n, q)], collapse = ", ")
  })) %>%
  select(sample_id, depth_m, N, qu_kPa, outlier_in)

fences <- clean %>%
  summarise(across(c(depth_m, N, qu_kPa), list(
    lower = ~ quantile(.x, 0.25) - 1.5 * IQR(.x),
    upper = ~ quantile(.x, 0.75) + 1.5 * IQR(.x))))

write_csv(outliers, "outputs/tables/Table4_3b_outliers.csv")

# ---- Plain-English summary -------------------------------------------------
cat("\n--- Table 4.3: Shapiro-Wilk normality tests ---\n")
print(as.data.frame(table4_3), row.names = FALSE)

cat("\nSummary:\n")
for (i in seq_len(nrow(table4_3))) {
  r <- table4_3[i, ]
  cat(sprintf("- %s: W = %.3f, p = %s -> %s.\n", r$Variable, r$W,
              format.pval(r$p_value, digits = 3, eps = 1e-4),
              ifelse(r$normal_at_5pct == "Yes",
                     "consistent with a normal distribution at 5%",
                     "not normal at 5%")))
}

cat(sprintf("\nOutlier fences (1.5 IQR): Depth %.2f to %.2f m; N %.2f to %.2f; qu %.2f to %.2f kPa.\n",
            fences$depth_m_lower, fences$depth_m_upper,
            fences$N_lower, fences$N_upper,
            fences$qu_kPa_lower, fences$qu_kPa_upper))
cat(nrow(outliers), "of", nrow(clean),
    "samples are outliers in at least one variable",
    sprintf("(Depth: %d, N: %d, qu: %d).", sum(grepl("Depth", outliers$outlier_in)),
            sum(grepl("\\bN\\b", outliers$outlier_in)),
            sum(grepl("qu", outliers$outlier_in))),
    "They are flagged only, not removed.\n\n")
print(as.data.frame(outliers), row.names = FALSE)
cat("\nSaved Table4_3_normality.csv, Table4_3b_outliers.csv, Fig4_3_qq_plots.png\n")

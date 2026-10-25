# 02_descriptive.R
# Descriptive statistics (Table 4.2), histograms with density curves (Fig 4.1),
# boxplots (Fig 4.2) and depth profiles of SPT-N and qu (Fig 4.2b).
# All 134 samples are used, including the 4 duplicated rows.

source("R/00_packages.R")
set.seed(123)

clean <- readRDS("data/clean.rds")

# Display labels with units, in a fixed order
var_labels <- c(depth_m = "Depth (m)",
                N       = "SPT-N (blows/300 mm)",
                qu_kPa  = "qu (kPa)")

# ragg uses system fonts directly, so Times New Roman works on Windows
theme_set(theme_bw(base_size = 11, base_family = "Times New Roman"))
save_fig <- function(plot, file, width, height) {
  ggsave(file.path("outputs/figures", file), plot, width = width,
         height = height, units = "in", dpi = 300, device = ragg::agg_png)
}

long <- clean %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = names(var_labels),
                           labels = var_labels))

# ---- Table 4.2: descriptive statistics -------------------------------------
# Skewness: moments::skewness. Kurtosis: excess kurtosis (normal = 0).
table4_2 <- long %>%
  group_by(Variable = variable) %>%
  summarise(n        = sum(!is.na(value)),
            Mean     = mean(value),
            Median   = median(value),
            SD       = sd(value),
            CV_pct   = 100 * SD / Mean,
            Min      = min(value),
            Max      = max(value),
            Q1       = quantile(value, 0.25),
            Q3       = quantile(value, 0.75),
            Skewness = moments::skewness(value),
            Kurtosis_excess = moments::kurtosis(value) - 3,
            .groups = "drop") %>%
  mutate(across(Mean:Kurtosis_excess, ~ round(.x, 2)))

write_csv(table4_2, "outputs/tables/Table4_2_descriptive.csv")

# ---- Fig 4.1: histograms with density curves -------------------------------
fig4_1 <- ggplot(long, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20,
                 fill = "grey80", colour = "grey30") +
  geom_density(colour = "black", linewidth = 0.7) +
  facet_wrap(~ variable, scales = "free", nrow = 1) +
  labs(x = NULL, y = "Density")

save_fig(fig4_1, "Fig4_1_histograms.png", width = 9, height = 3.2)

# ---- Fig 4.2: boxplots -----------------------------------------------------
fig4_2 <- ggplot(long, aes(x = "", y = value)) +
  geom_boxplot(fill = "grey85", outlier.shape = 21, outlier.fill = "white",
               width = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 2.5) +
  facet_wrap(~ variable, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL,
       caption = "Boxes show Q1 to Q3 with the median; x marks the mean.") +
  theme(axis.ticks.x = element_blank())

save_fig(fig4_2, "Fig4_2_boxplots.png", width = 7, height = 3.5)

# ---- Fig 4.2b: depth profiles (borehole-log style) -------------------------
profile <- clean %>%
  pivot_longer(c(N, qu_kPa), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = c("N", "qu_kPa"),
                           labels = var_labels[c("N", "qu_kPa")]))

# Very few samples lie below 25 m, so the trend line stops there
trend_max_depth <- 25

fig4_2b <- ggplot(profile, aes(x = value, y = depth_m)) +
  geom_point(shape = 21, fill = "grey60", colour = "black", size = 1.8) +
  geom_smooth(data = filter(profile, depth_m <= trend_max_depth),
              method = "loess", formula = y ~ x, orientation = "y",
              se = FALSE, colour = "black", linewidth = 0.7) +
  scale_y_reverse() +
  scale_x_continuous(position = "top") +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "Depth (m)",
       caption = paste0("Line: LOESS trend fitted to samples at depth ≤ ",
                        trend_max_depth, " m only (",
                        sum(clean$depth_m > trend_max_depth),
                        " samples lie deeper).")) +
  theme(strip.placement = "outside")

save_fig(fig4_2b, "Fig4_2b_depth_profiles.png", width = 7, height = 6)

# ---- Plain-English summary -------------------------------------------------
cat("\n--- Table 4.2: descriptive statistics ---\n")
print(as.data.frame(table4_2), row.names = FALSE)

describe_shape <- function(skew) {
  if (abs(skew) < 0.5) "roughly symmetric"
  else if (abs(skew) < 1) paste("moderately", ifelse(skew > 0, "right", "left"), "skewed")
  else paste("strongly", ifelse(skew > 0, "right", "left"), "skewed")
}

cat("\nSummary (n = ", nrow(clean), "):\n", sep = "")
for (i in seq_len(nrow(table4_2))) {
  r <- table4_2[i, ]
  cat(sprintf("- %s: median %s, mean %s, CV %s%%; the distribution is %s (skewness %s, excess kurtosis %s).\n",
              r$Variable, r$Median, r$Mean, r$CV_pct,
              describe_shape(r$Skewness), r$Skewness, r$Kurtosis_excess))
}
rho_N  <- cor(clean$depth_m, clean$N,      method = "spearman")
rho_qu <- cor(clean$depth_m, clean$qu_kPa, method = "spearman")
cat(sprintf("- With depth: Spearman rho = %.2f for SPT-N and %.2f for qu (see Fig 4.2b).\n",
            rho_N, rho_qu))
cat("Saved Table4_2_descriptive.csv, Fig4_1_histograms.png, Fig4_2_boxplots.png,",
    "Fig4_2b_depth_profiles.png\n")

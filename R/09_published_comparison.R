# 09_published_comparison.R
# Compare measured qu with qu predicted by published SPT-N correlations
# (Table 4.9, Fig 4.10, Fig 4.11). Correlations given for undrained shear
# strength (C, Cu or Su, kPa) are converted with qu = 2 Su.
# Coefficients are applied to field SPT-N as recorded (no N60 correction),
# and PI is not available, so every PI or soil-type variant is evaluated.
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

# ---- Published correlations ------------------------------------------------
# strength: "qu" if the correlation gives qu directly, "Su" if it gives C/Cu/Su.
# group: correlations with several variants (shown as bands in Fig 4.10).
K_this_study <- coef(lm(qu_kPa ~ 0 + N, data = clean))[["N"]]

linear <- function(k) function(N) k * N

correlations <- tribble(
  ~Correlation,                                   ~Source,                        ~group,       ~strength, ~original,              ~f,
  "Terzaghi & Peck (1967)",                       "Terzaghi & Peck 1967",         NA,           "qu",      "qu = 12.5N",           linear(12.5),
  "Sanglerat (1972), clay",                       "Sanglerat 1972",               NA,           "Su",      "C = 12.5N",            linear(12.5),
  "Sanglerat (1972), silty clay",                 "Sanglerat 1972",               NA,           "Su",      "C = 10N",              linear(10),
  "Stroud (1974), PI < 20",                       "Stroud 1974",                  "Stroud",     "Su",      "C = 6.5N",             linear(6.5),
  "Stroud (1974), PI 20-30",                      "Stroud 1974",                  "Stroud",     "Su",      "C = 4.5N",             linear(4.5),
  "Stroud (1974), PI > 30",                       "Stroud 1974",                  "Stroud",     "Su",      "C = 4.2N",             linear(4.2),
  "Sowers (1979), 3.75N",                         "Sowers 1979",                  "Sowers",     "Su",      "C = 3.75N",            linear(3.75),
  "Sowers (1979), 7.5N",                          "Sowers 1979",                  "Sowers",     "Su",      "C = 7.5N",             linear(7.5),
  "Sowers (1979), 12.5N",                         "Sowers 1979",                  "Sowers",     "Su",      "C = 12.5N",            linear(12.5),
  "Nixon (1982)",                                 "Nixon 1982",                   NA,           "Su",      "C = 12N",              linear(12),
  "Hara et al. (1971)",                           "Hara 1971",                    NA,           "Su",      "Cu = 29N^0.72",        function(N) 29 * N^0.72,
  "Sivrikaya & Togrol (2002), fine-grained",      "Sivrikaya & Togrol 2002",      NA,           "Su",      "Su = 4.45N",           linear(4.45),
  "Sivrikaya & Togrol (2002), clays",             "Sivrikaya & Togrol 2002",      NA,           "Su",      "Su = 4.75N",           linear(4.75),
  "Serajuddin & Chowdhury (1998)",                "Serajuddin & Chowdhury 1998",  NA,           "qu",      "qu = 16.8N",           linear(16.8),
  "This study",                                   "This study",                   NA,           "qu",      sprintf("qu = %.2fN", K_this_study), linear(K_this_study)
) %>%
  # qu = 2 Su for strength correlations
  mutate(qu_fun = map2(f, strength, function(f, s) {
    if (s == "Su") function(N) 2 * f(N) else f
  }))

# ---- Predictions and error metrics -----------------------------------------
preds <- correlations %>%
  select(Correlation, Source, group, qu_fun) %>%
  mutate(data = map(qu_fun, ~ clean %>% mutate(predicted = .x(N)))) %>%
  select(-qu_fun) %>%
  unnest(data)

table4_9 <- preds %>%
  group_by(Correlation) %>%
  summarise(n = n(),
            RMSE = sqrt(mean((predicted - qu_kPa)^2)),
            MAE = mean(abs(predicted - qu_kPa)),
            mean_bias = mean(predicted - qu_kPa),
            mean_pct_error = mean(100 * (predicted - qu_kPa) / qu_kPa),
            median_ratio_pred_meas = median(predicted / qu_kPa),
            pct_within_30 = 100 * mean(abs(predicted - qu_kPa) <= 0.3 * qu_kPa),
            .groups = "drop") %>%
  left_join(correlations %>% select(Correlation, original, strength), by = "Correlation") %>%
  mutate(qu_equation = map_chr(Correlation, function(cn) {
    r <- correlations[correlations$Correlation == cn, ]
    mult <- if (r$strength == "Su") 2 else 1
    if (cn == "Hara et al. (1971)") sprintf("qu = %gN^0.72", 29 * mult)
    else sprintf("qu = %.2fN", r$qu_fun[[1]](1))
  })) %>%
  select(Correlation, original, qu_equation, n, RMSE, MAE, mean_bias,
         mean_pct_error, median_ratio_pred_meas, pct_within_30) %>%
  arrange(RMSE) %>%
  mutate(rank = row_number(),
         across(c(RMSE, MAE, mean_bias, mean_pct_error, pct_within_30), ~ round(.x, 2)),
         median_ratio_pred_meas = round(median_ratio_pred_meas, 3))

write_csv(table4_9, "outputs/tables/Table4_9_published_comparison.csv")

# ---- Fig 4.10: all correlations over the data ------------------------------
N_grid <- seq(min(clean$N), max(clean$N), length.out = 300)
curves <- correlations %>%
  mutate(curve = map(qu_fun, ~ tibble(N = N_grid, qu = .x(N_grid)))) %>%
  select(Correlation, group, curve) %>%
  unnest(curve)

# Band labels give the original C coefficients (qu = 2C)
band_labels <- c(Stroud = "Stroud (1974): C = 4.2N to 6.5N",
                 Sowers = "Sowers (1979): C = 3.75N to 12.5N")
bands <- curves %>%
  filter(!is.na(group)) %>%
  group_by(group, N) %>%
  summarise(lo = min(qu), hi = max(qu), .groups = "drop") %>%
  mutate(group = band_labels[group])

line_levels <- c("This study", "Serajuddin & Chowdhury (1998)", "Terzaghi & Peck (1967)",
                 "Sanglerat (1972), clay", "Sanglerat (1972), silty clay", "Nixon (1982)",
                 "Hara et al. (1971)", "Sivrikaya & Togrol (2002), fine-grained",
                 "Sivrikaya & Togrol (2002), clays")
lines <- curves %>%
  filter(is.na(group)) %>%
  mutate(Correlation = factor(Correlation, levels = line_levels))

line_cols <- c("#000000", "#D55E00", "#0072B2", "#009E73", "#009E73", "#CC79A7",
               "#E69F00", "#56B4E9", "#56B4E9")
line_ltys <- c("solid", "longdash", "dashed", "dotdash", "dotted", "twodash",
               "solid", "dotdash", "dotted")
line_wds  <- c(1.3, rep(0.7, 8))

fig4_10 <- ggplot() +
  geom_ribbon(data = bands, aes(x = N, ymin = lo, ymax = hi, fill = group), alpha = 0.35) +
  geom_point(data = clean, aes(N, qu_kPa), shape = 21, fill = "white",
             colour = "grey30", size = 1.5) +
  geom_line(data = lines, aes(N, qu, colour = Correlation, linetype = Correlation,
                              linewidth = Correlation)) +
  scale_fill_manual(values = c("grey40", "#F0E442"), name = NULL) +
  scale_colour_manual(values = line_cols, name = NULL) +
  scale_linetype_manual(values = line_ltys, name = NULL) +
  scale_linewidth_manual(values = line_wds, name = NULL) +
  labs(x = "SPT-N (blows/300 mm)", y = "qu (kPa)",
       caption = "Su-based correlations converted with qu = 2Su; applied to field SPT-N.") +
  theme(legend.position = "right", legend.key.width = unit(1.4, "cm"),
        legend.text = element_text(size = 9))

save_fig(fig4_10, "Fig4_10_published_correlations.png", width = 10, height = 6)

# ---- Fig 4.11: measured vs predicted, one panel per correlation ------------
panel_order <- table4_9$Correlation
lims <- c(0, max(c(preds$qu_kPa, preds$predicted)) * 1.03)
panel_labels <- table4_9 %>%
  mutate(label = sprintf("RMSE %.0f kPa\nbias %+.0f kPa", RMSE, mean_bias),
         Correlation = factor(Correlation, levels = panel_order))

fig4_11 <- preds %>%
  mutate(Correlation = factor(Correlation, levels = panel_order)) %>%
  ggplot(aes(qu_kPa, predicted)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.5) +
  geom_point(shape = 21, fill = "grey60", colour = "black", size = 1, stroke = 0.3) +
  geom_text(data = panel_labels, aes(x = lims[1] + 15, y = lims[2] * 0.97, label = label),
            hjust = 0, vjust = 1, family = "Times New Roman", size = 2.6, lineheight = 0.9) +
  coord_equal(xlim = lims, ylim = lims) +
  facet_wrap(~ Correlation, ncol = 5, labeller = label_wrap_gen(width = 24)) +
  labs(x = "Measured qu (kPa)", y = "Predicted qu (kPa)",
       caption = "Panels ordered by RMSE (best first). Solid line: 1:1.") +
  theme(strip.text = element_text(size = 8))

save_fig(fig4_11, "Fig4_11_published_measured_vs_predicted.png", width = 10, height = 7.5)

# ---- Plain-English summary -------------------------------------------------
cat("\n--- Table 4.9: published correlations vs measured qu (n = ", nrow(clean), ") ---\n", sep = "")
print(as.data.frame(table4_9 %>% select(-original)), row.names = FALSE)

best_pub <- table4_9 %>% filter(Correlation != "This study") %>% slice(1)
sera <- table4_9 %>% filter(Correlation == "Serajuddin & Chowdhury (1998)")
cat("\nSummary:\n")
cat(sprintf("- Closest published correlation: %s (%s), RMSE %.1f kPa, mean bias %+.1f kPa, median predicted/measured %.2f.\n",
            best_pub$Correlation, best_pub$qu_equation, best_pub$RMSE,
            best_pub$mean_bias, best_pub$median_ratio_pred_meas))
cat(sprintf("- Serajuddin & Chowdhury (1998) over-predicts by %.1f kPa on average (mean %% error %+.1f%%, median ratio %.2f).\n",
            sera$mean_bias, sera$mean_pct_error, sera$median_ratio_pred_meas))
cat("- 'This study' is fitted to the same data, so its errors are in-sample (its 10-fold CV RMSE is in Table 4.8).\n")
cat("Saved Table4_9_published_comparison.csv, Fig4_10_published_correlations.png,",
    "Fig4_11_published_measured_vs_predicted.png\n")

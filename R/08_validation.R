# 08_validation.R
# Out-of-sample validation of the qu-N models:
#   - repeated 10-fold cross-validation, 10 repeats (Table 4.8)
#   - single 70/30 train/test split, for comparison with earlier work (Table 4.8b)
#   - bootstrap (2000 resamples) of K in qu = KN (Table 4.8c)
#   - predicted vs measured with 1:1 and +/-30% lines (Fig 4.8)
#   - qu = KN with 95% confidence band and log-log 95% prediction band (Fig 4.9)
# All metrics are on the qu scale. The log-log model is back-transformed with
# Duan's smearing correction estimated on the training data.
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

# ---- Model definitions -----------------------------------------------------
# Each function fits on `train` and returns qu-scale predictions for `test`
fit_predict <- list(
  "qu = KN" = function(train, test) {
    predict(lm(qu_kPa ~ 0 + N, data = train), test)
  },
  "Linear" = function(train, test) {
    predict(lm(qu_kPa ~ N, data = train), test)
  },
  "Polynomial" = function(train, test) {
    predict(lm(qu_kPa ~ N + I(N^2), data = train), test)
  },
  "Power (nls)" = function(train, test) {
    ll <- lm(log(qu_kPa) ~ log(N), data = train)
    m <- nls(qu_kPa ~ a * N^b, data = train,
             start = list(a = exp(coef(ll)[[1]]), b = coef(ll)[[2]]))
    predict(m, test)
  },
  "Log-log (smearing)" = function(train, test) {
    m <- lm(log(qu_kPa) ~ log(N), data = train)
    smear <- mean(exp(residuals(m)))
    exp(predict(m, test)) * smear
  }
)
model_levels <- names(fit_predict)

metrics <- function(obs, pred) {
  c(RMSE = sqrt(mean((obs - pred)^2)),
    MAE  = mean(abs(obs - pred)),
    # Predictive R2 on the held-out data: 1 - SSE / SST (held-out mean)
    R2   = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2))
}

# ---- Repeated 10-fold cross-validation -------------------------------------
set.seed(123)
folds <- rsample::vfold_cv(clean, v = 10, repeats = 10)

cv_preds <- map_dfr(seq_len(nrow(folds)), function(i) {
  split <- folds$splits[[i]]
  train <- rsample::analysis(split)
  test  <- rsample::assessment(split)
  imap_dfr(fit_predict, function(f, name) {
    tibble(repeat_id = folds$id[i], fold_id = folds$id2[i], Model = name,
           sample_id = test$sample_id, observed = test$qu_kPa,
           predicted = f(train, test))
  })
})

cv_folds <- cv_preds %>%
  group_by(Model, repeat_id, fold_id) %>%
  summarise(as_tibble_row(metrics(observed, predicted)), .groups = "drop")

# Paired comparison with qu = KN: every model is scored on the same folds
paired <- cv_folds %>%
  select(Model, repeat_id, fold_id, RMSE) %>%
  left_join(cv_folds %>% filter(Model == "qu = KN") %>%
              select(repeat_id, fold_id, RMSE_KN = RMSE),
            by = c("repeat_id", "fold_id")) %>%
  group_by(Model) %>%
  summarise(RMSE_diff_vs_KN = mean(RMSE - RMSE_KN),
            pct_folds_better_than_KN = 100 * mean(RMSE < RMSE_KN),
            .groups = "drop")

table4_8 <- cv_folds %>%
  left_join(paired, by = "Model") %>%
  group_by(Model, RMSE_diff_vs_KN, pct_folds_better_than_KN) %>%
  summarise(n_folds = n(),
            RMSE_mean = mean(RMSE), RMSE_SD = sd(RMSE),
            MAE_mean  = mean(MAE),  MAE_SD  = sd(MAE),
            R2_mean   = mean(R2),   R2_SD   = sd(R2),
            .groups = "drop") %>%
  mutate(Model = factor(Model, levels = model_levels)) %>%
  arrange(Model) %>%
  relocate(RMSE_diff_vs_KN, pct_folds_better_than_KN, .after = last_col()) %>%
  mutate(across(c(starts_with("RMSE"), starts_with("MAE")), ~ round(.x, 2)),
         across(starts_with("R2"), ~ round(.x, 4)),
         pct_folds_better_than_KN = round(pct_folds_better_than_KN, 0),
         rank_RMSE = rank(RMSE_mean, ties.method = "min"))

write_csv(table4_8, "outputs/tables/Table4_8_cv_results.csv")

# ---- 70/30 train/test split ------------------------------------------------
set.seed(123)
split_70 <- rsample::initial_split(clean, prop = 0.7)
train_70 <- rsample::training(split_70)
test_70  <- rsample::testing(split_70)

table4_8b <- imap_dfr(fit_predict, function(f, name) {
  bind_rows(
    as_tibble_row(metrics(train_70$qu_kPa, f(train_70, train_70))) %>% mutate(Set = "Train (70%)"),
    as_tibble_row(metrics(test_70$qu_kPa,  f(train_70, test_70)))  %>% mutate(Set = "Test (30%)")
  ) %>% mutate(Model = name, n = c(nrow(train_70), nrow(test_70)))
}) %>%
  select(Model, Set, n, RMSE, MAE, R2) %>%
  mutate(across(c(RMSE, MAE), ~ round(.x, 2)), R2 = round(R2, 4))

write_csv(table4_8b, "outputs/tables/Table4_8b_train_test_split.csv")

k_train_70 <- coef(lm(qu_kPa ~ 0 + N, data = train_70))[["N"]]

# ---- Bootstrap K in qu = KN ------------------------------------------------
set.seed(123)
boot_K <- boot::boot(clean, function(d, i) coef(lm(qu_kPa ~ 0 + N, data = d[i, ]))[["N"]],
                     R = 2000)
boot_ci <- boot::boot.ci(boot_K, type = "perc", conf = 0.95)

m_kn <- lm(qu_kPa ~ 0 + N, data = clean)
hc3_ci <- lmtest::coefci(m_kn, vcov. = sandwich::vcovHC(m_kn, type = "HC3"))
ols_ci <- confint(m_kn)

table4_8c <- tibble(
  Method = c("Bootstrap percentile (2000 resamples)", "HC3 robust", "Ordinary OLS"),
  K = coef(m_kn)[["N"]],
  SE = c(sd(boot_K$t[, 1]),
         sqrt(sandwich::vcovHC(m_kn, type = "HC3")[1, 1]),
         summary(m_kn)$coefficients[1, 2]),
  CI95_low  = c(boot_ci$percent[4], hc3_ci[1, 1], ols_ci[1, 1]),
  CI95_high = c(boot_ci$percent[5], hc3_ci[1, 2], ols_ci[1, 2])
) %>%
  mutate(across(K:CI95_high, ~ round(.x, 3)))

write_csv(table4_8c, "outputs/tables/Table4_8c_bootstrap_K.csv")

# ---- Fig 4.8: predicted vs measured ----------------------------------------
# Out-of-fold predictions averaged over the 10 CV repeats (each sample is
# predicted once per repeat by a model that did not see it)
oof <- cv_preds %>%
  group_by(Model, sample_id) %>%
  summarise(observed = first(observed), predicted = mean(predicted), .groups = "drop") %>%
  mutate(within_30 = abs(predicted - observed) <= 0.3 * observed,
         Model = factor(Model, levels = model_levels))

within_pct <- oof %>%
  group_by(Model) %>%
  summarise(pct_within_30 = 100 * mean(within_30), .groups = "drop")

lims <- c(0, max(c(oof$observed, oof$predicted)) * 1.05)
ref_lines <- tibble(x = lims) %>%
  mutate(`1:1` = x, `+30%` = 1.3 * x, `-30%` = 0.7 * x) %>%
  pivot_longer(-x, names_to = "line", values_to = "y")

fig4_8 <- ggplot(oof, aes(observed, predicted)) +
  geom_line(data = ref_lines, aes(x, y, linetype = line, group = line),
            colour = "black", linewidth = 0.5) +
  geom_point(shape = 21, fill = "grey60", colour = "black", size = 1.5) +
  geom_text(data = within_pct, aes(x = lims[1] + 10, y = lims[2] * 0.95,
                                   label = sprintf("%.0f%% within ±30%%", pct_within_30)),
            hjust = 0, family = "Times New Roman", size = 3.3) +
  scale_linetype_manual(values = c(`1:1` = "solid", `+30%` = "dashed", `-30%` = "dashed"),
                        breaks = c("1:1", "+30%"), labels = c("1:1 line", "±30%")) +
  coord_equal(xlim = lims, ylim = lims) +
  facet_wrap(~ Model, nrow = 2) +
  labs(x = "Measured qu (kPa)", y = "Predicted qu (kPa)", linetype = NULL,
       caption = "Out-of-fold predictions from 10-fold CV, averaged over 10 repeats.") +
  theme(legend.position = "inside", legend.position.inside = c(0.84, 0.25))

save_fig(fig4_8, "Fig4_8_predicted_vs_measured.png", width = 9, height = 6.5)

# ---- Fig 4.9: qu = KN with confidence and prediction bands -----------------
# Confidence band: K +/- t * SE_HC3(K) * N (robust SE because of heteroscedasticity)
# Prediction band: 95% prediction interval of the log-log model, back-transformed
# with exp() (quantiles transform directly, so no smearing is needed), so it widens with N
m_ll <- lm(log(qu_kPa) ~ log(N), data = clean)
se_k_hc3 <- sqrt(sandwich::vcovHC(m_kn, type = "HC3")[1, 1])
t_crit <- qt(0.975, df = df.residual(m_kn))

N_grid <- tibble(N = seq(min(clean$N), max(clean$N), length.out = 300))
pi_ll <- predict(m_ll, N_grid, interval = "prediction", level = 0.95)
bands <- N_grid %>%
  mutate(fit = coef(m_kn)[["N"]] * N,
         ci_low = (coef(m_kn)[["N"]] - t_crit * se_k_hc3) * N,
         ci_high = (coef(m_kn)[["N"]] + t_crit * se_k_hc3) * N,
         pi_low = exp(pi_ll[, "lwr"]),
         pi_high = exp(pi_ll[, "upr"]))

pi_obs <- exp(predict(m_ll, clean, interval = "prediction", level = 0.95))
coverage <- 100 * mean(clean$qu_kPa >= pi_obs[, "lwr"] & clean$qu_kPa <= pi_obs[, "upr"])

fig4_9 <- ggplot(bands, aes(N)) +
  geom_ribbon(aes(ymin = pi_low, ymax = pi_high, fill = "95% prediction band (log-log model)"),
              alpha = 0.6) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = "95% confidence band (HC3)"),
              alpha = 0.9) +
  geom_point(data = clean, aes(N, qu_kPa), shape = 21, fill = "white",
             colour = "black", size = 1.6) +
  geom_line(aes(y = fit, linetype = sprintf("qu = %.2f N", coef(m_kn)[["N"]])),
            linewidth = 0.8) +
  scale_fill_manual(values = c("95% prediction band (log-log model)" = "grey88",
                               "95% confidence band (HC3)" = "grey55")) +
  scale_linetype_manual(values = "solid") +
  guides(fill = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(x = "SPT-N (blows/300 mm)", y = "qu (kPa)", fill = NULL, linetype = NULL,
       caption = sprintf("%.0f%% of the %d samples lie inside the prediction band.",
                         coverage, nrow(clean))) +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.spacing.y = unit(0, "pt"),
        legend.background = element_rect(fill = "white", colour = "grey60"))

save_fig(fig4_9, "Fig4_9_KN_bands.png", width = 7, height = 5.5)

# ---- Plain-English summary -------------------------------------------------
cat("\n--- Table 4.8: repeated 10-fold CV (10 x 10 = 100 folds, qu scale) ---\n")
print(as.data.frame(table4_8), row.names = FALSE)
cat("\n--- Table 4.8b: 70/30 train/test split (n train =", nrow(train_70),
    ", n test =", nrow(test_70), ") ---\n")
print(as.data.frame(table4_8b), row.names = FALSE)
cat("\n--- Table 4.8c: K in qu = KN ---\n")
print(as.data.frame(table4_8c), row.names = FALSE)
cat("\nShare of out-of-fold predictions within +/-30% of measured qu:\n")
print(as.data.frame(within_pct %>% mutate(pct_within_30 = round(pct_within_30, 1))),
      row.names = FALSE)

best <- as.character(table4_8$Model[which.min(table4_8$RMSE_mean)])
cat("\nSummary:\n")
cat(sprintf("- Lowest mean CV RMSE: %s (%.2f kPa). Spread of mean CV RMSE across models: %.2f to %.2f kPa.\n",
            best, min(table4_8$RMSE_mean), min(table4_8$RMSE_mean), max(table4_8$RMSE_mean)))
cat(sprintf("- Fold-to-fold SD of RMSE is about %.0f kPa; on the same folds, %s beats qu = KN in %.0f%% of folds (mean RMSE difference %.2f kPa).\n",
            median(table4_8$RMSE_SD), best,
            table4_8$pct_folds_better_than_KN[table4_8$Model == best],
            table4_8$RMSE_diff_vs_KN[table4_8$Model == best]))
cat(sprintf("- 70/30 split: K fitted on the training set = %.3f; test RMSE for qu = KN = %.2f kPa.\n",
            k_train_70, table4_8b$RMSE[table4_8b$Model == "qu = KN" & table4_8b$Set == "Test (30%)"]))
cat(sprintf("- Bootstrap 95%% percentile CI for K: %.2f to %.2f (HC3: %.2f to %.2f).\n",
            table4_8c$CI95_low[1], table4_8c$CI95_high[1],
            table4_8c$CI95_low[2], table4_8c$CI95_high[2]))
cat(sprintf("- %.0f%% of samples fall inside the log-log 95%% prediction band.\n", coverage))
cat("Saved Table4_8_cv_results.csv, Table4_8b_train_test_split.csv, Table4_8c_bootstrap_K.csv,",
    "Fig4_8_predicted_vs_measured.png, Fig4_9_KN_bands.png\n")

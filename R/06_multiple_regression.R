# 06_multiple_regression.R
# Does depth add to SPT-N in predicting qu?
#   qu = KN                                 (N-only, through origin)
#   qu ~ N                                  (N-only, with intercept; nested baseline)
#   qu ~ N + depth_m
#   qu ~ N * depth_m
#   log(qu) ~ log(N)                        (N-only, log-log)
#   log(qu) ~ log(N) + log(depth_m)
# Coefficients with ordinary and robust HC3 standard errors (Table 4.6),
# fit comparison on the qu scale (Table 4.6b), VIF (Table 4.6c) and
# nested F tests against the N-only models (Table 4.7).
# All 134 samples are used, including the 4 duplicated rows.

source("R/00_packages.R")
set.seed(123)

clean <- readRDS("data/clean.rds")
n <- nrow(clean)

# ---- Fit models ------------------------------------------------------------
models <- list(
  "qu = KN"                   = lm(qu_kPa ~ 0 + N, data = clean),
  "qu ~ N"                    = lm(qu_kPa ~ N, data = clean),
  "qu ~ N + depth"            = lm(qu_kPa ~ N + depth_m, data = clean),
  "qu ~ N * depth"            = lm(qu_kPa ~ N * depth_m, data = clean),
  "ln qu ~ ln N"              = lm(log(qu_kPa) ~ log(N), data = clean),
  "ln qu ~ ln N + ln depth"   = lm(log(qu_kPa) ~ log(N) + log(depth_m), data = clean)
)
log_models <- c("ln qu ~ ln N", "ln qu ~ ln N + ln depth")

hc3 <- function(m) sandwich::vcovHC(m, type = "HC3")

# ---- Table 4.6: coefficients, ordinary and HC3 -----------------------------
coef_table <- function(m, name) {
  ols <- lmtest::coeftest(m)
  rob <- lmtest::coeftest(m, vcov. = hc3(m))
  ci_ols <- confint(m)
  ci_rob <- lmtest::coefci(m, vcov. = hc3(m))
  tibble(Model = name,
         Term = rownames(ols),
         Estimate = ols[, 1],
         SE_OLS = ols[, 2], t_OLS = ols[, 3], p_OLS = ols[, 4],
         CI95_low_OLS = ci_ols[, 1], CI95_high_OLS = ci_ols[, 2],
         SE_HC3 = rob[, 2], t_HC3 = rob[, 3], p_HC3 = rob[, 4],
         CI95_low_HC3 = ci_rob[, 1], CI95_high_HC3 = ci_rob[, 2])
}

term_labels <- c("(Intercept)" = "Intercept", "depth_m" = "Depth",
                 "N:depth_m" = "N x Depth", "log(N)" = "ln N",
                 "log(depth_m)" = "ln Depth", "N_c" = "N (centred)",
                 "depth_c" = "Depth (centred)", "N_c:depth_c" = "N x Depth (centred)")

table4_6 <- imap_dfr(models, coef_table) %>%
  mutate(Term = coalesce(term_labels[Term], Term),
         across(c(Estimate, starts_with("SE"), starts_with("t_"), starts_with("CI95")),
                ~ round(.x, 4)),
         across(starts_with("p_"), ~ signif(.x, 3)))

write_csv(table4_6, "outputs/tables/Table4_6_mlr_coefficients.csv")

# ---- Table 4.6b: fit comparison on the qu scale ----------------------------
# Same definitions as script 05: R2 = 1 - SSE/SST (centred); AIC/BIC from the
# Gaussian log-likelihood of qu-scale residuals. Log models are back-transformed
# with exp(fitted) (median prediction). R2_own_scale is the model's own R2
# (uncentred for qu = KN, ln qu scale for the log models).
sst <- sum((clean$qu_kPa - mean(clean$qu_kPa))^2)
fit_stats <- function(m, name) {
  pred <- if (name %in% log_models) exp(fitted(m)) else fitted(m)
  res <- clean$qu_kPa - pred
  sse <- sum(res^2)
  p <- length(coef(m))
  r2 <- 1 - sse / sst
  loglik <- -n / 2 * (log(2 * pi) + log(sse / n) + 1)
  tibble(Model = name, n_coef = p,
         R2_own_scale = summary(m)$r.squared,
         Adj_R2_own_scale = summary(m)$adj.r.squared,
         R2_qu = r2,
         Adj_R2_qu = 1 - (1 - r2) * (n - 1) / (n - p),
         RMSE_qu = sqrt(mean(res^2)),
         MAE_qu = mean(abs(res)),
         AIC_qu = -2 * loglik + 2 * (p + 1),
         BIC_qu = -2 * loglik + log(n) * (p + 1),
         # qu = KN has no intercept, so its variance is regressed on N explicitly
         BP_p = if (name == "qu = KN") lmtest::bptest(m, varformula = ~ N, data = clean)$p.value
                else lmtest::bptest(m)$p.value)
}

table4_6b <- imap_dfr(models, fit_stats) %>%
  mutate(across(c(starts_with("R2"), starts_with("Adj")), ~ round(.x, 4)),
         across(c(RMSE_qu, MAE_qu, AIC_qu, BIC_qu), ~ round(.x, 2)),
         BP_p = signif(BP_p, 3))

write_csv(table4_6b, "outputs/tables/Table4_6b_mlr_fit.csv")

# ---- Table 4.6c: VIF -------------------------------------------------------
# The interaction model is refitted with centred N and depth: the product term
# inflates raw VIFs by construction, centring removes that without changing the fit.
clean_c <- clean %>% mutate(N_c = N - mean(N), depth_c = depth_m - mean(depth_m))
vif_list <- list(
  "qu ~ N + depth"                  = car::vif(models[["qu ~ N + depth"]]),
  "qu ~ N * depth (centred)"        = suppressMessages(
    car::vif(lm(qu_kPa ~ N_c * depth_c, data = clean_c), type = "terms")),
  "ln qu ~ ln N + ln depth"         = car::vif(models[["ln qu ~ ln N + ln depth"]])
)
table4_6c <- imap_dfr(vif_list, function(v, name) {
  if (is.matrix(v)) {
    tibble(Model = name, Term = rownames(v), VIF = v[, "GVIF"])
  } else {
    tibble(Model = name, Term = names(v), VIF = unname(v))
  }
}) %>%
  mutate(Term = coalesce(term_labels[Term], Term),
         VIF = round(VIF, 3))

write_csv(table4_6c, "outputs/tables/Table4_6c_vif.csv")

# ---- Table 4.7: nested F tests ---------------------------------------------
# Classical partial F test (anova) and robust Wald F test using HC3
nested <- tribble(
  ~Comparison,                                     ~small,           ~big,
  "Add depth to qu ~ N",                           "qu ~ N",         "qu ~ N + depth",
  "Add depth and N x depth to qu ~ N",             "qu ~ N",         "qu ~ N * depth",
  "Add N x depth to qu ~ N + depth",               "qu ~ N + depth", "qu ~ N * depth",
  "Add intercept and depth to qu = KN",            "qu = KN",        "qu ~ N + depth",
  "Add ln depth to ln qu ~ ln N",                  "ln qu ~ ln N",   "ln qu ~ ln N + ln depth"
)

table4_7 <- pmap_dfr(nested, function(Comparison, small, big) {
  a <- anova(models[[small]], models[[big]])
  w <- lmtest::waldtest(models[[small]], models[[big]], vcov = hc3, test = "F")
  tibble(Comparison, Reduced = small, Full = big,
         df_added = a$Df[2],
         RSS_reduced = a$RSS[1], RSS_full = a$RSS[2],
         F_OLS = a$F[2], p_OLS = a$`Pr(>F)`[2],
         F_HC3 = w$F[2], p_HC3 = w$`Pr(>F)`[2])
}) %>%
  mutate(across(c(RSS_reduced, RSS_full), ~ round(.x, 2)),
         across(c(F_OLS, F_HC3), ~ round(.x, 3)),
         across(c(p_OLS, p_HC3), ~ signif(.x, 3)),
         significant_at_5pct_HC3 = ifelse(p_HC3 < 0.05, "Yes", "No"))

write_csv(table4_7, "outputs/tables/Table4_7_nested_F_test.csv")

# ---- Plain-English summary -------------------------------------------------
fmt_p <- function(p) format.pval(p, digits = 3, eps = 1e-4)

cat("\n--- Table 4.6: coefficients (ordinary and HC3 robust) ---\n")
print(as.data.frame(table4_6 %>%
                      select(Model, Term, Estimate, SE_OLS, p_OLS, SE_HC3, p_HC3,
                             CI95_low_HC3, CI95_high_HC3)), row.names = FALSE)

cat("\n--- Table 4.6b: fit (qu scale unless stated) ---\n")
print(as.data.frame(table4_6b %>%
                      select(Model, R2_own_scale, R2_qu, Adj_R2_qu, RMSE_qu, AIC_qu, BIC_qu, BP_p)),
      row.names = FALSE)

cat("\n--- Table 4.6c: VIF ---\n")
print(as.data.frame(table4_6c), row.names = FALSE)

cat("\n--- Table 4.7: nested F tests ---\n")
print(as.data.frame(table4_7 %>%
                      select(Comparison, df_added, F_OLS, p_OLS, F_HC3, p_HC3,
                             significant_at_5pct_HC3)), row.names = FALSE)

cat("\nSummary:\n")
for (i in seq_len(nrow(table4_7))) {
  r <- table4_7[i, ]
  cat(sprintf("- %s: ordinary F p = %s, robust HC3 p = %s -> %s.\n",
              r$Comparison, fmt_p(r$p_OLS), fmt_p(r$p_HC3),
              ifelse(r$p_HC3 < 0.05, "significant improvement",
                     "no significant improvement")))
}
cat(sprintf("- Largest VIF = %.2f (values below 5 mean collinearity is not a concern).\n",
            max(table4_6c$VIF)))
cat("Saved Table4_6_mlr_coefficients.csv, Table4_6b_mlr_fit.csv, Table4_6c_vif.csv,",
    "Table4_7_nested_F_test.csv\n")

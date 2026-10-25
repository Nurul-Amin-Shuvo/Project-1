# 05_simple_models.R
# Simple regression models of qu on SPT-N (Table 4.5, Fig 4.6):
#   (a) linear          qu = a + bN
#   (b) through origin  qu = KN
#   (c) polynomial      qu = a + bN + cN^2
#   (d) logarithmic     qu = a + b ln N
#   (e) power           qu = aN^b, fitted by log-log lm and by nls
# All fit statistics are computed on the original qu scale so the models
# compare fairly. Breusch-Pagan test on the linear model (Table 4.5b).
# All 134 samples are used, including the 4 duplicated rows.

source("R/00_packages.R")
set.seed(123)

clean <- readRDS("data/clean.rds")
n <- nrow(clean)

theme_set(theme_bw(base_size = 11, base_family = "Times New Roman"))
save_fig <- function(plot, file, width, height) {
  ggsave(file.path("outputs/figures", file), plot, width = width,
         height = height, units = "in", dpi = 300, device = ragg::agg_png)
}

# ---- Fit models ------------------------------------------------------------
m_lin    <- lm(qu_kPa ~ N, data = clean)
m_origin <- lm(qu_kPa ~ 0 + N, data = clean)
m_poly   <- lm(qu_kPa ~ N + I(N^2), data = clean)
m_log    <- lm(qu_kPa ~ log(N), data = clean)
m_loglog <- lm(log(qu_kPa) ~ log(N), data = clean)
m_nls    <- nls(qu_kPa ~ a * N^b, data = clean,
                start = list(a = exp(coef(m_loglog)[[1]]),
                             b = coef(m_loglog)[[2]]))

# Predictions on the original qu scale, as functions of N
# (log-log back-transform is exp(fitted), i.e. the median, without bias correction)
predict_qu <- list(
  "(a) Linear"            = function(N) predict(m_lin,    data.frame(N = N)),
  "(b) Through origin"    = function(N) predict(m_origin, data.frame(N = N)),
  "(c) Polynomial"        = function(N) predict(m_poly,   data.frame(N = N)),
  "(d) Logarithmic"       = function(N) predict(m_log,    data.frame(N = N)),
  "(e1) Power (log-log)"  = function(N) exp(predict(m_loglog, data.frame(N = N))),
  "(e2) Power (nls)"      = function(N) predict(m_nls,    data.frame(N = N))
)

# Number of estimated coefficients (excluding the residual SD)
n_coef <- c(2, 1, 3, 2, 2, 2)

f2 <- function(x) formatC(x, format = "f", digits = 3)
signed <- function(x) ifelse(x < 0, paste("-", f2(abs(x))), paste("+", f2(x)))
equations <- c(
  sprintf("qu = %s %s N", f2(coef(m_lin)[1]), signed(coef(m_lin)[2])),
  sprintf("qu = %s N", f2(coef(m_origin)[1])),
  sprintf("qu = %s %s N %s N^2", f2(coef(m_poly)[1]), signed(coef(m_poly)[2]),
          signed(coef(m_poly)[3])),
  sprintf("qu = %s %s ln(N)", f2(coef(m_log)[1]), signed(coef(m_log)[2])),
  sprintf("qu = %s N^%s", f2(exp(coef(m_loglog)[1])), f2(coef(m_loglog)[2])),
  sprintf("qu = %s N^%s", f2(coef(m_nls)[["a"]]), f2(coef(m_nls)[["b"]]))
)

# ---- F tests ---------------------------------------------------------------
# lm models: the model's own overall F test. For (b) this tests K = 0
# (no intercept); for (e1) it is on the log scale.
# nls has no built-in F test: approximate F against the intercept-only model.
lm_F <- function(m) {
  f <- summary(m)$fstatistic
  c(F = unname(f[1]), p = unname(pf(f[1], f[2], f[3], lower.tail = FALSE)))
}
sse_nls <- sum(residuals(m_nls)^2)
sst <- sum((clean$qu_kPa - mean(clean$qu_kPa))^2)
F_nls <- ((sst - sse_nls) / (2 - 1)) / (sse_nls / (n - 2))
F_tests <- rbind(lm_F(m_lin), lm_F(m_origin), lm_F(m_poly), lm_F(m_log),
                 lm_F(m_loglog),
                 c(F = F_nls, p = pf(F_nls, 1, n - 2, lower.tail = FALSE)))
F_basis <- c("qu scale", "qu scale, H0: K = 0 (uncentred)", "qu scale",
             "qu scale", "ln(qu) scale", "qu scale, approximate")

# ---- Fit statistics on the original qu scale -------------------------------
# R2 = 1 - SSE/SST (centred, also for the no-intercept model)
# AIC/BIC from the Gaussian log-likelihood of the qu-scale residuals,
# with k = coefficients + 1 (residual SD); this equals AIC() for the qu-scale lm/nls fits
fit_stats <- function(pred, p) {
  res <- clean$qu_kPa - pred
  sse <- sum(res^2)
  r2 <- 1 - sse / sst
  loglik <- -n / 2 * (log(2 * pi) + log(sse / n) + 1)
  k <- p + 1
  c(R2 = r2,
    Adj_R2 = 1 - (1 - r2) * (n - 1) / (n - p),
    RMSE = sqrt(mean(res^2)),
    MAE = mean(abs(res)),
    AIC = -2 * loglik + 2 * k,
    BIC = -2 * loglik + log(n) * k)
}

stats <- map2_dfr(predict_qu, n_coef, ~ as_tibble_row(fit_stats(.x(clean$N), .y)))

table4_5 <- tibble(Model = names(predict_qu), Equation = equations) %>%
  bind_cols(stats) %>%
  mutate(F_statistic = F_tests[, "F"],
         F_p_value = F_tests[, "p"],
         F_basis = F_basis) %>%
  mutate(across(c(R2, Adj_R2), ~ round(.x, 4)),
         across(c(RMSE, MAE, AIC, BIC, F_statistic), ~ round(.x, 2)),
         F_p_value = signif(F_p_value, 3))

write_csv(table4_5, "outputs/tables/Table4_5_model_comparison.csv")

# ---- Breusch-Pagan test ----------------------------------------------------
# Linear model as requested; the log-log model is added for comparison
bp_lin <- lmtest::bptest(m_lin)
bp_loglog <- lmtest::bptest(m_loglog)
table4_5b <- tibble(
  Model = c("(a) Linear", "(e1) Power (log-log)"),
  BP_statistic = round(c(bp_lin$statistic, bp_loglog$statistic), 3),
  df = c(bp_lin$parameter, bp_loglog$parameter),
  p_value = signif(c(bp_lin$p.value, bp_loglog$p.value), 3),
  constant_variance_at_5pct = ifelse(c(bp_lin$p.value, bp_loglog$p.value) > 0.05,
                                     "Yes", "No"))

write_csv(table4_5b, "outputs/tables/Table4_5b_breusch_pagan.csv")

# ---- Fig 4.6: data with all fitted curves ----------------------------------
N_grid <- seq(min(clean$N), max(clean$N), length.out = 300)
curves <- imap_dfr(predict_qu, ~ tibble(Model = .y, N = N_grid, qu = .x(N_grid))) %>%
  mutate(Model = factor(Model, levels = names(predict_qu)))

model_cols <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")
model_ltys <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")

fig4_6 <- ggplot() +
  geom_point(data = clean, aes(x = N, y = qu_kPa),
             shape = 21, fill = "grey75", colour = "grey30", size = 1.8) +
  geom_line(data = curves, aes(x = N, y = qu, colour = Model, linetype = Model),
            linewidth = 0.8) +
  scale_colour_manual(values = model_cols) +
  scale_linetype_manual(values = model_ltys) +
  labs(x = "SPT-N (blows/300 mm)", y = "qu (kPa)", colour = NULL, linetype = NULL) +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(colour = "grey60"),
        legend.key.width = unit(1.6, "cm"))

save_fig(fig4_6, "Fig4_6_simple_models.png", width = 7, height = 5.5)

# ---- Plain-English summary -------------------------------------------------
cat("\n--- Table 4.5: model comparison (statistics on the qu scale, n = ", n, ") ---\n", sep = "")
print(as.data.frame(table4_5 %>% select(Model, Equation, R2, Adj_R2, RMSE, MAE, AIC, BIC)),
      row.names = FALSE)
cat("\nF tests:\n")
print(as.data.frame(table4_5 %>% select(Model, F_statistic, F_p_value, F_basis)),
      row.names = FALSE)

cat("\n--- Table 4.5b: Breusch-Pagan test ---\n")
print(as.data.frame(table4_5b), row.names = FALSE)

best_aic <- table4_5$Model[which.min(table4_5$AIC)]
best_rmse <- table4_5$Model[which.min(table4_5$RMSE)]
cat("\nSummary:\n")
cat("- Lowest AIC:", best_aic, "| lowest RMSE:", best_rmse, "\n")
cat(sprintf("- Breusch-Pagan on the linear model: p = %s -> %s.\n",
            format.pval(bp_lin$p.value, digits = 3, eps = 1e-4),
            ifelse(bp_lin$p.value > 0.05, "no evidence of changing variance",
                   "residual spread changes with N (heteroscedastic)")))
cat(sprintf("- Breusch-Pagan on the log-log model: p = %s -> %s.\n",
            format.pval(bp_loglog$p.value, digits = 3, eps = 1e-4),
            ifelse(bp_loglog$p.value > 0.05, "no evidence of changing variance",
                   "residual spread still changes with N")))
cat("Saved Table4_5_model_comparison.csv, Table4_5b_breusch_pagan.csv, Fig4_6_simple_models.png\n")

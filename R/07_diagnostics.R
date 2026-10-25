# 07_diagnostics.R
# Residual diagnostics for the main model qu = KN and the log-log model
# ln qu = ln a + b ln N (Fig 4.7a, Fig 4.7b), Breusch-Pagan and Durbin-Watson
# tests and Cook's distance with the 4/n cut-off (Table 4.7c, Table 4.7d),
# and sensitivity refits of qu = KN (Table 4.7b):
#   (1) all 134 samples
#   (2) without the 4 duplicated rows
#   (3) without the 5 stiff samples at 9.20 m (rows 123, 128, 132, 133, 134)
#   (4) without points above the Cook's distance cut-off

source("R/00_packages.R")
set.seed(123)

clean <- readRDS("data/clean.rds") %>%
  mutate(sample_id = row_number(), .before = 1)
n <- nrow(clean)
cook_cut <- 4 / n

theme_set(theme_bw(base_size = 10, base_family = "Times New Roman"))
save_fig <- function(plot, file, width, height) {
  ggsave(file.path("outputs/figures", file), plot, width = width,
         height = height, units = "in", dpi = 300, device = ragg::agg_png)
}

hc3 <- function(m) sandwich::vcovHC(m, type = "HC3")

# ---- Main and log-log models -----------------------------------------------
m_kn     <- lm(qu_kPa ~ 0 + N, data = clean)
m_loglog <- lm(log(qu_kPa) ~ log(N), data = clean)

diag_data <- function(m) {
  tibble(sample_id = clean$sample_id,
         fitted    = fitted(m),
         resid     = residuals(m),
         std_resid = rstandard(m),
         leverage  = hatvalues(m),
         cooks     = cooks.distance(m)) %>%
    mutate(influential = cooks > cook_cut)
}

# ---- Four-panel residual diagnostics ---------------------------------------
four_panel <- function(m, title, fitted_lab, resid_lab) {
  d <- diag_data(m)
  pts <- function(p) {
    p + geom_point(aes(fill = influential), shape = 21, colour = "black", size = 1.5) +
      scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "black"), guide = "none")
  }
  lbl <- ggrepel::geom_text_repel(data = ~ filter(.x, influential), aes(label = sample_id),
                                  family = "Times New Roman", size = 2.6,
                                  min.segment.length = 0.2, box.padding = 0.3, seed = 123)
  smooth <- geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                        colour = "black", linewidth = 0.6)

  p1 <- pts(ggplot(d, aes(fitted, resid))) +
    geom_hline(yintercept = 0, linetype = "dashed") + smooth + lbl +
    labs(title = "Residuals vs fitted", x = fitted_lab, y = resid_lab)

  p2 <- ggplot(d, aes(sample = std_resid)) +
    stat_qq_line(colour = "black") +
    stat_qq(shape = 21, fill = "grey70", colour = "black", size = 1.5) +
    labs(title = "Normal Q-Q", x = "Theoretical quantiles",
         y = "Standardised residuals")

  p3 <- pts(ggplot(d, aes(fitted, sqrt(abs(std_resid))))) + smooth +
    labs(title = "Scale-location", x = fitted_lab,
         y = expression(sqrt("|Standardised residuals|")))

  p4 <- pts(ggplot(d, aes(leverage, std_resid))) +
    geom_hline(yintercept = 0, linetype = "dashed") + lbl +
    labs(title = "Residuals vs leverage", x = "Leverage",
         y = "Standardised residuals")

  (p1 + p2) / (p3 + p4) +
    patchwork::plot_annotation(
      title = title,
      caption = sprintf("Black points: Cook's distance > 4/n = %.3f, labelled by sample number.",
                        cook_cut),
      theme = theme(plot.title = element_text(family = "Times New Roman", size = 12)))
}

fig4_7a <- four_panel(m_kn, "Residual diagnostics: qu = KN",
                      "Fitted qu (kPa)", "Residual (kPa)")
fig4_7b <- four_panel(m_loglog, "Residual diagnostics: ln qu = ln a + b ln N",
                      "Fitted ln qu", "Residual (ln qu)")

save_fig(fig4_7a, "Fig4_7a_diagnostics_KN.png", width = 8, height = 7)
save_fig(fig4_7b, "Fig4_7b_diagnostics_loglog.png", width = 8, height = 7)

# ---- Table 4.7c: Breusch-Pagan and Durbin-Watson ---------------------------
# qu = KN has no intercept, so its BP variance regression uses N explicitly.
# Durbin-Watson checks autocorrelation in the order the samples appear in the file.
bp_kn <- lmtest::bptest(m_kn, varformula = ~ N, data = clean)
bp_ll <- lmtest::bptest(m_loglog)
dw_kn <- lmtest::dwtest(m_kn)
dw_ll <- lmtest::dwtest(m_loglog)

table4_7c <- tibble(
  Model = c("qu = KN", "ln qu ~ ln N"),
  BP_statistic = c(bp_kn$statistic, bp_ll$statistic),
  BP_p = c(bp_kn$p.value, bp_ll$p.value),
  DW_statistic = c(dw_kn$statistic, dw_ll$statistic),
  DW_p = c(dw_kn$p.value, dw_ll$p.value),
  n_cooks_above_cutoff = c(sum(cooks.distance(m_kn) > cook_cut),
                           sum(cooks.distance(m_loglog) > cook_cut))) %>%
  mutate(constant_variance_at_5pct = ifelse(BP_p > 0.05, "Yes", "No"),
         no_autocorrelation_at_5pct = ifelse(DW_p > 0.05, "Yes", "No"),
         across(c(BP_statistic, DW_statistic), ~ round(.x, 3)),
         across(c(BP_p, DW_p), ~ signif(.x, 3)))

write_csv(table4_7c, "outputs/tables/Table4_7c_diagnostic_tests.csv")

# ---- Table 4.7d: influential points ----------------------------------------
table4_7d <- bind_rows(
  diag_data(m_kn) %>% mutate(Model = "qu = KN"),
  diag_data(m_loglog) %>% mutate(Model = "ln qu ~ ln N")
) %>%
  filter(influential) %>%
  left_join(clean, by = "sample_id") %>%
  select(Model, sample_id, depth_m, N, qu_kPa, std_resid, leverage, cooks) %>%
  mutate(across(c(std_resid, leverage, cooks), ~ round(.x, 4))) %>%
  arrange(Model, desc(cooks))

write_csv(table4_7d, "outputs/tables/Table4_7d_influential_points.csv")

# ---- Table 4.7b: sensitivity refits of qu = KN -----------------------------
dup_ids   <- clean$sample_id[duplicated(clean %>% select(-sample_id))]
stiff_ids <- c(123, 128, 132, 133, 134)
cook_ids  <- clean$sample_id[cooks.distance(m_kn) > cook_cut]

# Sanity check that the stiff rows are the 9.20 m samples
stopifnot(all(clean$depth_m[stiff_ids] == 9.20))

cases <- list(
  "(1) All samples"                    = integer(0),
  "(2) Without 4 duplicated rows"      = dup_ids,
  "(3) Without 5 stiff samples (9.20 m)" = stiff_ids,
  "(4) Without Cook's D > 4/n"         = cook_ids
)

table4_7b <- imap_dfr(cases, function(drop, name) {
  d <- filter(clean, !sample_id %in% drop)
  m <- lm(qu_kPa ~ 0 + N, data = d)
  ci <- lmtest::coefci(m, vcov. = hc3(m))
  res <- residuals(m)
  tibble(Case = name,
         n = nrow(d),
         Removed = if (length(drop)) paste(drop, collapse = ", ") else "none",
         K = coef(m)[["N"]],
         K_CI95_low_HC3 = ci[1, 1],
         K_CI95_high_HC3 = ci[1, 2],
         SE_HC3 = sqrt(hc3(m)[1, 1]),
         R2 = 1 - sum(res^2) / sum((d$qu_kPa - mean(d$qu_kPa))^2),
         RMSE = sqrt(mean(res^2)))
}) %>%
  mutate(K_change_pct = 100 * (K - K[1]) / K[1],
         across(c(K, K_CI95_low_HC3, K_CI95_high_HC3, SE_HC3, RMSE, K_change_pct),
                ~ round(.x, 3)),
         R2 = round(R2, 4))

write_csv(table4_7b, "outputs/tables/Table4_7b_sensitivity.csv")

# ---- Plain-English summary -------------------------------------------------
fmt_p <- function(p) format.pval(p, digits = 3, eps = 1e-4)

cat("\n--- Table 4.7c: diagnostic tests ---\n")
print(as.data.frame(table4_7c), row.names = FALSE)
cat("\n--- Table 4.7d: influential points (Cook's D > 4/n =", round(cook_cut, 4), ") ---\n")
print(as.data.frame(table4_7d), row.names = FALSE)
cat("\n--- Table 4.7b: sensitivity of K in qu = KN ---\n")
print(as.data.frame(table4_7b %>% select(-Removed)), row.names = FALSE)

cat("\nSummary:\n")
for (i in 1:2) {
  r <- table4_7c[i, ]
  cat(sprintf("- %s: Breusch-Pagan p = %s (%s); Durbin-Watson = %.2f, p = %s (%s); %d points above the Cook's cut-off.\n",
              r$Model, fmt_p(r$BP_p),
              ifelse(r$BP_p > 0.05, "constant spread", "spread changes with N"),
              r$DW_statistic, fmt_p(r$DW_p),
              ifelse(r$DW_p > 0.05, "no autocorrelation", "residuals autocorrelated in file order"),
              r$n_cooks_above_cutoff))
}
k_range <- range(table4_7b$K)
cat(sprintf("- K ranges from %.2f to %.2f kPa per blow across the four cases (largest change %.1f%% from the full-data value of %.2f).\n",
            k_range[1], k_range[2], max(abs(table4_7b$K_change_pct)), table4_7b$K[1]))
overlap <- all(table4_7b$K >= table4_7b$K_CI95_low_HC3[1] &
               table4_7b$K <= table4_7b$K_CI95_high_HC3[1])
cat(sprintf("- %s\n", ifelse(overlap,
  "Every refitted K lies inside the full-data robust 95% CI, so K is stable.",
  "At least one refitted K lies outside the full-data robust 95% CI.")))
cat("Saved Fig4_7a_diagnostics_KN.png, Fig4_7b_diagnostics_loglog.png, Table4_7b_sensitivity.csv,",
    "Table4_7c_diagnostic_tests.csv, Table4_7d_influential_points.csv\n")

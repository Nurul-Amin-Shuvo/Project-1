# 04_correlation.R
# Pearson and Spearman correlations with p-values (Table 4.4),
# correlation heatmap (Fig 4.4) and scatter-plot matrix (Fig 4.5).
# All 134 samples are used, including the 4 duplicated rows.

source("R/00_packages.R")
set.seed(123)

clean <- readRDS("data/clean.rds")

var_labels <- c(depth_m = "Depth (m)",
                N       = "SPT-N (blows/300 mm)",
                qu_kPa  = "qu (kPa)")

theme_set(theme_bw(base_size = 11, base_family = "Times New Roman"))
save_fig <- function(plot, file, width, height) {
  ggsave(file.path("outputs/figures", file), plot, width = width,
         height = height, units = "in", dpi = 300, device = ragg::agg_png)
}

sig_stars <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
}

# ---- Table 4.4: pairwise correlations --------------------------------------
# Spearman p-values use the normal approximation (exact = FALSE) because of ties
pairs <- combn(names(var_labels), 2, simplify = FALSE)

table4_4 <- map_dfr(pairs, function(v) {
  x <- clean[[v[1]]]
  y <- clean[[v[2]]]
  pear <- cor.test(x, y, method = "pearson")
  spear <- cor.test(x, y, method = "spearman", exact = FALSE)
  tibble(Variable_1   = var_labels[[v[1]]],
         Variable_2   = var_labels[[v[2]]],
         n            = length(x),
         Pearson_r    = round(unname(pear$estimate), 3),
         Pearson_CI95 = sprintf("%.3f to %.3f", pear$conf.int[1], pear$conf.int[2]),
         Pearson_p    = signif(pear$p.value, 3),
         Spearman_rho = round(unname(spear$estimate), 3),
         Spearman_p   = signif(spear$p.value, 3))
})

write_csv(table4_4, "outputs/tables/Table4_4_correlation.csv")

# ---- Fig 4.4: correlation heatmap (Pearson and Spearman) -------------------
cor_long <- expand_grid(method = c("Pearson", "Spearman"),
                        row = names(var_labels), col = names(var_labels)) %>%
  mutate(test = pmap(list(method, row, col), function(m, a, b) {
    if (a == b) return(list(r = 1, p = NA_real_))
    ct <- cor.test(clean[[a]], clean[[b]], method = tolower(m), exact = FALSE)
    list(r = unname(ct$estimate), p = ct$p.value)
  }),
  r = map_dbl(test, "r"),
  p = map_dbl(test, "p"),
  row = factor(var_labels[row], levels = var_labels),
  col = factor(var_labels[col], levels = rev(var_labels)),
  label = ifelse(is.na(p), sprintf("%.2f", r),
                 paste0(sprintf("%.2f", r), sig_stars(p))))

fig4_4 <- ggplot(cor_long, aes(x = row, y = col, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = label), family = "Times New Roman", size = 4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  facet_wrap(~ method) +
  coord_equal() +
  labs(x = NULL, y = NULL,
       caption = "* p < 0.05, ** p < 0.01, *** p < 0.001 (n = 134)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())

save_fig(fig4_4, "Fig4_4_correlation_heatmap.png", width = 8, height = 4.2)

# ---- Fig 4.5: scatter-plot matrix ------------------------------------------
# Upper panels show both Pearson r and Spearman rho
cor_both <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  r   <- cor(x, y, method = "pearson")
  rho <- cor(x, y, method = "spearman")
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, family = "Times New Roman", size = 4.5,
             label = sprintf("Pearson r = %.2f\nSpearman ρ = %.2f", r, rho)) +
    xlim(0, 1) + ylim(0, 1) +
    theme(panel.grid = element_blank())
}

points_lm <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_point(shape = 21, fill = "grey60", colour = "black", size = 1.3) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "black", linewidth = 0.6)
}

fig4_5 <- GGally::ggpairs(
  clean %>% rename(!!!setNames(names(var_labels), var_labels)),
  columns = 1:3,
  upper = list(continuous = cor_both),
  lower = list(continuous = points_lm),
  diag  = list(continuous = GGally::wrap("densityDiag", colour = "black"))
)

save_fig(fig4_5, "Fig4_5_scatter_matrix.png", width = 7.5, height = 7)

# ---- Plain-English summary -------------------------------------------------
strength <- function(r) {
  a <- abs(r)
  if (a < 0.1) "negligible" else if (a < 0.3) "weak" else if (a < 0.5) "moderate"
  else if (a < 0.7) "strong" else "very strong"
}

cat("\n--- Table 4.4: correlations (n = 134) ---\n")
print(as.data.frame(table4_4), row.names = FALSE)

cat("\nSummary:\n")
for (i in seq_len(nrow(table4_4))) {
  r <- table4_4[i, ]
  cat(sprintf("- %s vs %s: Pearson r = %.2f (%s, p = %s), Spearman rho = %.2f (%s, p = %s).\n",
              r$Variable_1, r$Variable_2,
              r$Pearson_r, strength(r$Pearson_r), format.pval(r$Pearson_p, digits = 3, eps = 1e-4),
              r$Spearman_rho, strength(r$Spearman_rho), format.pval(r$Spearman_p, digits = 3, eps = 1e-4)))
}
cat("Strength bands: <0.1 negligible, 0.1-0.3 weak, 0.3-0.5 moderate, 0.5-0.7 strong, >=0.7 very strong.\n")
cat("Saved Table4_4_correlation.csv, Fig4_4_correlation_heatmap.png, Fig4_5_scatter_matrix.png\n")

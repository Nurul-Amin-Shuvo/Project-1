# Project 1 results summary

**Statistical Analysis of SPT-N Values and Shear Strength Parameters in Dhaka**

134 samples (field SPT-N, depth and unconfined compressive strength qu) from Kawlar, Ashkona, Purbachal, 300 ft road and MRT Line 6 (Uttara to Motijheel). All analyses use all 134 rows. The 4 duplicated rows and the 5 stiff samples at 9.20 m are kept, and are removed only in the sensitivity checks (Table 4.7b).

Reproduce everything with `Rscript R/99_run_all.R` (about 30 s). The R version and package versions are in `outputs/sessionInfo.txt`.

## Key findings

- **Recommended equation: qu = 11.6 N** (kPa, field SPT-N). Robust 95% CI for K: 10.7–12.5. Bootstrap 95% CI: 10.7–12.4.
- **Cross-validation winner:** the polynomial has the lowest CV RMSE (34.7 kPa, against 37.4 kPa for qu = KN). It beats KN in only 64% of folds, and its advantage almost disappears (0.6 kPa) when the 5 stiff samples at 9.20 m are left out of scoring. **qu = KN comes second and is the recommended model.**
- Depth adds no reliable information once SPT-N is known.
- Terzaghi & Peck (1967), qu = 12.5N, is the closest published correlation. Serajuddin & Chowdhury (1998) over-predicts by about 51 kPa (+48%).

## Tables (`outputs/tables`)

| File | Interpretation |
|---|---|
| `Table4_1_data_summary.csv` | All 134 samples are complete, with no missing values; 4 rows are exact duplicates and are kept.<br>Depth 1.95–61.45 m, SPT-N 2–41, qu 22.5–589.5 kPa.<br>Mean qu 119 kPa, against a median of 93 kPa, already shows a right skew. |
| `Table4_2_descriptive.csv` | All three variables are strongly right-skewed (skewness 1.5–2.9), with CV of 65–78%.<br>Depth has heavy tails (excess kurtosis 14.6) because of a few deep samples (27–61 m).<br>Median SPT-N 8 and median qu 93 kPa: soft-to-medium fine-grained soils dominate. |
| `Table4_3_normality.csv` | Shapiro-Wilk rejects normality for raw depth, N and qu (all p < 1e-9).<br>ln qu passes (p = 0.053). ln N narrowly fails (p = 0.029), mainly because N values are repeated whole numbers.<br>This supports log-log models and Spearman correlations alongside Pearson. |
| `Table4_3b_outliers.csv` | 16 of 134 samples fall outside the 1.5 IQR limits in at least one variable (depth 4, N 9, qu 12).<br>The largest are the 5 stiff samples at 9.20 m (N 28–41, qu 401–589 kPa).<br>All are flagged, not removed. |
| `Table4_4_correlation.csv` | SPT-N vs qu is very strong: Pearson r = 0.91 (95% CI 0.88–0.94), Spearman ρ = 0.89, both p < 1e-46.<br>Depth vs N is moderate (r ≈ 0.33). Depth vs qu is weak and depends on the method (r = 0.16, p = 0.058; ρ = 0.30, p < 0.001).<br>That disagreement comes from the high-qu stiff samples at a shallow depth. |
| `Table4_5_model_comparison.csv` | On the full data, the polynomial fits best (R² 0.855, RMSE 35.4 kPa, AIC 1344).<br>qu = 11.57N fits almost as well as linear (R² 0.830, RMSE 38.3 kPa) and has the lowest BIC after the polynomial. The linear intercept is not significant (p = 0.34).<br>The nls power exponent is 1.03, close to 1, so qu rises roughly in proportion to N. The logarithmic model is worst (R² 0.655). |
| `Table4_5b_breusch_pagan.csv` | The linear model has uneven spread (BP = 24.9, p = 6e-7): errors grow with N.<br>The log-log model has constant spread (p = 0.40).<br>So qu-scale models need robust (HC3) standard errors for valid tests and intervals. |
| `Table4_6_mlr_coefficients.csv` | In qu ~ N + depth, the depth term is −1.98 kPa/m (robust p = 0.012). It is negative, which is not physically sensible.<br>With robust SEs, the N × depth interaction is not significant (p = 0.097), and ln depth in the log model is not significant (p = 0.23).<br>Robust SEs are 1.5–2.6 times the ordinary ones for the qu-scale models. |
| `Table4_6b_mlr_fit.csv` | Adding depth lowers RMSE from 38.2 to 35.5 kPa. The interaction model reaches 33.2 kPa.<br>The log-log models have the highest RMSE but are the only ones with constant spread.<br>The fit gains come mainly from the 5 stiff samples (see Table 4.7). |
| `Table4_6c_vif.csv` | All VIFs are 1.13–1.66 (the interaction model is centred first).<br>N and depth are only weakly related, so overlap between predictors (collinearity) is not a concern.<br>The depth coefficients are therefore not distorted by collinearity. |
| `Table4_7_nested_F_test.csv` | Ordinary partial F-tests say depth helps (p < 1e-4), but they assume constant spread.<br>With robust Wald tests, only the simple depth term survives (p = 0.012). The interaction (p = 0.097) and ln depth (p = 0.23) do not.<br>Without the 5 stiff samples, depth is not significant (robust p = 0.15): depth does not reliably improve prediction. |
| `Table4_7b_sensitivity.csv` | K = 11.57 (all 134), 11.59 (without duplicates), 10.76 (without stiff samples), 11.32 (without Cook's outliers).<br>The largest change is −7%, and every refit lies inside the full-data robust 95% CI.<br>K is stable. Report qu = 11.6N with a sensitivity range of 10.8–11.6. |
| `Table4_7c_diagnostic_tests.csv` | qu = KN fails Breusch-Pagan (p = 4e-6) and Durbin-Watson (DW = 1.41, p = 0.0003) when the data are in file order.<br>The log-log model passes both (BP p = 0.40, DW p = 0.15).<br>The DW failure most likely reflects samples grouped by site and the leftover curvature, not time-series correlation. |
| `Table4_7d_influential_points.csv` | 11 points exceed Cook's D > 4/n for KN. Sample 134 (N 41, qu 589) has D = 0.83, and samples 10 and 11 (deep, N = 26.5, low qu) are next.<br>The log-log model has 13 points above the cut-off, but all are small (max D = 0.07).<br>No single sample controls the log-log fit. |
| `Table4_8_cv_results.csv` | Repeated 10-fold CV (10 × 10 folds) RMSE: polynomial 34.7, KN 37.4, linear 37.8, nls 38.4, log-log 40.4 kPa.<br>The polynomial beats KN in 64% of folds, by 2.8 kPa. That is small next to the fold-to-fold SD of about 12 kPa, and it comes from the stiff samples.<br>**CV winner: polynomial, narrowly. Recommended: qu = KN (second, with 1 coefficient).** |
| `Table4_8b_train_test_split.csv` | 70/30 split (93/41 samples): test RMSE 41.2 (log-log) to 45.2 kPa (nls). For KN, test R² is 0.659 against 0.866 on training.<br>A single split ranks the models differently from CV, which shows how much it depends on which samples land in the test set.<br>Use the CV figures; this split is for comparison with the bachelor work only. |
| `Table4_8c_bootstrap_K.csv` | Bootstrap (2000 resamples) 95% percentile CI for K: 10.66–12.36. HC3 robust: 10.68–12.47.<br>The two agree closely. The ordinary CI (11.05–12.09) is about 40% too narrow.<br>This confirms that robust or bootstrap intervals should be reported. |
| `Table4_9_published_comparison.csv` | Best published correlation: Terzaghi & Peck (1967) qu = 12.5N (RMSE 40.0 kPa, bias +7.5 kPa, median ratio 1.02). Stroud PI < 20 (13N) is next.<br>Serajuddin & Chowdhury (1998), qu = 16.8N, over-predicts by 51 kPa on average (+48%, median ratio 1.37).<br>Sanglerat, Nixon, Sowers 12.5N and Hara over-predict by 76–185%. |

## Figures (`outputs/figures`)

| File | Interpretation |
|---|---|
| `Fig3_1_study_area.png` | Map of the 18 locations: 14 MRT Line 6 stations from Uttara to Motijheel, plus Kawlar, Ashkona, Purbachal and the 300 ft road.<br>Points are coloured by project (Aga Khan Academy, 100 ft Khal project, MRT Line 6), and the MRT alignment is drawn as a line.<br>Kawlar and Ashkona (Aga Khan Academy, confirmed) are placed at their OpenStreetMap locality positions in Dakshinkhan. **Check:** Purbachal and the 300 ft road (marked *) are approximate, and their project assignment is unconfirmed. |
| `Fig4_1_histograms.png` | Histograms with density curves for depth, SPT-N and qu.<br>All three are right-skewed, with long upper tails.<br>qu has a secondary peak near 230 kPa and a few values above 400 kPa. |
| `Fig4_2_boxplots.png` | Boxplots show the mean above the median for all three variables.<br>Outliers sit only on the high side.<br>The qu outliers above 250 kPa are mostly the stiff 9.20 m samples. |
| `Fig4_2b_depth_profiles.png` | Borehole-log-style profiles of N and qu against depth.<br>Both increase gently down to about 25 m, with large scatter at any given depth. The stiff cluster at 9.20 m stands out.<br>The trend line stops at 25 m because only 4 samples lie deeper. |
| `Fig4_3_qq_plots.png` | Q-Q plots of raw depth, N and qu curve strongly upward (right skew).<br>The ln N and ln qu plots are nearly straight.<br>This supports the log transform for qu. |
| `Fig4_4_correlation_heatmap.png` | Pearson and Spearman heatmaps side by side, with significance stars.<br>N–qu is very strong (0.91 / 0.89) under both methods.<br>Depth–qu is the only pair where the methods disagree (0.16 vs 0.30). |
| `Fig4_5_scatter_matrix.png` | Scatter-plot matrix with straight-line fits, density curves and both correlation coefficients.<br>The N–qu panel is tight and close to straight, with spread that grows with N.<br>The depth panels show weak, scattered relationships. |
| `Fig4_6_simple_models.png` | Six fitted curves over the data.<br>The linear, through-origin and nls power curves almost overlap. The polynomial bends up above N ≈ 30.<br>The logarithmic curve goes negative at N = 2, and the back-transformed log-log curve runs too low at high N. |
| `Fig4_7a_diagnostics_KN.png` | Residual diagnostics for qu = KN: the scatter widens as fitted qu increases.<br>A U-shaped residual trend is pulled up by samples 123, 128, 132, 133 and 134.<br>Heavy tails come from samples 21 (N 8, qu 269), 10 and 11. |
| `Fig4_7b_diagnostics_loglog.png` | Residual diagnostics for the log-log model: even spread and a near-straight Q-Q plot.<br>No single point has much influence.<br>Its statistical assumptions are the best satisfied of any model. |
| `Fig4_8_predicted_vs_measured.png` | Out-of-fold CV predictions against measured qu, with 1:1 and ±30% lines.<br>70% of predictions are within ±30% for KN, and 77% for the polynomial.<br>All models predict too low for the highest qu values and too high for deep samples 10 and 11. |
| `Fig4_9_KN_bands.png` | qu = 11.57N with its robust 95% confidence band and a 95% prediction band from the log-log model.<br>The prediction band widens with N, matching the uneven spread.<br>96% of samples fall inside it, close to the nominal 95%. |
| `Fig4_10_published_correlations.png` | Published correlations drawn over the data, with Stroud and Sowers as bands.<br>The data follow This study and Terzaghi & Peck, inside the Stroud band.<br>Sanglerat, Nixon and Hara lie far above all the data. |
| `Fig4_11_published_measured_vs_predicted.png` | Measured vs predicted panels for all 15 correlations, ordered by RMSE.<br>The best panels cluster along the 1:1 line. The Sivrikaya & Toğrol and Stroud PI > 30 panels sit below it.<br>Serajuddin & Chowdhury and the worse panels lie consistently above the 1:1 line. |

## Caveats

- Published correlation coefficients were used as supplied, without checking them against the original papers. Several were developed with N60, and this data uses field N.
- No PI data: Stroud and Sowers variants are evaluated separately.
- Fig 3.1: Purbachal and 300 ft road locations and their project assignment need confirming. All non-MRT points are locality positions, not exact borehole locations.

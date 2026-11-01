# Final Project: Monthly Atmospheric CO₂ at Mauna Loa (1958–2026)

STA6938 Time Series. Data: NOAA Global Monitoring Laboratory, 822 monthly means, March 1958 – August 2026.

**Final model:** ARIMA(0,1,1)×(0,1,1)₁₂, θ̂₁ = 0.389, Θ̂₁ = 0.870. **Headline forecast:** CO₂ first exceeds 440 ppm in April 2030.

## Folder contents

| Path | What it is |
|---|---|
| `report/main.pdf` | **The report** (5 pages; text, values and interpretation only, with no tables or figures, as the professor requires) |
| `report/main.tex`, `report/references.bib` | LaTeX source of the report |
| `CO2_Report_LaTeX.zip` | The LaTeX report as a zip, ready to upload to Overleaf |
| `figures/` | 13 figures for the slides: PNG (300 dpi) and PDF (vector) |
| `tables/` | 15 tables for the slides: CSV (Excel/PowerPoint) and `.tex` (Beamer, `booktabs`) |
| `R/co2_analysis.R` | The full analysis; regenerates every figure, table and number |
| `results/key_values.txt` | Every number quoted in the report |
| `data/CO2_MaunaLoa_monthly.csv` | The data |

## Figures (suggested slide order)

| File | Shows | Report section |
|---|---|---|
| `fig01_co2_series` | The series with its 2×12 moving-average trend | 2. Description |
| `fig02_classical_decomposition` | Observed, trend, seasonal and remainder components | 2. Description |
| `fig03_seasonal_indices` | Seasonal index by month (peak May, trough October) | 2. Description |
| `fig04_seasonal_amplitude_by_year` | Seasonal amplitude by year (why no transformation) | 3. Transformation |
| `fig05_acf_original_and_differences` | ACF of X, ∇X and ∇₁₂X (why both differences) | 3. Stationarity |
| `fig06_annual_growth_lag12_difference` | 12-month change and decade averages (accelerating growth) | 2–3 |
| `fig07_regular_and_seasonal_difference` | The stationary series ∇∇₁₂X | 3. Stationarity |
| `fig08_acf_pacf_differenced` | ACF and PACF of ∇∇₁₂X (identification) | 4. Model building |
| `fig09_regression_residuals` | Residuals of the trend + dummies regression (estimation vs differencing) | 3. Stationarity |
| `fig10_final_model_residual_diagnostics` | Residuals, residual ACF, histogram, Q–Q plot | 6. Adequacy |
| `fig11_ljung_box_pvalues` | Ljung–Box p-values by lag, final vs alternative model | 6. Adequacy |
| `fig12_holdout_forecast_comparison` | 36-month hold-out: ARIMA vs Holt–Winters vs regression | 7. Forecasting |
| `fig13_forecast_to_2030` | Forecast to December 2030 with 80% and 95% intervals | 7. Forecasting |

## Tables

| File | Contents |
|---|---|
| `tab01_descriptive_statistics` | n, range, mean, SD, min/max, average growth |
| `tab02_growth_by_decade` | Mean level, annual increase and seasonal amplitude by decade |
| `tab03_seasonal_indices` | Monthly seasonal indices |
| `tab04_adf_tests` | ADF tests (single mean and trend, lags 0–2) |
| `tab05_regression_trend_seasonal` | Quadratic trend + seasonal dummy regression |
| `tab06_white_noise_differenced` | Ljung–Box white-noise tests on ∇∇₁₂X |
| `tab07_candidate_models` | Nine candidate models with AIC, SBC and Ljung–Box |
| `tab08_parameter_estimates` | ML and CLS estimates, z-tests and 95% CIs |
| `tab09_residual_diagnostics` | Independence, constant variance and normality tests |
| `tab10_outlying_months` | The five months with standardized residuals beyond ±3 |
| `tab11_sensitivity_estimates` | Estimates with the outlying months set to missing |
| `tab12_holt_winters_parameters` | Holt–Winters α, β, γ, level and slope |
| `tab13_holdout_accuracy` | Hold-out RMSE, MAE, MAPE and interval coverage |
| `tab14_forecasts_selected` / `tab14_forecasts_full` | Forecasts and prediction intervals to December 2030 |
| `tab15_regression_autocorrelated_errors` | Regression with ARMA errors (Note 12): AR coefficients near a unit root |

To use a table in Beamer, put `\usepackage{booktabs}` and `\usepackage[utf8]{inputenc}` in the preamble and `\input{tables/tab13_holdout_accuracy.tex}` on the slide.

## Reproducing everything

```bash
cd "CO2 Project"
LC_ALL=C.UTF-8 Rscript R/co2_analysis.R      # figures/, tables/, results/ (UTF-8 locale for the × symbol)
cd report && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

R packages: `forecast`, `tseries`, `urca`, `nortest`, `ggplot2`, `patchwork`, `scales`, `kableExtra`, `ggrepel`.

## Before submitting

- Replace **Student Name** in `report/main.tex` (title block) with your name.
- Check the instructor line and the date (December 6, 2026) in the title block.

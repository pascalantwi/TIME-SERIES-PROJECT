# STA6938 Time Series — Final Project Report

**Student:** _To be completed_
**Dataset:** CSX Corporation daily stock prices, January 3, 2000 – October 19, 2026 (`CSX_daily.csv`, 6,739 trading days)
**Software used:** R 4.3.3 (`forecast`, `tseries`, `urca`, `nortest`, `ggplot2`, `patchwork`, `kableExtra`)

| Deliverable | Weight | Due | Status |
|---|---|---|---|
| Project report (≤ 5 pages, excluding output and graphs) | 15% | December 6, 11:59 PM (Canvas) | Draft ready: `CSX_Time_Series_Report.pdf` |
| Project presentation (≤ 15 minutes, slides) | 5% | December 3, in class | To be completed |

## Files

| File | What it is |
|---|---|
| `CSX_Time_Series_Report.pdf` | The report (text, 8 tables, 7 figures, R code appendix) |
| `CSX_Time_Series_Report.Rmd` | R Markdown source that reproduces the whole report from the data |
| `CSX_daily.csv` | The dataset |
| `figures/` | Every figure as a 300-dpi PNG (for slides) and as a vector PDF |
| Note PDFs and `project.pdf` | Course notes and the assignment sheet |

## Report sections (all required items covered)

1. **Topic and Motivation** — CSX, a major U.S. freight railroad, and why its stock is interesting.
2. **Description of the Series** — data source, adjusted close, Table 1 descriptive statistics, Figure 1.
3. **Stationarity of the Series** — log transformation, plot and ACF, ADF tests at lags 0–2 (Table 2), estimation vs. differencing (Note 4), white-noise tests (Table 3).
4. **Model Selection and Model-Building Procedure** — ACF/PACF identification (Figure 3), ARMA candidates by AIC/SBC (Table 4), initial ARIMA(0,1,1) with drift.
5. **Hypothesis Tests and Confidence Intervals** — parameter z-tests, 95% CIs, invertibility (Table 5).
6. **Model Adequacy Checking and Remedial Methods** — residual diagnostics; remedy: subset ARIMA(0,1,6) with MA lags 1 and 6 (Table 6, Table 7, Figures 4–6).
7. **Forecasting** — out-of-sample comparison with random walk, SES and Holt (Table 8); 20-day price forecasts with 80%/95% prediction intervals (Table 9, Figure 7).
8. **Conclusions** — findings and limitations.

## To be completed

- Replace `Student Name` in the first lines of `CSX_Time_Series_Report.Rmd` with your name.
- Personalize Section 1 (Topic and Motivation) with your own reason for choosing CSX.
- Confirm the data source (the file format matches Yahoo Finance downloads).
- Presentation slides.

## Rebuilding the report

In R or RStudio, with the working directory set to this folder:

```r
rmarkdown::render("CSX_Time_Series_Report.Rmd")
```

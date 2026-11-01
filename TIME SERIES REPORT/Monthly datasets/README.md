# Monthly Datasets for the Final Project

Five real, public monthly time series that fit the project requirements in `project.pdf`. Each is a large sample from an official source, and none is used in the course notes. Every file has two columns: `Month` (YYYY-MM) and the value.

| File | Series | Months | Period | Trend | Seasonality | Source |
|---|---|---|---|---|---|---|
| `CO2_MaunaLoa_monthly.csv` | Atmospheric CO2 at Mauna Loa (ppm) | 822 | Mar 1958 – Aug 2026 | Strong, upward | Strong (12-month cycle) | NOAA Global Monitoring Laboratory (listed on Data.gov) |
| `US_CPI_monthly.csv` | U.S. Consumer Price Index, all urban consumers (1982–84 = 100) | 1,363 | Jan 1913 – Aug 2026 | Strong, upward | Weak | U.S. Bureau of Labor Statistics |
| `Global_Temperature_Anomaly_NASA_monthly.csv` | Global surface temperature anomaly (°C) | 1,752 | Jan 1880 – Dec 2025 | Upward since ~1970 | None (anomalies are deseasonalized) | NASA GISS (GISTEMP) |
| `Natural_Gas_HenryHub_monthly.csv` | Henry Hub natural gas spot price (USD per MMBtu) | 357 | Jan 1997 – Sep 2026 | None (wanders) | Weak | U.S. Energy Information Administration |
| `WTI_Crude_Oil_monthly.csv` | WTI crude oil spot price (USD per barrel) | 489 | Jan 1986 – Sep 2026 | None (wanders) | Weak | U.S. Energy Information Administration |

The files were downloaded from the [GitHub "datasets" collection](https://github.com/datasets), which republishes the official agency data. No values are missing, except that `US_CPI_monthly.csv` has no October 2025 row; to my knowledge the BLS did not publish that month because of the 2025 government shutdown. If you use the CPI series, fill that month in, e.g., as the average of September and November 2025.

## Which course methods each series supports

| Series | Transformation (Note 4) | Decomposition / regression (Note 4) | Differencing (Note 4) | Exponential smoothing (Note 5) | AR / MA / ARMA / ARIMA (Notes 7–9) |
|---|---|---|---|---|---|
| CO2 | Not needed | Yes: trend + seasonal | Lag 1 and lag 12 | Holt-Winters | ARIMA on the differenced series |
| US CPI | Log | Trend only | Log, then difference (inflation rate) | Holt | ARIMA on monthly inflation |
| Global temperature | Not needed | Trend only | Lag 1 | Holt | ARIMA(p,1,q) |
| Natural gas | Log | No clear components | Lag 1 | Simple ES | ARIMA on log returns |
| WTI oil | Log | No clear components | Lag 1 | Simple ES | ARIMA on log returns |

## Recommendation

1. **CO2 at Mauna Loa (best fit for the course).** It has a clear trend and a clear 12-month seasonal pattern, so it uses almost every tool in the notes: classical decomposition, regression on trend and seasonal dummies, lag-1 and lag-12 differencing, Holt-Winters smoothing, then ACF/PACF identification and an ARIMA model for the differenced series. With 822 months it is a large sample, and the residuals are much better behaved than stock returns.
2. **U.S. CPI / inflation.** It covers more than a century of data and is easy to motivate (the cost of living). Taking logs and differencing gives the monthly inflation rate, which suits AR and ARMA models well.
3. **Natural gas price.** It has a good economic story (energy prices and weather), but like the CSX stock price it behaves almost like a random walk.

## Not downloadable from this session, but a good option

**U.S. Rail Freight Carloads**, monthly from January 2000. It is published by the Bureau of Transportation Statistics, listed on [Data.gov](https://catalog.data.gov/dataset/freight-rail-traffic-carloads) and available on [FRED](https://fred.stlouisfed.org/series/RAILFRTCARLOADS). It links directly to the CSX topic and has seasonality. On FRED, click **Download → CSV** and upload the file to this folder. If FRED offers both versions, choose the *not seasonally adjusted* one, so that the seasonal pattern is still in the data for you to analyze.

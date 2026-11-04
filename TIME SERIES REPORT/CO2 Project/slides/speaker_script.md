# Speaker script: Modeling and Forecasting Monthly Atmospheric CO2 at Mauna Loa

Pascal Antwi, STA6938 Time Series. About 15 minutes; timing cues in brackets.

## Slide 1: Title

[0:00 -- 0:20] Good morning. My name is Pascal Antwi. My project is titled Modeling and Forecasting Monthly Atmospheric CO2 at Mauna Loa, 1958 to 2026. I will show how the series was made stationary, how a seasonal ARIMA model was identified, tested and validated, and what it forecasts up to 2030.

## Slide 2: Roadmap

[0:20 - 0:50] The presentation follows the Box-Jenkins approach. I begin with the background and the data, then examine trend, seasonality and stationarity. Next I identify and select a model, report the estimates with hypothesis tests and confidence intervals, check model adequacy and apply remedial measures, and finally evaluate and present the forecasts. Three questions guide the analysis: how fast CO2 is rising, whether the rise is accelerating, and what levels should be expected by 2030.

## Slide 3: Background and motivation

[0:50 - 2:00] Carbon dioxide is the principal driver of anthropogenic climate change, and international climate targets are framed around how fast its concentration rises and when it crosses particular thresholds, such as 440 parts per million. A rigorous statistical characterization and forecast of the Mauna Loa record is therefore of direct scientific and policy relevance. The series is also well suited to time series methodology, because it combines a stochastic trend, a pronounced seasonal cycle and short-term serial dependence. The question is not whether CO2 is rising, which is well established, but how fast it is rising, whether the rise is accelerating, and what levels to expect, each with proper statistical uncertainty. The four objectives on the right follow directly from that question.

## Slide 4: The data: Mauna Loa monthly CO2

[2:00 - 2:40] The data are the monthly mean CO2 concentrations published by the NOAA Global Monitoring Laboratory, which are freely available online. The record contains 822 consecutive months, from March 1958 to August 2026, with no missing values, which is a large sample for time series modeling. Over this period the concentration rose from 315.7 to 427.6 parts per million, an increase of about 112 parts per million, or 35 percent.

## Slide 5: Step 1: the time plot

[2:40 - 3:20] Following standard practice, the first step was the time plot. Two features are immediately visible. First, there is a persistent upward trend, so the mean of the series depends on time. Second, there is a regular twelve-month seasonal cycle with almost the same shape every year. The series is therefore clearly non-stationary, and both the trend and the seasonal component must be addressed before an ARMA-type model can be fitted.

## Slide 6: Decomposition: trend plus a stable seasonal cycle

[3:20 - 4:15] A classical additive decomposition with a two-by-twelve moving-average filter separates these components. The trend accounts for about 97 percent of the variance. The seasonal component peaks in May and reaches its minimum in October, with an amplitude of about 6.3 parts per million, reflecting the uptake of carbon by Northern Hemisphere vegetation during the growing season. An important modeling decision follows from this panel: the seasonal amplitude rose only from about 6.0 to 7.0 parts per million while the level rose by 35 percent. The seasonal swings do not grow in proportion to the level, so an additive model is appropriate and no log or Box-Cox transformation is required.

## Slide 7: The growth rate is accelerating

[4:15 - 5:00] This figure shows the twelve-month change, which is the annual growth in concentration. The decade averages rise from 0.86 parts per million per year in the 1960s to 2.54 in the 2020s, so the growth rate has roughly tripled. In other words, CO2 is not only increasing but increasing faster. The largest twelve-month increases occurred in 2016 and 2024, both strong El Nino years, which will matter later when we examine the residuals and the forecasts.

## Slide 8: Achieving stationarity

[5:00 - 6:00] The sample ACF confirms the non-stationarity. For the original series it decays extremely slowly, from 0.996, and every autocorrelation stays far outside the shaded 95 percent confidence band of plus or minus 0.068, the classic signature of a unit root. After a first difference, the autocorrelations at the seasonal lags 12, 24 and 36 remain around 0.9, which shows that a seasonal difference is also needed. After the seasonal difference alone, the ACF still decays slowly, so a regular difference is needed as well. The Augmented Dickey-Fuller test agrees: for the seasonally adjusted series the unit-root hypothesis cannot be rejected, whereas after both differences it is rejected decisively. The stationary series used for modeling is therefore W sub t, the regular and seasonal difference of the data.

## Slide 9: Estimation versus differencing

[6:00 - 6:45] As an alternative to differencing, I estimated the trend and seasonality deterministically, by regressing the series on a quadratic trend and eleven monthly dummy variables. The fit appears excellent, with an R-squared of 0.9994. However, the residuals are far from white noise: their lag-one autocorrelation is 0.915 and the Durbin-Watson statistic is 0.17, far below the value of 2 expected under independence. Even when the errors are modeled as an ARMA process, the autoregressive root is essentially one. The deterministic approach is therefore inadequate, and differencing is the appropriate treatment.

## Slide 10: Identification from the ACF and PACF of W_t

[6:45 - 7:45] The Ljung-Box test shows that the differenced series is not white noise, so its dependence must be modeled. I identify the orders from the overall pattern of the two plots, using the shaded 95 percent confidence band of plus or minus 0.069 as a guide to individual lags, and the Ljung-Box test for joint significance. The ACF is dominated by two large spikes, at lag 1 and the seasonal lag 12, and then drops to near zero, so it cuts off. The two borderline values at lags 14 and 23 are isolated, and about two or three such values are expected by chance among 48 lags, so I do not treat them as structure. The PACF instead decays gradually, both over lags 1 to 3 and across the seasonal lags 12, 24 and 36, so it tails off. An ACF that cuts off together with a PACF that tails off is the signature of moving-average behavior, here at both the non-seasonal and the seasonal level. This identifies the multiplicative seasonal ARIMA model with one regular and one seasonal moving-average term, shown in the box.

## Slide 11: Model selection: AIC and SBC

[7:45 - 8:30] To confirm the identification, nine candidate models were estimated by maximum likelihood, including models with additional autoregressive, moving-average and seasonal autoregressive terms. The identified two-parameter model has the lowest Schwarz Bayesian criterion. One richer model has a slightly lower AIC, but SBC penalizes additional parameters more heavily, which is appropriate in a large sample. Pure autoregressive and seasonal autoregressive alternatives fit clearly worse, consistent with the identification on the previous slide.

## Slide 12: Four model families compared

[8:30 - 9:30] I also compared four model families: ARMA, ARIMA, regression with autocorrelated errors, and seasonal ARIMA, choosing the best order within each family by SBC. Only the seasonal ARIMA leaves white-noise residuals; the other three retain significant autocorrelation at lag 12. The best non-seasonal ARIMA even places its autoregressive roots exactly on the unit circle, at periods of 6 and 12 months, in effect imitating a seasonal difference. The seasonal ARIMA captures the full structure with only two parameters.

## Slide 13: Estimates, hypothesis tests and confidence intervals

[9:30 - 10:30] This slide reports the hypothesis tests and confidence intervals. The non-seasonal coefficient theta is 0.389 and the seasonal coefficient capital theta is 0.870. Both are highly significant, with z-statistics of 10.6 and 51.3, and their 95 percent confidence intervals exclude zero. Conditional least squares gives very similar values, so the conclusions do not depend on the estimation method. All roots lie outside the unit circle, so the model is invertible and suitable for forecasting. In interpretation, about 39 percent of each monthly shock is reversed the following month, and the large seasonal coefficient indicates a very stable seasonal pattern. Because its confidence interval excludes one, seasonal differencing is justified.

## Slide 14: Model adequacy: residual diagnostics

[10:30 - 11:20] Model adequacy was checked on the residuals. The residual ACF shows no systematic pattern; only lag 3 falls just outside the 95 percent confidence band, fewer than the number expected by chance over 36 lags. The residuals also pass the Ljung-Box test at lags 6, 12, 18 and 24, so they are consistent with white noise, and the Ljung-Box test on the squared residuals shows no change in variance over time. Normality is only marginally rejected by the Shapiro-Wilk test, at a p-value of 0.048. This is driven by five months with standardized residuals beyond plus or minus three, including April 2016, during a strong El Nino.

## Slide 15: Remedial measures

[11:20 - 12:10] Two remedial measures were considered. First, a small residual autocorrelation at lag 3 suggested adding an autoregressive term. The extended model has a significant AR coefficient, but its SBC is higher and, out of sample, it forecasts less accurately and its 95 percent intervals cover only 83 percent of outcomes. It was therefore rejected in favor of the parsimonious model. Second, the five outlying months were treated as missing and the model re-estimated. The estimates were essentially unchanged and normality was no longer rejected, so the inference is robust.

## Slide 16: Forecast evaluation

[12:10 - 13:00] Forecast performance was evaluated with a rolling-origin design: each model was re-estimated every August from 2000 to 2022 and used to forecast the following twelve months, giving 276 out-of-sample forecasts. The seasonal ARIMA achieved an RMSE of about 0.55 parts per million, with 95 percent interval coverage of 96 percent, close to the nominal level. It is the most accurate of the four model families and essentially tied with Holt-Winters exponential smoothing. On the most recent three-year hold-out, its intervals covered 97 percent of the observations.

## Slide 17: Forecasts to 2030

[13:00 - 13:50] Finally, the model estimated on all 822 months was used to forecast to the end of 2030. Monthly CO2 is forecast to exceed 440 parts per million for the first time in April 2030, with a 95 percent prediction interval of roughly 438 to 444. Concentrations continue to rise by about 2.4 parts per million per year, and the prediction intervals widen with the horizon, reflecting the accumulated uncertainty of a doubly differenced process.

## Slide 18: Conclusions

[13:50 - 14:50] To conclude, monthly CO2 at Mauna Loa is non-stationary, with an accelerating trend and a stable seasonal cycle, and both a regular and a seasonal difference are required. The seasonal ARIMA model with one regular and one seasonal moving-average term was selected by identification and SBC. Both parameters are highly significant, the model is invertible, and its residuals are white noise. It outperforms the other model families out of sample. The results show that CO2 is rising at about 2.4 parts per million per year and is forecast to pass 440 parts per million in April 2030. The main limitation is that a univariate model cannot anticipate El Nino-driven surges; incorporating climate indices or emissions data is a natural extension.

## Slide 19: Thank you

[14:50 - 15:00] Thank you for your attention. I would be happy to take any questions.

## Backup: residual ACF of the four model families

Backup, use if asked why the other model families were rejected: the residual ACF of the best model in each family. Only the seasonal ARIMA, panel d, has no significant spikes at the seasonal lags 12, 24 and 36.

## Backup: 36-month hold-out forecasts

Backup, use if asked about the recent forecasts: forecasts made in August 2023 for the following three years. All methods under-predicted 2024, the year of record growth during the 2023-24 El Nino, which no univariate model can anticipate.

## Backup: Ljung-Box $p$-values by lag

Backup, use if asked about residual independence: Ljung-Box p-values for lags up to 36 for the final and the alternative model. The final model is above 0.05 from lag 6 onwards.

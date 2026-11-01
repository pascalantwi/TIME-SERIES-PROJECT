# =============================================================================
# STA6938 Time Series -- Final Project
# Monthly atmospheric CO2 at Mauna Loa (NOAA GML), March 1958 -- August 2026
#
# Run from the "CO2 Project" folder:   Rscript R/co2_analysis.R
# Outputs: figures/ (PNG 300 dpi + PDF), tables/ (CSV + LaTeX), results/
# =============================================================================

suppressPackageStartupMessages({
  library(forecast); library(tseries); library(urca); library(nortest)
  library(ggplot2); library(patchwork); library(scales); library(kableExtra); library(ggrepel)
})
set.seed(6938)

# ---- Helpers -----------------------------------------------------------------
f  <- function(x, d = 2) formatC(x, format = "f", digits = d, big.mark = ",")
fp <- function(p) ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3))
ym <- function(t) format(as.Date(paste(floor(t + 1e-6), round((t - floor(t + 1e-6)) * 12) + 1, 1, sep = "-")), "%b %Y")

save_fig <- function(p, name, w = 9, h = 5) {
  ggsave(file.path("figures", paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path("figures", paste0(name, ".pdf")), p, width = w, height = h)
}
save_tab <- function(df, name, align = NULL) {
  write.csv(df, file.path("tables", paste0(name, ".csv")), row.names = FALSE)
  tex <- kbl(df, format = "latex", booktabs = TRUE, align = align, linesep = "", escape = TRUE, row.names = FALSE)
  writeLines(as.character(tex), file.path("tables", paste0(name, ".tex")))
}
values <- list()
keep <- function(name, value) values[[name]] <<- value

# ---- Chart style -------------------------------------------------------------
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e7e6e2"
c_blue <- "#2a78d6"; c_orange <- "#eb6834"; c_aqua <- "#1baf7a"; c_red <- "#e34948"; c_grey <- "#8f8d86"

theme_slide <- function(base = 13) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = grid_col, linewidth = 0.35),
      axis.title = element_text(colour = ink2, size = rel(0.85)),
      axis.text = element_text(colour = ink2, size = rel(0.8)),
      plot.title = element_text(face = "bold", colour = ink, size = rel(0.95)),
      plot.title.position = "plot",
      legend.position = "top", legend.justification = "left",
      legend.title = element_blank(), legend.text = element_text(colour = ink2),
      plot.margin = margin(6, 12, 6, 6)
    )
}

acf_panel <- function(x, title, lag.max = 48, partial = FALSE, seasonal = TRUE) {
  a <- if (partial) pacf(x, lag.max = lag.max, plot = FALSE) else acf(x, lag.max = lag.max, plot = FALSE)
  df <- data.frame(lag = round(as.numeric(a$lag) * frequency(x)), r = as.numeric(a$acf))
  df <- df[df$lag > 0, ]
  df$seas <- seasonal & df$lag %% 12 == 0
  b <- 1.96 / sqrt(length(x))
  ggplot(df, aes(lag, r)) +
    geom_hline(yintercept = 0, colour = ink2, linewidth = 0.3) +
    geom_hline(yintercept = c(-b, b), colour = c_orange, linetype = "dashed", linewidth = 0.5) +
    geom_segment(aes(xend = lag, yend = 0, colour = seas), linewidth = 0.8, show.legend = FALSE) +
    geom_point(aes(colour = seas), size = 1.3, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = c_blue, `TRUE` = c_red)) +
    scale_x_continuous(breaks = seq(0, lag.max, 12)) +
    labs(title = title, x = "Lag (months)", y = if (partial) "PACF" else "ACF") + theme_slide()
}

lb <- function(x, lags = c(6, 12, 18, 24), fitdf = 0)
  t(sapply(lags, function(h) { b <- Box.test(x, h, "Ljung-Box", fitdf = fitdf); c(b$statistic, b$p.value) }))

# ---- 1. Data -----------------------------------------------------------------
raw <- read.csv("data/CO2_MaunaLoa_monthly.csv")
x <- ts(raw$CO2_ppm, start = c(1958, 3), frequency = 12)
n <- length(x); tm <- as.numeric(time(x))
dates <- as.Date(paste0(raw$Month, "-01"))
keep("n", n); keep("start", ym(tm[1])); keep("end", ym(tm[n]))
keep("first", x[1]); keep("last", x[n]); keep("mean", mean(x)); keep("sd", sd(x))
keep("min", min(x)); keep("min_date", ym(tm[which.min(x)])); keep("max", max(x)); keep("max_date", ym(tm[which.max(x)]))

# ---- 2. Classical decomposition (Note 4) ------------------------------------
dec <- decompose(x, type = "additive")             # 2x12 moving-average trend filter
seas_idx <- tapply(as.numeric(dec$seasonal), cycle(x), mean)
names(seas_idx) <- month.abb
keep("seas_peak", names(which.max(seas_idx))); keep("seas_peak_val", max(seas_idx))
keep("seas_trough", names(which.min(seas_idx))); keep("seas_trough_val", min(seas_idx))
keep("seas_amp", diff(range(seas_idx)))
detr <- x - dec$trend
yr <- floor(tm + 1e-6)
amp_year <- tapply(as.numeric(detr), yr, function(v) if (sum(!is.na(v)) == 12) diff(range(v)) else NA)
amp_year <- amp_year[!is.na(amp_year)]
keep("amp_min", min(amp_year)); keep("amp_max", max(amp_year))
keep("amp_1960s", mean(amp_year[names(amp_year) %in% 1960:1969]))
keep("amp_2010s", mean(amp_year[names(amp_year) %in% 2010:2019]))
keep("var_trend_share", var(na.omit(dec$trend)) / var(x))
rem <- na.omit(dec$random); keep("rem_sd", sd(rem))
keep("lambda", BoxCox.lambda(x))

# Growth rates from the lag-12 difference
d12 <- diff(x, lag = 12)
dec_lab <- paste0(10 * floor(floor(as.numeric(time(d12)) + 1e-6) / 10), "s")
growth_dec <- tapply(as.numeric(d12), dec_lab, mean)
level_dec <- tapply(as.numeric(x), paste0(10 * floor(yr / 10), "s"), mean)
keep("growth_1960s", growth_dec["1960s"]); keep("growth_2010s", growth_dec["2010s"]); keep("growth_2020s", growth_dec["2020s"])
keep("growth_all", mean(d12))

# ---- 3. Regression approach (Note 4): quadratic trend + seasonal dummies -------
tt <- seq_len(n); M <- factor(cycle(x), labels = month.abb[c(1:12)])
reg <- lm(as.numeric(x) ~ tt + I(tt^2) + M)
reg_res <- resid(reg)
keep("reg_r2", summary(reg)$r.squared); keep("reg_sigma", summary(reg)$sigma)
keep("reg_b1", coef(reg)["tt"]); keep("reg_b2", coef(reg)["I(tt^2)"])
keep("reg_res_acf1", acf(reg_res, 1, plot = FALSE)$acf[2])
keep("reg_lb12", Box.test(reg_res, 12, "Ljung-Box")$statistic)

# ---- 4. Stationarity: ACF and ADF (Notes 3-4) --------------------------------
d1 <- diff(x); dd <- diff(d12)
sa <- x - dec$seasonal
adf_row <- function(series, label, type, lags) {
  t(sapply(lags, function(L) {
    u <- ur.df(series, type = type, lags = L)
    c(label = label, type = type, lags = L, tau = u@teststat[1], cv1 = u@cval[1, 1], cv5 = u@cval[1, 2])
  }))
}
adf <- rbind(
  adf_row(sa, "Seasonally adjusted X_t", "drift", 0:2),
  adf_row(sa, "Seasonally adjusted X_t", "trend", 0:2),
  adf_row(x,  "X_t (13 augmentation lags)", "trend", 13),
  adf_row(d12, "Lag-12 difference", "drift", 0:2),
  adf_row(dd, "Regular + seasonal difference", "drift", 0:2)
)
adf <- data.frame(adf, stringsAsFactors = FALSE)
for (k in c("tau", "cv1", "cv5")) adf[[k]] <- as.numeric(adf[[k]])
adf$Decision <- ifelse(adf$tau < adf$cv5, "Reject H0 (stationary)", "Fail to reject (unit root)")
keep("adf_sa_drift_range", range(adf$tau[1:3])); keep("adf_sa_trend_range", range(adf$tau[4:6]))
keep("adf_raw13", adf$tau[7]); keep("adf_d12_range", range(adf$tau[8:10])); keep("adf_dd_range", range(adf$tau[11:13]))
keep("adf_cv5_drift", adf$cv5[1]); keep("adf_cv5_trend", adf$cv5[4])

acf_x <- acf(x, 36, plot = FALSE)$acf; acf_d1 <- acf(d1, 36, plot = FALSE)$acf; acf_d12 <- acf(d12, 36, plot = FALSE)$acf
keep("acf_x_1", acf_x[2]); keep("acf_x_36", acf_x[37])
keep("acf_d1_12", acf_d1[13]); keep("acf_d1_24", acf_d1[25]); keep("acf_d1_36", acf_d1[37])
keep("acf_d12_1", acf_d12[2]); keep("acf_d12_2", acf_d12[3]); keep("acf_d12_3", acf_d12[4])

# Unit-root check through the model: lag-12 difference with an AR(1) term and drift
g_d0 <- Arima(x, order = c(1, 0, 1), seasonal = c(0, 1, 1), include.drift = TRUE, method = "ML")
keep("phi_d0", coef(g_d0)["ar1"]); keep("phi_d0_se", sqrt(diag(g_d0$var.coef))["ar1"])

# ---- 5. White-noise tests on the differenced series (Note 6) -----------------
wn_dd <- lb(dd)
keep("lb_dd", wn_dd)
acf_dd <- acf(dd, 48, plot = FALSE)$acf[-1]; pacf_dd <- pacf(dd, 48, plot = FALSE)$acf
keep("acf_dd_1", acf_dd[1]); keep("acf_dd_11", acf_dd[11]); keep("acf_dd_12", acf_dd[12]); keep("acf_dd_13", acf_dd[13])
keep("pacf_dd_1", pacf_dd[1]); keep("pacf_dd_12", pacf_dd[12]); keep("pacf_dd_24", pacf_dd[24]); keep("pacf_dd_36", pacf_dd[36])
keep("bound_dd", 1.96 / sqrt(length(dd))); keep("n_dd", length(dd))

# ---- 6. Candidate ARIMA models (Notes 7-9), ML estimation ---------------------
cands <- list(c(0,1,1,0,1,1), c(1,1,0,0,1,1), c(1,1,1,0,1,1), c(0,1,2,0,1,1), c(2,1,0,0,1,1),
              c(0,1,1,1,1,1), c(0,1,1,0,1,2), c(0,1,1,1,1,0), c(1,1,1,1,1,1))
cand <- do.call(rbind, lapply(cands, function(o) {
  fit <- Arima(x, order = o[1:3], seasonal = o[4:6], method = "ML")
  k <- sum(o[c(1, 3, 4, 6)])
  data.frame(Model = sprintf("ARIMA(%d,%d,%d)(%d,%d,%d)[12]", o[1], o[2], o[3], o[4], o[5], o[6]),
             Parameters = k, LogLik = as.numeric(logLik(fit)), AIC = AIC(fit), SBC = BIC(fit),
             LB24_p = Box.test(resid(fit), 24, "Ljung-Box", fitdf = k)$p.value)
}))
cand <- cand[order(cand$SBC), ]
cand$dSBC <- cand$SBC - min(cand$SBC)
keep("cand", cand)

# ---- 7. Final model ARIMA(0,1,1)(0,1,1)[12]; alternative ARIMA(1,1,1)(0,1,1)[12] ------
# Course sign convention: MA terms enter as (1 - theta B), so theta = -(R's ma coefficient);
# AR terms enter as (1 - phi B), the same sign R reports.
course_sign <- function(fit) { b <- coef(fit); ma <- grepl("ma", names(b)); b[ma] <- -b[ma]; b }
par_label <- c(ar1 = "phi_1 (AR(1))", ma1 = "theta_1 (MA(1))", sma1 = "Theta_1 (seasonal MA(1))")
est_tab <- function(fit, model, method) {
  b <- course_sign(fit); se <- sqrt(diag(fit$var.coef)); z <- b / se
  data.frame(Model = model, Method = method, Parameter = par_label[names(b)],
             Estimate = b, SE = se, z = z, p = 2 * pnorm(-abs(z)), Lower95 = b - 1.96 * se, Upper95 = b + 1.96 * se)
}
fit_ml  <- Arima(x, order = c(0, 1, 1), seasonal = c(0, 1, 1), method = "ML")
fit_css <- Arima(x, order = c(0, 1, 1), seasonal = c(0, 1, 1), method = "CSS")
fit_alt <- Arima(x, order = c(1, 1, 1), seasonal = c(0, 1, 1), method = "ML")
est <- rbind(est_tab(fit_ml, "Final: ARIMA(0,1,1)(0,1,1)[12]", "Maximum likelihood"),
             est_tab(fit_css, "Final: ARIMA(0,1,1)(0,1,1)[12]", "Conditional least squares"),
             est_tab(fit_alt, "Alternative: ARIMA(1,1,1)(0,1,1)[12]", "Maximum likelihood"))
keep("est", est); keep("sigma2_ml", fit_ml$sigma2); keep("sigma2_css", fit_css$sigma2); keep("sigma2_alt", fit_alt$sigma2)
keep("loglik", as.numeric(logLik(fit_ml))); keep("aic", AIC(fit_ml)); keep("sbc", BIC(fit_ml))
keep("aic_alt", AIC(fit_alt)); keep("sbc_alt", BIC(fit_alt))
th <- course_sign(fit_ml); th_alt <- course_sign(fit_alt)
keep("root_ma", 1 / th["ma1"]); keep("root_sma", (1 / th["sma1"])^(1 / 12))
keep("root_ar_alt", 1 / th_alt["ar1"]); keep("root_ma_alt", 1 / th_alt["ma1"])
# Autocorrelations of the differenced series implied by the fitted model
keep("rho1_model", -th[1] / (1 + th[1]^2)); keep("rho12_model", -th[2] / (1 + th[2]^2))
keep("rho13_model", th[1] * th[2] / ((1 + th[1]^2) * (1 + th[2]^2)))
# pi weights of the invertible form (1 - theta B)^-1 (1 - Theta B^12)^-1 for the differenced series
keep("pi_1_3", th[1]^(1:3)); keep("pi_12", th[2])

# ---- 8. Diagnostics (Notes 7-9) -----------------------------------------------
e <- as.numeric(resid(fit_ml)); z <- e / sqrt(fit_ml$sigma2)
e_alt <- as.numeric(resid(fit_alt))
keep("acf_final_1_6", round(acf(e, 6, plot = FALSE)$acf[2:7], 3))
keep("acf_alt_1_6", round(acf(e_alt, 6, plot = FALSE)$acf[2:7], 3))
keep("lb_short_final", round(sapply(3:5, function(h) Box.test(e, h, "Ljung-Box", fitdf = 2)$p.value), 4))
keep("n_acf_out", sum(abs(acf(e, 36, plot = FALSE)$acf[-1]) > 1.96 / sqrt(n)))
diag_tab <- function(e, label, k = 2) {
  b <- lb(e, fitdf = k); s <- lb(e^2, c(12, 24))
  data.frame(Residuals = label,
             Test = c(paste0("Ljung-Box Q(", c(6, 12, 18, 24), ")"), "Ljung-Box Q(12), squared", "Ljung-Box Q(24), squared",
                      "Shapiro-Wilk", "Anderson-Darling", "Jarque-Bera"),
             Statistic = c(b[, 1], s[, 1], shapiro.test(e)$statistic, ad.test(e)$statistic, jarque.bera.test(e)$statistic),
             p = c(b[, 2], s[, 2], shapiro.test(e)$p.value, ad.test(e)$p.value, jarque.bera.test(e)$p.value))
}
out_idx <- which(abs(z) > 3)
diag_all <- diag_tab(e, "All residuals")
diag_alt <- diag_tab(e_alt, "Alternative model", k = 3)
keep("diag_alt", diag_alt)
keep("diag_all", diag_all)
keep("skew", mean(z^3)); keep("exkurt", mean(z^4) - 3)
keep("n_out", length(out_idx)); keep("out_dates", ym(tm[out_idx])); keep("out_z", z[out_idx])
keep("var_first_half", var(e[1:floor(n / 2)])); keep("var_second_half", var(e[(floor(n / 2) + 1):n]))

# Remedy: treat the outlying months as missing and re-estimate (Kalman filter handles NA)
x_na <- x; x_na[out_idx] <- NA
fit_rob <- Arima(x_na, order = c(0, 1, 1), seasonal = c(0, 1, 1), method = "ML")
e_rob <- as.numeric(resid(fit_rob)); e_rob <- e_rob[!is.na(e_rob)]
diag_rob <- diag_tab(e_rob, "Outlying months set to missing")
keep("diag_rob", diag_rob); keep("est_rob", course_sign(fit_rob)); keep("est_rob_se", sqrt(diag(fit_rob$var.coef)))
keep("sigma2_rob", fit_rob$sigma2)
lb_curve <- function(e, k, model) data.frame(lag = (k + 1):36, model = model,
  p = sapply((k + 1):36, function(h) Box.test(e, h, "Ljung-Box", fitdf = k)$p.value))
lb_lags <- rbind(lb_curve(e, 2, "Final: ARIMA(0,1,1)(0,1,1)[12]"), lb_curve(e_alt, 3, "Alternative: ARIMA(1,1,1)(0,1,1)[12]"))
fin_lb <- lb_lags[grepl("Final", lb_lags$model), ]
keep("lb_min_p_final", min(fin_lb$p)); keep("lb_min_lag_final", fin_lb$lag[which.min(fin_lb$p)])

# ---- 9. Holt-Winters (Note 5) ---------------------------------------------------
hw <- HoltWinters(x, seasonal = "additive")
keep("hw_alpha", hw$alpha); keep("hw_beta", hw$beta); keep("hw_gamma", hw$gamma)
keep("hw_level", hw$coefficients["a"]); keep("hw_slope", hw$coefficients["b"])
hw_res <- residuals(hw); keep("hw_rmse_in", sqrt(mean(hw_res^2))); keep("arima_rmse_in", sqrt(mean(e[-(1:13)]^2)))

# ---- 10. Out-of-sample comparison: last 36 months -------------------------------
H <- 36
tr <- window(x, end = tm[n - H]); te <- window(x, start = tm[n - H + 1])
fc_ar <- forecast(Arima(tr, order = c(0, 1, 1), seasonal = c(0, 1, 1), method = "ML"), h = H, level = 95)
fc_ai <- forecast(Arima(tr, order = c(1, 1, 1), seasonal = c(0, 1, 1), method = "ML"), h = H, level = 95)
fc_hw <- forecast(HoltWinters(tr, seasonal = "additive"), h = H, level = 95)
ntr <- length(tr); ttr <- seq_len(ntr); Mtr <- factor(cycle(tr), levels = 1:12)
reg_tr <- lm(as.numeric(tr) ~ ttr + I(ttr^2) + Mtr)
nd <- data.frame(ttr = (ntr + 1):(ntr + H), Mtr = factor(cycle(te), levels = 1:12))
pr <- predict(reg_tr, nd, interval = "prediction", level = 0.95)
fc_sn <- snaive(tr, h = H, level = 95)
acc <- function(fc, lo, hi) c(RMSE = sqrt(mean((te - fc)^2)), MAE = mean(abs(te - fc)),
                              MAPE = 100 * mean(abs((te - fc) / te)), Coverage = 100 * mean(te >= lo & te <= hi))
acc_tab <- rbind(
  "Final: ARIMA(0,1,1)(0,1,1)[12]" = acc(fc_ar$mean, fc_ar$lower, fc_ar$upper),
  "Alternative: ARIMA(1,1,1)(0,1,1)[12]" = acc(fc_ai$mean, fc_ai$lower, fc_ai$upper),
  "Holt-Winters (additive)" = acc(fc_hw$mean, fc_hw$lower, fc_hw$upper),
  "Regression: quadratic trend + seasonal dummies" = acc(pr[, 1], pr[, 2], pr[, 3]),
  "Seasonal naive" = acc(fc_sn$mean, fc_sn$lower, fc_sn$upper)
)
keep("acc", acc_tab); keep("test_start", ym(tm[n - H + 1])); keep("test_end", ym(tm[n]))

# ---- 11. Forecasts to December 2030 ------------------------------------------------
h_end <- (2030 - 2026) * 12 + (12 - 8)
fc <- forecast(fit_ml, h = h_end, level = c(80, 95))
tf <- as.numeric(time(fc$mean))
fc_df <- data.frame(Month = ym(tf), Forecast = as.numeric(fc$mean),
                    Lower80 = as.numeric(fc$lower[, 1]), Upper80 = as.numeric(fc$upper[, 1]),
                    Lower95 = as.numeric(fc$lower[, 2]), Upper95 = as.numeric(fc$upper[, 2]))
fc_hw_full <- forecast(hw, h = h_end, level = 95)
fc_df$HoltWinters <- as.numeric(fc_hw_full$mean)
keep("fc_first", fc_df[1, ]); keep("fc_may2027", fc_df[fc_df$Month == "May 2027", ])
keep("fc_may2030", fc_df[fc_df$Month == "May 2030", ]); keep("fc_dec2030", fc_df[fc_df$Month == "Dec 2030", ])
keep("fc_440", fc_df$Month[which(fc_df$Forecast > 440)[1]])
keep("fc_year_means", tapply(fc_df$Forecast, substr(fc_df$Month, 5, 8), mean))
keep("hw_dec2030", fc_df$HoltWinters[nrow(fc_df)])
keep("pi_width_1", fc_df$Upper95[1] - fc_df$Lower95[1]); keep("pi_width_end", fc_df$Upper95[nrow(fc_df)] - fc_df$Lower95[nrow(fc_df)])

# =============================================================================
# TABLES (CSV + LaTeX booktabs, for the slides)
# =============================================================================
save_tab(data.frame(
  Statistic = c("Observations (months)", "First month", "Last month", "Mean (ppm)", "Standard deviation (ppm)",
                "Minimum (ppm)", "Maximum (ppm)", "Average annual increase (ppm/yr)"),
  Value = c(n, ym(tm[1]), ym(tm[n]), f(mean(x)), f(sd(x)),
            paste0(f(min(x)), " (", ym(tm[which.min(x)]), ")"), paste0(f(max(x)), " (", ym(tm[which.max(x)]), ")"),
            f(mean(d12)))), "tab01_descriptive_statistics", "lr")

amp_dec <- tapply(amp_year, paste0(10 * floor(as.numeric(names(amp_year)) / 10), "s"), mean)
decs <- names(level_dec)
save_tab(data.frame(Decade = decs, `Mean CO2 (ppm)` = f(level_dec[decs]),
                    `Annual increase (ppm/yr)` = f(growth_dec[decs]),
                    `Seasonal amplitude (ppm)` = ifelse(is.na(amp_dec[decs]), "--", f(amp_dec[decs])), check.names = FALSE),
         "tab02_growth_by_decade", "lrrr")

save_tab(data.frame(Month = month.abb, `Seasonal index (ppm)` = f(seas_idx, 3), check.names = FALSE),
         "tab03_seasonal_indices", "lr")

save_tab(data.frame(Series = adf$label, Type = ifelse(adf$type == "drift", "Single mean", "Trend"), Lags = adf$lags,
                    Tau = f(adf$tau), `5% critical value` = f(adf$cv5), `1% critical value` = f(adf$cv1),
                    Decision = adf$Decision, check.names = FALSE), "tab04_adf_tests", "llrrrrl")

rc <- coef(summary(reg))
save_tab(data.frame(Term = c("Intercept", "t", "t^2", paste("Month:", month.abb[2:12])),
                    Estimate = f(rc[, 1], 4), SE = f(rc[, 2], 4), t = f(rc[, 3]), p = fp(rc[, 4])),
         "tab05_regression_trend_seasonal", "lrrrr")

wn <- data.frame(Lag = c(6, 12, 18, 24), Q = f(wn_dd[, 1]), p = fp(wn_dd[, 2]))
save_tab(wn, "tab06_white_noise_differenced", "rrr")

save_tab(data.frame(Model = cand$Model, Parameters = cand$Parameters, `Log-likelihood` = f(cand$LogLik),
                    AIC = f(cand$AIC), SBC = f(cand$SBC), `Delta SBC` = f(cand$dSBC),
                    `Ljung-Box Q(24) p` = fp(cand$LB24_p), check.names = FALSE), "tab07_candidate_models", "lrrrrrr")

save_tab(data.frame(Model = est$Model, Method = est$Method, Parameter = est$Parameter, Estimate = f(est$Estimate, 4), SE = f(est$SE, 4),
                    z = f(est$z), p = fp(est$p), `95% CI` = paste0("[", f(est$Lower95, 3), ", ", f(est$Upper95, 3), "]"),
                    check.names = FALSE), "tab08_parameter_estimates", "lllrrrrc")

save_tab(data.frame(Test = diag_all$Test,
                    `Final: statistic` = f(diag_all$Statistic, 3), `Final: p` = fp(diag_all$p),
                    `Final, outliers missing: statistic` = f(diag_rob$Statistic, 3), `Final, outliers missing: p` = fp(diag_rob$p),
                    `Alternative: statistic` = f(diag_alt$Statistic, 3), `Alternative: p` = fp(diag_alt$p),
                    check.names = FALSE), "tab09_residual_diagnostics", "lrrrrrr")

save_tab(data.frame(Month = ym(tm[out_idx]), `Observed (ppm)` = f(x[out_idx]),
                    `Standardized residual` = f(z[out_idx])), "tab10_outlying_months", "lrr")

rb <- course_sign(fit_rob)
save_tab(data.frame(Data = c("All months", "Outlying months set to missing"),
                    theta_1 = f(c(th["ma1"], rb["ma1"]), 4), Theta_1 = f(c(th["sma1"], rb["sma1"]), 4),
                    sigma2 = f(c(fit_ml$sigma2, fit_rob$sigma2), 4)),
         "tab11_sensitivity_estimates", "lrrr")

save_tab(data.frame(Parameter = c("alpha (level)", "beta (trend)", "gamma (seasonal)", "Final level (ppm)", "Final slope (ppm/month)"),
                    Value = f(c(hw$alpha, hw$beta, hw$gamma, hw$coefficients["a"], hw$coefficients["b"]), 4)),
         "tab12_holt_winters_parameters", "lr")

save_tab(data.frame(Method = rownames(acc_tab), RMSE = f(acc_tab[, "RMSE"], 3), MAE = f(acc_tab[, "MAE"], 3),
                    `MAPE (%)` = f(acc_tab[, "MAPE"], 3), `95% PI coverage (%)` = f(acc_tab[, "Coverage"], 1),
                    check.names = FALSE), "tab13_holdout_accuracy", "lrrrr")

write.csv(transform(fc_df, Forecast = round(Forecast, 2), Lower80 = round(Lower80, 2), Upper80 = round(Upper80, 2),
                    Lower95 = round(Lower95, 2), Upper95 = round(Upper95, 2), HoltWinters = round(HoltWinters, 2)),
          "tables/tab14_forecasts_full.csv", row.names = FALSE)
sel <- fc_df$Month %in% c("Sep 2026", "Dec 2026", paste(rep(c("May", "Dec"), 4), rep(2027:2030, each = 2)))
save_tab(data.frame(Month = fc_df$Month[sel], Forecast = f(fc_df$Forecast[sel]),
                    `80% PI` = paste0("[", f(fc_df$Lower80[sel]), ", ", f(fc_df$Upper80[sel]), "]"),
                    `95% PI` = paste0("[", f(fc_df$Lower95[sel]), ", ", f(fc_df$Upper95[sel]), "]"),
                    `Holt-Winters` = f(fc_df$HoltWinters[sel]), check.names = FALSE), "tab14_forecasts_selected", "lrccr")

# =============================================================================
# FIGURES (PNG 300 dpi + PDF, 16:9-friendly for slides)
# =============================================================================
df <- data.frame(Date = dates, co2 = as.numeric(x), trend = as.numeric(dec$trend))
p <- ggplot(df, aes(Date)) +
  geom_line(aes(y = co2, colour = "Monthly mean"), linewidth = 0.45) +
  geom_line(aes(y = trend, colour = "Trend (2x12 moving average)"), linewidth = 1) +
  scale_colour_manual(values = c("Monthly mean" = c_blue, "Trend (2x12 moving average)" = c_orange)) +
  annotate("label", x = dates[n] - 3000, y = x[n] + 4, label = paste0(f(x[n]), " ppm (", ym(tm[n]), ")"),
           size = 3.6, colour = ink, fill = "white", label.size = 0) +
  annotate("label", x = dates[1] + 3000, y = x[1] - 3, label = paste0(f(x[1]), " ppm (", ym(tm[1]), ")"),
           size = 3.6, colour = ink, fill = "white", label.size = 0) +
  scale_y_continuous(expand = expansion(mult = 0.06)) +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
  labs(x = NULL, y = expression(CO[2] ~ "(ppm)")) + theme_slide()
save_fig(p, "fig01_co2_series")

dd_df <- rbind(
  data.frame(Date = dates, value = as.numeric(x), panel = "Observed"),
  data.frame(Date = dates, value = as.numeric(dec$trend), panel = "Trend"),
  data.frame(Date = dates, value = as.numeric(dec$seasonal), panel = "Seasonal"),
  data.frame(Date = dates, value = as.numeric(dec$random), panel = "Remainder"))
dd_df$panel <- factor(dd_df$panel, levels = c("Observed", "Trend", "Seasonal", "Remainder"))
p <- ggplot(dd_df, aes(Date, value)) + geom_line(colour = c_blue, linewidth = 0.4, na.rm = TRUE) +
  facet_wrap(~panel, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
  labs(x = NULL, y = NULL) + theme_slide() +
  theme(strip.placement = "outside", strip.text = element_text(face = "bold", colour = ink2))
save_fig(p, "fig02_classical_decomposition", h = 6.5)

si <- data.frame(Month = factor(month.abb, levels = month.abb), idx = as.numeric(seas_idx))
p <- ggplot(si, aes(Month, idx, fill = idx > 0)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = 0.4) +
  geom_text(aes(label = f(idx), vjust = ifelse(idx > 0, -0.4, 1.3)), size = 3.6, colour = ink) +
  scale_fill_manual(values = c(`TRUE` = c_blue, `FALSE` = c_orange)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(x = NULL, y = "Seasonal index (ppm)") + theme_slide()
save_fig(p, "fig03_seasonal_indices")

ad <- data.frame(Year = as.numeric(names(amp_year)), amp = as.numeric(amp_year))
p <- ggplot(ad, aes(Year, amp)) +
  geom_line(colour = c_blue, linewidth = 0.6) + geom_point(colour = c_blue, size = 1.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = c_orange, linewidth = 0.8) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "Within-year range of detrended CO2 (ppm)") + theme_slide()
save_fig(p, "fig04_seasonal_amplitude_by_year")

p <- (acf_panel(x, "(a) Original series", 48) + coord_cartesian(ylim = c(-0.2, 1))) /
     (acf_panel(d1, "(b) First difference", 48) + coord_cartesian(ylim = c(-1, 1))) /
     (acf_panel(d12, "(c) Lag-12 (seasonal) difference", 48) + coord_cartesian(ylim = c(-0.2, 1)))
save_fig(p, "fig05_acf_original_and_differences", h = 8)

gdf <- data.frame(Date = dates[-(1:12)], g = as.numeric(d12), dec = dec_lab)
gm <- aggregate(g ~ dec, gdf, mean); gm$start <- as.Date(paste0(substr(gm$dec, 1, 4), "-01-01"))
gm$start[1] <- min(gdf$Date); gm$end <- c(gm$start[-1], max(gdf$Date))
p <- ggplot(gdf, aes(Date, g)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = 0.3) +
  geom_line(aes(colour = "Lag-12 difference (ppm/yr)"), linewidth = 0.4) +
  geom_segment(data = gm, aes(x = start, xend = end, y = g, yend = g, colour = "Decade average"), linewidth = 1.2) +
  geom_label(data = gm[gm$dec != "1950s", ], aes(x = start + (end - start) / 2, y = g + 0.35, label = f(g)), size = 3.4,
             colour = ink, fill = "white", label.size = 0) +
  scale_colour_manual(values = c("Lag-12 difference (ppm/yr)" = c_blue, "Decade average" = c_orange)) +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
  labs(x = NULL, y = "Increase over 12 months (ppm)") + theme_slide()
save_fig(p, "fig06_annual_growth_lag12_difference")

p <- ggplot(data.frame(Date = dates[-(1:13)], v = as.numeric(dd)), aes(Date, v)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = 0.3) +
  geom_line(colour = c_blue, linewidth = 0.4) +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
  labs(x = NULL, y = expression(nabla * nabla[12] * X[t] ~ "(ppm)")) + theme_slide()
save_fig(p, "fig07_regular_and_seasonal_difference")

p <- (acf_panel(dd, "(a) Sample ACF") + coord_cartesian(ylim = c(-0.55, 0.3))) /
     (acf_panel(dd, "(b) Sample PACF", partial = TRUE) + coord_cartesian(ylim = c(-0.55, 0.3)))
save_fig(p, "fig08_acf_pacf_differenced", h = 6.5)

rdf <- data.frame(Date = dates, r = reg_res)
p <- (ggplot(rdf, aes(Date, r)) + geom_hline(yintercept = 0, colour = ink2, linewidth = 0.3) +
        geom_line(colour = c_blue, linewidth = 0.4) +
        scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
        labs(title = "(a) Residuals of quadratic trend + seasonal dummies", x = NULL, y = "ppm") + theme_slide()) |
     (acf_panel(ts(reg_res, frequency = 12), "(b) ACF of regression residuals", 36, seasonal = FALSE) +
        coord_cartesian(ylim = c(-0.2, 1)))
save_fig(p, "fig09_regression_residuals")


zdf <- data.frame(Date = dates, z = z)
odf <- zdf[out_idx, ]; odf$lab <- format(odf$Date, "%b %Y")
p1 <- ggplot(zdf, aes(Date, z)) + geom_hline(yintercept = 0, colour = ink2, linewidth = 0.3) +
  geom_hline(yintercept = c(-3, 3), colour = c_orange, linetype = "dashed", linewidth = 0.5) +
  geom_line(colour = c_blue, linewidth = 0.35) +
  geom_point(data = odf, colour = c_red, size = 2) +
  geom_text_repel(data = odf, aes(label = lab), size = 3, colour = ink, min.segment.length = 0, seed = 1,
                  nudge_y = ifelse(odf$z > 0, 0.6, -0.6), segment.colour = ink2) +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(title = "(a) Standardized residuals", x = NULL, y = "Residual / SD") + theme_slide()
p2 <- acf_panel(ts(e, frequency = 12), "(b) ACF of residuals", 36) + coord_cartesian(ylim = c(-0.2, 0.2))
p3 <- ggplot(zdf, aes(z)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.25, fill = c_blue, colour = "white", linewidth = 0.2) +
  stat_function(fun = dnorm, colour = c_orange, linewidth = 0.8) +
  labs(title = "(c) Histogram and N(0,1) density", x = "Standardized residual", y = "Density") + theme_slide()
p4 <- ggplot(zdf, aes(sample = z)) + stat_qq(colour = c_blue, size = 0.9) + stat_qq_line(colour = c_orange, linewidth = 0.8) +
  labs(title = "(d) Normal Q-Q plot", x = "Normal quantiles", y = "Sample quantiles") + theme_slide()
save_fig((p1 | p2) / (p3 | p4), "fig10_final_model_residual_diagnostics", h = 6.5)

p <- ggplot(lb_lags, aes(lag, p, colour = model, shape = model)) +
  geom_hline(yintercept = 0.05, colour = c_orange, linetype = "dashed", linewidth = 0.6) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(c_blue, c_grey)) + scale_shape_manual(values = c(16, 1)) +
  annotate("text", x = 36, y = 0.05, label = "p = 0.05", hjust = 1, vjust = -0.6, size = 3.6, colour = ink2) +
  scale_x_continuous(breaks = seq(0, 36, 6)) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Lag h (months)", y = "Ljung-Box p-value") + theme_slide()
save_fig(p, "fig11_ljung_box_pvalues")

hdf <- data.frame(Date = dates[(n - H + 1):n], Actual = as.numeric(te),
                  ARIMA = as.numeric(fc_ar$mean), lo = as.numeric(fc_ar$lower), hi = as.numeric(fc_ar$upper),
                  HW = as.numeric(fc_hw$mean), Reg = pr[, 1])
hist_tr <- data.frame(Date = dates[(n - H - 35):(n - H)], y = as.numeric(x[(n - H - 35):(n - H)]))
p <- ggplot(hdf, aes(Date)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = "ARIMA 95% prediction interval"), alpha = 0.2) +
  geom_line(data = hist_tr, aes(Date, y), colour = c_grey, linewidth = 0.6) +
  geom_line(aes(y = Actual, colour = "Actual"), linewidth = 0.9) +
  geom_line(aes(y = ARIMA, colour = "ARIMA(0,1,1)(0,1,1)[12]"), linewidth = 0.8) +
  geom_line(aes(y = HW, colour = "Holt-Winters"), linewidth = 0.8, linetype = "longdash") +
  geom_line(aes(y = Reg, colour = "Regression"), linewidth = 0.9, linetype = "dotted") +
  geom_vline(xintercept = dates[n - H + 1], colour = ink2, linewidth = 0.3) +
  scale_colour_manual(values = c("Actual" = ink, "ARIMA(0,1,1)(0,1,1)[12]" = c_blue, "Holt-Winters" = c_aqua, "Regression" = c_orange),
                      breaks = c("Actual", "ARIMA(0,1,1)(0,1,1)[12]", "Holt-Winters", "Regression")) +
  scale_fill_manual(values = c("ARIMA 95% prediction interval" = c_blue)) +
  guides(colour = guide_legend(order = 1, override.aes = list(linetype = c("solid", "solid", "longdash", "dotted"))),
         fill = guide_legend(order = 2)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = NULL, y = expression(CO[2] ~ "(ppm)")) + theme_slide() +
  theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"), legend.margin = margin(0, 0, 0, 0))
save_fig(p, "fig12_holdout_forecast_comparison")

fdf <- data.frame(Date = seq(dates[n], by = "month", length.out = h_end + 1)[-1], fc_df)
br <- rbind(data.frame(Date = dates[n], Forecast = x[n], Lower80 = x[n], Upper80 = x[n], Lower95 = x[n], Upper95 = x[n]),
            fdf[, c("Date", "Forecast", "Lower80", "Upper80", "Lower95", "Upper95")])
hist_df <- data.frame(Date = dates[dates >= as.Date("2018-01-01")], y = as.numeric(x[dates >= as.Date("2018-01-01")]))
lastlab <- fdf[nrow(fdf), ]
p <- ggplot() +
  geom_ribbon(data = br, aes(Date, ymin = Lower95, ymax = Upper95, fill = "95% prediction interval"), alpha = 0.18) +
  geom_ribbon(data = br, aes(Date, ymin = Lower80, ymax = Upper80, fill = "80% prediction interval"), alpha = 0.32) +
  geom_line(data = hist_df, aes(Date, y, colour = "Observed"), linewidth = 0.7) +
  geom_line(data = br, aes(Date, Forecast, colour = "Forecast"), linewidth = 0.8) +
  geom_hline(yintercept = 440, colour = ink2, linewidth = 0.3, linetype = "dotted") +
  annotate("text", x = as.Date("2018-03-01"), y = 440, label = "440 ppm", vjust = -0.5, hjust = 0, size = 3.4, colour = ink2) +
  annotate("text", x = lastlab$Date, y = lastlab$Forecast - 1.6, label = paste0(f(lastlab$Forecast), "\nDec 2030"),
           hjust = -0.1, size = 3.3, colour = ink) +
  scale_colour_manual(values = c("Observed" = ink, "Forecast" = c_blue), breaks = c("Observed", "Forecast")) +
  scale_fill_manual(values = c("80% prediction interval" = c_blue, "95% prediction interval" = c_blue),
                    breaks = c("80% prediction interval", "95% prediction interval")) +
  guides(fill = guide_legend(override.aes = list(alpha = c(0.45, 0.2)))) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = expansion(mult = c(0.01, 0.1))) +
  labs(x = NULL, y = expression(CO[2] ~ "(ppm)")) + theme_slide()
save_fig(p, "fig13_forecast_to_2030")

# ---- Key values used in the report ----------------------------------------------
sink("results/key_values.txt")
for (k in names(values)) { cat("\n==", k, "==\n"); print(values[[k]]) }
sink()
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
cat("Done:", length(list.files("figures")), "figure files,", length(list.files("tables")), "table files\n")

"""Typeset the analysis tables (../tables/*.tex, written by R/co2_analysis.R) for the
report appendix: Greek symbols, true minus signs, ARIMA notation and compact headers.

Run from the report folder:  python3 make_appendix_tables.py
"""
import pathlib
import re

SRC = pathlib.Path("../tables")
OUT = pathlib.Path("tables")
OUT.mkdir(exist_ok=True)

WORDS = [
    ("theta\\_1", r"$\theta_1$"), ("Theta\\_1", r"$\Theta_1$"), ("phi\\_1", r"$\phi_1$"),
    ("W\\_t", r"$W_t$"), ("sigma2", r"$\hat\sigma^2_Z$"), ("Ljung-Box", "Ljung--Box"),
    ("Shapiro-Wilk", "Shapiro--Wilk"), ("Anderson-Darling", "Anderson--Darling"),
    ("Jarque-Bera", "Jarque--Bera"), ("Holt-Winters", "Holt--Winters"),
    ("alpha (level)", r"$\alpha$ (level)"), ("beta (trend)", r"$\beta$ (trend)"),
    ("gamma (seasonal)", r"$\gamma$ (seasonal)"), ("t\\textasciicircum{}2", r"$t^2$"),
    ("< 0.001", r"$<$\,0.001"), ("Delta SBC", r"$\Delta$SBC"),
    ("X\\_t", r"$X_t$"), ("& Tau &", r"& $\tau$ &"), ("H0", r"$H_0$"), ("Mean CO2", r"Mean CO$_2$"),
    ("R-squared", r"$R^2$"), ("Observed..ppm.", "Observed (ppm)"), ("Standardized.residual", "Standardized residual"),
    ("Maximum likelihood", "ML"), ("Conditional least squares", "CLS"), ("Lag & Q & p", r"Lag & $Q$ & $p$"),
]

HEADERS = {
    "tab09_residual_diagnostics": (
        r"Test & Final: statistic & Final: p & Final, outliers missing: statistic & Final, outliers missing: p & Alternative: statistic & Alternative: p\\",
        r" & \multicolumn{2}{c}{Final model} & \multicolumn{2}{c}{Outlying months missing} & \multicolumn{2}{c}{Alternative model}\\"
        "\n" r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}"
        "\n" r"Test & Statistic & $p$ & Statistic & $p$ & Statistic & $p$\\"),
    "tab16_model_families": (
        r"Family & Model & Parameters & Residual ACF lag 12 & Ljung-Box Q(12) p & Ljung-Box Q(24) p & Hold-out RMSE & Hold-out 95\% coverage (\%)\\",
        r" & & & & \multicolumn{2}{c}{Ljung--Box $p$} & \multicolumn{2}{c}{Hold-out}\\"
        "\n" r"\cmidrule(lr){5-6}\cmidrule(lr){7-8}"
        "\n" r"Family & Best model & $k$ & ACF lag 12 & $Q(12)$ & $Q(24)$ & RMSE & Coverage (\%)\\"),
}


def arima(m):
    p, s = m.group(1), m.group(2)
    return rf"ARIMA$({p})\times({s})_{{12}}$"


for src in sorted(SRC.glob("*.tex")):
    tex = src.read_text(encoding="utf-8").strip()
    if src.stem in HEADERS:
        old, new = HEADERS[src.stem]
        assert old in tex, src.stem
        tex = tex.replace(old, new)
    for a, b in WORDS:
        tex = tex.replace(a, b)
    tex = re.sub(r"ARIMA\((\d,\d,\d)\)×\((\d,\d,\d)\)\[12\]", arima, tex)
    tex = re.sub(r"(?<![\w.$])-(?=\d)", "$-$", tex)            # true minus signs
    tex = re.sub(r"(?<=[ (])Q\((\d+)\)", r"$Q(\1)$", tex)     # Q(12) -> $Q(12)$
    tex = re.sub(r"& t & p", r"& $t$ & p", tex)
    tex = re.sub(r"& z & p", r"& $z$ & p", tex)
    tex = re.sub(r"& p\\\\", r"& $p$\\\\", tex)
    tex = tex.replace(r"& $z$ & p &", r"& $z$ & $p$ &")
    tex = re.sub(r"(?m)^t & ", r"$t$ & ", tex)
    tex = tex.replace(r"\begin{tabular}[t]", r"\begin{tabular}")
    (OUT / src.name).write_text(tex + "\n", encoding="utf-8")
    print("wrote", OUT / src.name)

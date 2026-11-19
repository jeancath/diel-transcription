"""Cosinor rhythm detection.

For each gene the log2 expression across a time course is fitted with

    y = sum_s b_s * I(series = s) + A*cos(2*pi*t/T) + B*sin(2*pi*t/T)

Per-series intercepts absorb between-replicate / between-platform offsets, so a
gene is called rhythmic on the shape of its time course, not on a batch shift.
Significance is an F-test of the two harmonic terms against the series-only
model; p values are corrected across genes with Benjamini-Hochberg.

Amplitude is reported in log2 units as the peak-to-trough swing (2*sqrt(A^2+B^2))
and phase as the time of the fitted peak, folded into [0, T).
"""
import numpy as np
import pandas as pd
from scipy import stats


def _design(t, series, period):
    t = np.asarray(t, float)
    levels = pd.Index(pd.unique(pd.Series(series)))
    S = np.column_stack([(np.asarray(series) == lv).astype(float) for lv in levels])
    H = np.column_stack([np.cos(2 * np.pi * t / period), np.sin(2 * np.pi * t / period)])
    return S, H


def cosinor(expr, meta, period=24.0):
    """expr: genes x samples (log2). meta: index=samples, columns time, series."""
    meta = meta.loc[expr.columns]
    S, H = _design(meta["time"].values, meta["series"].values, period)
    X = np.hstack([S, H])
    Y = expr.to_numpy(float)
    ok = np.isfinite(Y).all(axis=1)
    Yf = Y[ok]
    n, p_full, p_red = X.shape[0], X.shape[1], S.shape[1]
    df_n, df_d = p_full - p_red, n - p_full
    if df_d <= 0:
        raise ValueError("not enough samples for this design")

    def rss(design):
        coef, *_ = np.linalg.lstsq(design, Yf.T, rcond=None)
        resid = Yf.T - design @ coef
        return (resid ** 2).sum(axis=0), coef

    rss_f, coef_f = rss(X)
    rss_r, _ = rss(S)
    with np.errstate(divide="ignore", invalid="ignore"):
        F = ((rss_r - rss_f) / df_n) / (rss_f / df_d)
    pval = stats.f.sf(F, df_n, df_d)
    A, B = coef_f[-2], coef_f[-1]
    amp = 2 * np.sqrt(A ** 2 + B ** 2)
    phase = (np.arctan2(B, A) * period / (2 * np.pi)) % period   # time of fitted peak
    out = pd.DataFrame(index=pd.Index(expr.index[ok], name="gene"))
    out["amp_log2"] = amp
    out["phase_h"] = phase
    out["p"] = np.where(np.isfinite(pval), pval, 1.0)
    out["q"] = bh(out["p"].values)
    out["mean_log2"] = Yf.mean(axis=1)
    out["period_used"] = period
    return out


def bh(p):
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n)
    ranked = p[order] * n / (np.arange(n) + 1)
    q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return np.clip(q, 0, 1)


def best_period(expr, meta, periods=np.arange(20.0, 28.5, 0.5)):
    """Period that minimises total residual sum of squares across all genes.

    Used only as a sensitivity check on the fixed 24 h assumption; for
    free-running (LL) series the true period need not be exactly 24 h.
    """
    meta = meta.loc[expr.columns]
    Y = expr.to_numpy(float)
    Y = Y[np.isfinite(Y).all(axis=1)]
    best, best_rss = None, np.inf
    for T in periods:
        S, H = _design(meta["time"].values, meta["series"].values, T)
        X = np.hstack([S, H])
        coef, *_ = np.linalg.lstsq(X, Y.T, rcond=None)
        r = ((Y.T - X @ coef) ** 2).sum()
        if r < best_rss:
            best, best_rss = T, r
    return float(best)

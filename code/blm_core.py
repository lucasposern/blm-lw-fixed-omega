"""
blm_core.py
===========
Shared helpers for the tau/Omega independence report.

Implements:
  - sample covariance (annualised from monthly excess returns)
  - Ledoit-Wolf shrinkage in the report's analytic form
        Sigma_LW = (1 - alpha) * Sigma_s + alpha * mu_s * I,   mu_s = tr(Sigma_s)/n
  - effective view absorption  c_eff = tau*v / (tau*v + omega)
  - single-view Black-Litterman posterior  mu_BL
  - Idzorek omega calibration  omega = (1-c)/c * tau * p' Sigma p

All formulas match the equations in the superseded old draft
(tau_omega_independence_BLM.tex, removed from the tree; see git history).
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from sklearn.covariance import LedoitWolf

MONTHS_PER_YEAR = 12


def download_monthly_returns(tickers, start, end):
    """Monthly simple returns from Yahoo Finance adjusted close.

    Returns a DataFrame indexed by month-end, one column per ticker.
    A near-constant monthly risk-free rate has a negligible effect on the
    covariance, so we work with total returns here; the report computes the
    covariance on excess returns, which shifts the matrices only marginally.
    """
    import yfinance as yf

    px = yf.download(
        list(tickers), start=start, end=end,
        interval="1mo", auto_adjust=True, progress=False,
    )["Close"]
    px = px[list(tickers)].dropna(how="all")
    returns = px.pct_change().dropna(how="any")
    return returns


def sample_cov(returns: pd.DataFrame, annualize: bool = True) -> np.ndarray:
    """Sample covariance (ddof = 1), annualised by 12 for monthly data."""
    cov = np.cov(returns.values, rowvar=False, ddof=1)
    if annualize:
        cov = cov * MONTHS_PER_YEAR
    return cov


def ledoit_wolf(returns: pd.DataFrame, annualize: bool = True):
    """Analytic Ledoit-Wolf shrinkage in the report's form.

    Returns (Sigma_LW, alpha, mu_s). The intensity alpha comes from the
    sklearn analytic estimator; the shrunk matrix is rebuilt with the
    explicit convex combination so it matches Eq. (lw) in the paper exactly.
    """
    X = returns.values
    alpha = float(LedoitWolf().fit(X).shrinkage_)

    Sigma_s = sample_cov(returns, annualize=annualize)
    n = Sigma_s.shape[0]
    mu_s = float(np.trace(Sigma_s) / n)
    Sigma_lw = (1.0 - alpha) * Sigma_s + alpha * mu_s * np.eye(n)
    return Sigma_lw, alpha, mu_s


def view_variance(p: np.ndarray, Sigma: np.ndarray) -> float:
    """v = p' Sigma p  for a (relative) view vector p."""
    return float(p @ Sigma @ p)


def idzorek_omega(c: float, tau: float, v_std: float) -> float:
    """omega = (1-c)/c * tau * v   (He-Litterman / Idzorek calibration)."""
    return (1.0 - c) / c * tau * v_std


def c_eff(v: float, omega: float, tau: float) -> float:
    """Effective view absorption  c_eff = tau*v / (tau*v + omega)."""
    return tau * v / (tau * v + omega)


def mu_bl_single_view(Sigma, pi, p, Q, omega, tau):
    """Black-Litterman posterior mean for a single (relative) view.

    mu* = [(tau Sigma)^-1 + p Omega^-1 p']^-1 [(tau Sigma)^-1 pi + p Omega^-1 Q]
    """
    P = p.reshape(1, -1)
    tSinv = np.linalg.inv(tau * Sigma)
    Oinv = np.array([[1.0 / omega]])
    M = tSinv + P.T @ Oinv @ P
    b = tSinv @ pi + P.T @ Oinv @ np.array([Q])
    return np.linalg.solve(M, b)

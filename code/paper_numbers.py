"""
paper_numbers.py
================

Reproduces every number quoted in the standalone paper

    "Black-Litterman with Ledoit-Wolf Shrinkage:
     Fixed-Omega Leakage and the Critical Correlation Threshold"
    (main.tex)

from the frozen offline data in data/, using the thesis pipeline in
core/ with two conventions pinned for exactness (revision 2026-07 v2):

  * Sigma_LW = (1 - alpha) * Sigma_s + alpha * mu_s * I, with alpha from
    sklearn's analytic estimator but applied to the ddof=1 sample matrix.
    (sklearn's own covariance output normalises by T and would rescale
    every closed form in the paper by (T-1)/T.)
  * Weights via w = (delta * Sigma)^{-1} mu_BL exactly as stated in the
    paper (add_M_to_Sigma=False), NOT the He-Litterman Sigma_post = Sigma + M.
    Relative views then keep the portfolios exactly fully invested.

Sections
--------
  [1] Universe, estimators, condition numbers, shrinkage intensity
  [2] Prior shift Delta pi under Idzorek-Omega (linear-in-alpha channel)
  [3] Critical-correlation table rho*(alpha, eps)  (symmetric case)
  [4] Exact per-pair Delta c_eff for all 28 pairs (heterogeneous, exact)
  [5] Clean 2x2 example IWF/IWD (both shrinkage designs)
  [6] Out-of-sample backtest Jan 2019 - Dec 2024
      (sample / LW with recomputed Omega / LW with frozen Omega,
       multi-view absorption diagnostics for Remark 1, Memmel test,
       long-only twins via cvxpy if available)

Run:  python code/paper_numbers.py     (from the repo root)

Author: Lucas Posern
"""
from __future__ import annotations

import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, str(ROOT / "core"))

from covariance import sample_cov                               # noqa: E402
from views import ViewSet                                        # noqa: E402
from black_litterman import run_blm, implied_risk_aversion      # noqa: E402
from portfolio_metrics import performance_summary               # noqa: E402

# ----------------------------------------------------------------------------
# Parameters (identical to thesis chapter 5.1)
# ----------------------------------------------------------------------------
TICKERS = ["AGG", "BWX", "IWF", "IWD", "IWO", "IWN", "EFA", "EEM"]
TRAIN = ("2008-01", "2018-12")
TEST = ("2019-01", "2024-12")
TAU = 0.05
C_CONF = 0.50
Q_VIEW = 0.03            # IWF > IWD by +3% p.a. (single view of the paper)
PAIR = ("IWF", "IWD")

DATA = ROOT.parent / "data"


def load_all():
    ret = pd.read_csv(DATA / "etf_returns_monthly.csv",
                      index_col=0, parse_dates=True)[TICKERS].dropna()
    rf = pd.read_csv(DATA / "rf_monthly.csv",
                     index_col=0, parse_dates=True).iloc[:, 0]
    rf.index = pd.to_datetime(rf.index).to_period("M").to_timestamp("M")
    rf = rf.reindex(ret.index).ffill()
    mkt = pd.read_csv(DATA / "spx_returns_monthly.csv",
                      index_col=0, parse_dates=True).iloc[:, 0]
    w_eq = pd.read_csv(DATA / "eq_weights_train_start.csv",
                       index_col=0)["w_eq"].reindex(TICKERS)
    tr = slice(*TRAIN)
    te = slice(*TEST)
    excess = ret.sub(rf, axis=0)
    return (excess.loc[tr], excess.loc[te], mkt.loc[tr], w_eq)


def kappa(S: np.ndarray) -> float:
    return float(np.linalg.cond(S))


def memmel_test(r1: pd.Series, r2: pd.Series) -> tuple[float, float]:
    """Memmel (2003) corrected Jobson-Korkie test on monthly Sharpe ratios."""
    s1 = r1.mean() / r1.std(ddof=1)
    s2 = r2.mean() / r2.std(ddof=1)
    rho = float(np.corrcoef(r1, r2)[0, 1])
    T = len(r1)
    var = (2 * (1 - rho) + 0.5 * (s1 ** 2 + s2 ** 2
           - s1 * s2 * (1 + rho ** 2))) / T
    z = (s1 - s2) / math.sqrt(var)
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return z, p


def main() -> None:
    excess_train, excess_test, mkt_train, w_eq = load_all()
    T, n = excess_train.shape
    assets = TICKERS

    # ------------------------------------------------------------------
    # [1] Estimators (annualized, on excess returns)
    #     Convention: alpha from sklearn's analytic estimator, applied to
    #     the ddof=1 sample matrix, so Sigma_LW = (1-a) Sigma_s + a mu_s I
    #     holds EXACTLY (paper Section 4).
    # ------------------------------------------------------------------
    Sigma_s = sample_cov(excess_train, annualize=True, frequency="M")
    from sklearn.covariance import LedoitWolf
    alpha = float(LedoitWolf().fit(excess_train.values).shrinkage_)
    mu_s_ann = float(np.trace(Sigma_s.values) / n)
    Sigma_lw = pd.DataFrame(
        (1 - alpha) * Sigma_s.values + alpha * mu_s_ann * np.eye(n),
        index=Sigma_s.index, columns=Sigma_s.columns)
    delta = implied_risk_aversion(mkt_train, annualize=True, frequency="M")

    print("=" * 72)
    print("[1] UNIVERSE AND ESTIMATORS (excess returns, annualized)")
    print("=" * 72)
    print(f"T = {T} months ({TRAIN[0]}..{TRAIN[1]}), n = {n}")
    print(f"delta (implied, S&P 500 excess)  = {delta:.4f}")
    print(f"alpha (analytic LW, sklearn)     = {alpha:.4f}")
    print(f"mu_s = tr(Sigma_s)/n annualized  = {mu_s_ann:.6f}"
          f"   (monthly {mu_s_ann / 12:.6f})")
    print(f"kappa(Sigma_s)  = {kappa(Sigma_s.values):7.1f}")
    print(f"kappa(Sigma_lw) = {kappa(Sigma_lw.values):7.1f}"
          f"   (reduction factor "
          f"{kappa(Sigma_s.values) / kappa(Sigma_lw.values):.1f})")

    # ------------------------------------------------------------------
    # [2] Prior shift under Idzorek-Omega
    # ------------------------------------------------------------------
    d_pi = delta * (Sigma_lw.values - Sigma_s.values) @ w_eq.values
    print()
    print("=" * 72)
    print("[2] PRIOR SHIFT Delta pi = alpha*delta*(mu_s I - Sigma_s) w_eq "
          "(annualized)")
    print("=" * 72)
    for a, dp in zip(assets, d_pi):
        print(f"  {a}: {100 * dp:+.3f} % p.a.")
    print(f"max |Delta pi| = {100 * np.max(np.abs(d_pi)):.3f} % p.a.")

    # ------------------------------------------------------------------
    # [3] rho*(alpha, eps) table (symmetric case)
    # ------------------------------------------------------------------
    # Illustrative parameters as printed in the paper table; same order of
    # magnitude as the monthly quantities in [1]. The table depends only on
    # the ratio mu_s / sigma^2 (and c).
    MU_S_ILL, SIG2_ILL, c = 4e-3, 3.5e-3, C_CONF

    def rho_star(a: float, eps: float) -> float:
        """Closed form, positive branch (Delta c_eff = +eps)."""
        return 1 - (c * a * MU_S_ILL * ((1 - c) - eps)) / (
            SIG2_ILL * (eps * (1 - c * a) + c * (1 - c) * a))

    def rho_star_neg(a: float, eps: float) -> float:
        """Negative branch (Delta c_eff = -eps); exists iff denominator > 0."""
        den = SIG2_ILL * (c * (1 - c) * a - eps * (1 - c * a))
        if den <= 0:
            return float("nan")
        return 1 - (c * a * MU_S_ILL * ((1 - c) + eps)) / den

    def dceff_sym(rho: float, a: float) -> float:
        """Exact Delta c_eff, symmetric case, omega frozen on sample."""
        v_std = 2 * SIG2_ILL * (1 - rho)
        v_lw = (1 - a) * v_std + 2 * a * MU_S_ILL
        return c * (1 - c) * (v_lw - v_std) / (c * v_lw + (1 - c) * v_std)

    alphas = [0.04, 0.10, 0.166, 0.30]
    epss = [0.01, 0.02, 0.05, 0.10]
    print()
    print("=" * 72)
    print("[3] rho*(alpha, eps)  [mu_s=4e-3, sigma^2=3.5e-3, c=0.5]")
    print("=" * 72)
    print("eps    " + "".join(f"a={a:<8}" for a in alphas))
    for eps in epss:
        row = [rho_star(a, eps) for a in alphas]
        for a, rs in zip(alphas, row):
            if -0.999 < rs < 1:
                assert abs(dceff_sym(rs, a) - eps) < 1e-10, (a, eps)
        print(f"{eps:<7}" + "".join(f"{v:<10.3f}" for v in row))
    print("(each entry verified: Delta c_eff(rho*, alpha) == eps to 1e-10)")
    print(f"negative branch (paper Sec. 6.4): rho**(0.30, 1%) = "
          f"{rho_star_neg(0.30, 0.01):.3f}, "
          f"rho**(0.041, 1%) = {rho_star_neg(0.041, 0.01):.1f} (outside [-1,1])")

    # ------------------------------------------------------------------
    # [4] Exact per-pair Delta c_eff (heterogeneous, no approximation)
    # ------------------------------------------------------------------
    # Practitioner design: shrink the FULL 8x8 covariance once, then evaluate
    # every relative view p = e_i - e_j on the same pair of matrices,
    # one view at a time. Delta c_eff is exact from v_std, v_lw.
    Ss, Sl = Sigma_s.values, Sigma_lw.values
    print()
    print("=" * 72)
    print("[4] EXACT Delta c_eff PER PAIR (full-universe shrinkage, c=0.5)")
    print("=" * 72)
    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            p = np.zeros(n)
            p[i], p[j] = 1.0, -1.0
            v_std = float(p @ Ss @ p)
            v_lw = float(p @ Sl @ p)
            dc = c * (1 - c) * (v_lw - v_std) / (c * v_lw + (1 - c) * v_std)
            rho = Ss[i, j] / np.sqrt(Ss[i, i] * Ss[j, j])
            k = np.sqrt(Ss[i, i] / Ss[j, j])
            rows.append((dc, assets[i], assets[j], rho, k))
    rows.sort(reverse=True)
    print(f"{'pair':<10}{'rho':>7}{'k=si/sj':>9}{'dc_eff':>10}")
    for dc, a, b, rho, k in rows:
        flag = "  <-- > 5pp" if dc > 0.05 else ""
        print(f"{a}/{b:<5}{rho:>7.3f}{k:>9.2f}{100 * dc:>+9.2f}pp{flag}")
    print(f"\npairs with Delta c_eff > 5pp: "
          f"{sum(1 for x in rows if x[0] > 0.05)} of {len(rows)}; "
          f"negative pairs: {sum(1 for x in rows if x[0] < 0)} "
          f"(all pairs satisfy v_std < 2 mu_s in this universe)")

    # ------------------------------------------------------------------
    # [5] Clean 2x2 example IWF/IWD
    # ------------------------------------------------------------------
    i, j = assets.index(PAIR[0]), assets.index(PAIR[1])

    def run_2x2(S2_s: np.ndarray, S2_lw: np.ndarray, tag: str) -> None:
        p = np.array([1.0, -1.0])
        w2 = np.array([0.5, 0.5])
        v_std = float(p @ S2_s @ p)
        v_lw = float(p @ S2_lw @ p)
        omega = (1 - C_CONF) / C_CONF * TAU * v_std   # frozen on sample
        ce_s = TAU * v_std / (TAU * v_std + omega)
        ce_l = TAU * v_lw / (TAU * v_lw + omega)
        pi = delta * S2_s @ w2

        def mu_bl(S: np.ndarray) -> np.ndarray:
            M = np.linalg.inv(TAU * S) + np.outer(p, p) / omega
            b = np.linalg.inv(TAU * S) @ pi + p * Q_VIEW / omega
            return np.linalg.solve(M, b)

        mu_std, mu_lw = mu_bl(S2_s), mu_bl(S2_lw)
        # return-channel isolation: same (sample) Sigma^-1 in both scenarios
        w_std = np.linalg.solve(S2_s, mu_std) / delta
        w_lw = np.linalg.solve(S2_s, mu_lw) / delta

        rho_s = S2_s[0, 1] / np.sqrt(S2_s[0, 0] * S2_s[1, 1])
        rho_l = S2_lw[0, 1] / np.sqrt(S2_lw[0, 0] * S2_lw[1, 1])
        qpp = Q_VIEW - float(p @ pi)
        dce = ce_l - ce_s

        print(f"--- {tag} ---")
        print(f"Sigma_s  = [[{S2_s[0, 0]:.5f}, {S2_s[0, 1]:.5f}], "
              f"[{S2_s[1, 0]:.5f}, {S2_s[1, 1]:.5f}]]   rho = {rho_s:.3f}")
        print(f"Sigma_lw = [[{S2_lw[0, 0]:.5f}, {S2_lw[0, 1]:.5f}], "
              f"[{S2_lw[1, 0]:.5f}, {S2_lw[1, 1]:.5f}]]   rho = {rho_l:.3f}")
        print(f"sigma_i = {np.sqrt(S2_s[0, 0]):.4f}, "
              f"sigma_j = {np.sqrt(S2_s[1, 1]):.4f}, "
              f"k = {np.sqrt(S2_s[0, 0] / S2_s[1, 1]):.3f}")
        print(f"v_std = {v_std:.5f}   v_lw = {v_lw:.5f}   omega = {omega:.4e}")
        print(f"c_eff: {100 * ce_s:.2f}% -> {100 * ce_l:.2f}%   "
              f"Delta c_eff = {100 * dce:+.2f}pp")
        print(f"q - p'pi = {qpp:.4f}")
        print(f"mu_BL_std = ({100 * mu_std[0]:.2f}%, {100 * mu_std[1]:.2f}%)"
              f"   mu_BL_lw = ({100 * mu_lw[0]:.2f}%, {100 * mu_lw[1]:.2f}%)")
        print(f"w_BL_std  = ({100 * w_std[0]:+.1f}%, {100 * w_std[1]:+.1f}%)"
              f"   w_BL_lw  = ({100 * w_lw[0]:+.1f}%, {100 * w_lw[1]:+.1f}%)")
        print(f"Delta mu = ({100 * (mu_lw[0] - mu_std[0]):+.3f}pp, "
              f"{100 * (mu_lw[1] - mu_std[1]):+.3f}pp)")
        print(f"Delta w  = ({100 * (w_lw[0] - w_std[0]):+.1f}pp, "
              f"{100 * (w_lw[1] - w_std[1]):+.1f}pp)")
        print(f"closed form 0.5*dc*(q-p'pi)     = "
              f"{100 * 0.5 * dce * qpp:+.3f}pp")
        print(f"closed form dc*(q-p'pi)/(d*v_k) = "
              f"{100 * dce * qpp / (delta * v_std):+.1f}pp")
        print()

    print()
    print("=" * 72)
    print("[5] 2x2 EXAMPLE IWF/IWD (annualized, omega frozen on sample)")
    print("=" * 72)

    # Design A (paper example): extract the 2x2 sample block, shrink it with
    # its own target mu_s^(2) = tr(S2)/2, alpha from the full universe.
    S2_s = Ss[np.ix_([i, j], [i, j])]
    mu2 = np.trace(S2_s) / 2
    S2_lw_A = (1 - alpha) * S2_s + alpha * mu2 * np.eye(2)
    run_2x2(S2_s, S2_lw_A, "Design A: shrink extracted 2x2 (own target)")

    # Design B (practitioner check): shrink full 8x8, then extract the block.
    # The universe-wide target mu_s exceeds the pair target mu_s^(2), so the
    # spread variance widens further and the shift grows (paper Sec. 6.6).
    S2_lw_B = Sl[np.ix_([i, j], [i, j])]
    run_2x2(S2_s, S2_lw_B, "Design B: extract from full-universe shrinkage")
    print(f"mu_s (universe) = {mu_s_ann:.5f}  >  mu_s^(2) (pair) = {mu2:.5f}")

    # ------------------------------------------------------------------
    # [6] Out-of-sample backtest (thesis view set)
    #     Three BLM variants: sample / LW with recomputed Idzorek-Omega /
    #     LW with Omega frozen on the sample estimator (the leakage case).
    #     Weights: w = (delta * Sigma)^{-1} mu_BL  (paper Section 2).
    # ------------------------------------------------------------------
    print("=" * 72)
    print("[6] BACKTEST Jan 2019 - Dec 2024 "
          "(constant weights, monthly rebalancing)")
    print("=" * 72)

    def build_views(Sigma_used: pd.DataFrame):
        vs = ViewSet(asset_names=assets)
        vs.add_relative(["IWF"], ["IWD"], expected_diff=0.03,
                        confidence=0.50, label="IWF>IWD +3%")
        vs.add_relative(["EEM"], ["EFA"], expected_diff=0.04,
                        confidence=0.50, label="EEM>EFA +4%")
        vs.add_relative(["IWF"], ["AGG"], expected_diff=0.06,
                        confidence=0.50, label="IWF>AGG +6%")
        return vs.build(Sigma=Sigma_used, tau=TAU, omega_method="idzorek")

    def blm_run(Sigma_used: pd.DataFrame, omega_from: pd.DataFrame):
        """BLM with Omega calibrated on omega_from (frozen if != Sigma_used)."""
        P, Q, _ = build_views(Sigma_used)
        _, _, Omega = build_views(omega_from)
        res = run_blm(Sigma_used, w_eq, P, Q, Omega, delta=delta, tau=TAU,
                      add_M_to_Sigma=False)   # w = (delta Sigma)^-1 mu_BL
        return res

    res_s = blm_run(Sigma_s, Sigma_s)
    res_l = blm_run(Sigma_lw, Sigma_lw)      # Omega recomputed (clean)
    res_f = blm_run(Sigma_lw, Sigma_s)       # Omega frozen on sample (leakage)

    for tag, res in [("sample", res_s), ("LW, recomputed Omega", res_l),
                     ("LW, frozen Omega", res_f)]:
        w = res["w_BL"]
        print(f"BLM weights ({tag}): "
              + ", ".join(f"{a} {100 * x:+.1f}%" for a, x in w.items()))
        print(f"   gross {np.abs(w.values).sum():.2f}, "
              f"net {w.values.sum():.4f}")
    print(f"max |mu_BL_lw - mu_BL_s| = "
          f"{100 * np.max(np.abs(res_l['mu_BL'].values - res_s['mu_BL'].values)):.3f} % p.a.; "
          f"prior channel max = "
          f"{100 * np.max(np.abs(res_l['pi'].values - res_s['pi'].values)):.3f} % p.a.")

    # multi-view absorption diagnostics (paper Remark 1)
    P, Q, Om_s = build_views(Sigma_s)
    _, _, Om_l = build_views(Sigma_lw)
    for tag, S, Om, res in [("sample", Sigma_s, Om_s, res_s),
                            ("LW", Sigma_lw, Om_l, res_l)]:
        Smat = P @ (TAU * S.values) @ P.T
        absorb = Smat @ np.linalg.inv(Smat + Om)
        resid = P @ res["mu_BL"].values - (
            (1 - C_CONF) * P @ res["pi"].values + C_CONF * Q)
        print(f"Remark-1 diagnostics ({tag}): absorption diag = "
              f"{np.round(np.diag(absorb), 4)}, residuals (% p.a.) = "
              f"{np.round(100 * resid, 3)}")

    ports = {
        "Market equilibrium": excess_test @ w_eq.values,
        "BLM (sample)": excess_test @ res_s["w_BL"].values,
        "BLM (LW, Omega recomputed)": excess_test @ res_l["w_BL"].values,
        "BLM (LW, Omega frozen)": excess_test @ res_f["w_BL"].values,
    }

    # long-only twins (cvxpy QP, as in thesis section 5.1)
    try:
        from optimize import mv_weights
        for name, (S, mu) in {
            "BLM (sample, long-only)": (Sigma_s, res_s["mu_BL"]),
            "BLM (LW, long-only)": (Sigma_lw, res_l["mu_BL"]),
        }.items():
            w_lo = mv_weights(mu, S, delta=delta, long_only=True)
            ports[name] = excess_test @ (
                w_lo.values if hasattr(w_lo, "values") else np.asarray(w_lo))
    except Exception as e:  # pragma: no cover
        print(f"(long-only skipped: {type(e).__name__}: {e})")

    print()
    summary = performance_summary(ports, frequency=12)
    with pd.option_context("display.float_format", "{:8.4f}".format,
                           "display.width", 120):
        print(summary.to_string())

    z, pval = memmel_test(ports["BLM (sample)"],
                          ports["BLM (LW, Omega recomputed)"])
    print(f"\nSharpe difference sample vs LW (recomputed Omega): "
          f"Memmel z = {z:.2f}, p = {pval:.2f}  -> not significant")


if __name__ == "__main__":
    main()

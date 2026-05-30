"""
example_positive_iwf_iwd.py
===========================
Positive-drift example of Section "Real Data".

IWF (US Large-Cap Growth) vs IWD (US Large-Cap Value), a highly correlated
style pair (rho ~ 0.92). Freezing omega on the sample estimator and switching
to Ledoit-Wolf raises the effective absorption c_eff from 50% to ~59%.

Downloads live monthly data from Yahoo Finance (2008-2018) and runs the full
chain  Delta c_eff -> Delta mu -> Delta w.  Numbers match the boxed result in
the report to within rounding / data-vintage noise.

The weight step deliberately uses Sigma_s^-1 for both scenarios so that only
the view-absorption channel (not the change in Sigma^-1) shows up.
"""

import numpy as np
from blm_core import (
    download_monthly_returns, sample_cov, ledoit_wolf,
    view_variance, idzorek_omega, c_eff, mu_bl_single_view,
)

# ---- Inputs (report Section "Positive Drift") ----
TICKERS = ["IWF", "IWD"]
START, END = "2008-01-01", "2018-12-31"
DELTA = 2.44
TAU = 0.05
C = 0.50
W_EQ = np.array([0.50, 0.50])
Q = 0.03                      # IWF beats IWD by 3% p.a.
p = np.array([1.0, -1.0])     # relative view  p = e_IWF - e_IWD


def main():
    returns = download_monthly_returns(TICKERS, START, END)
    print(f"Window: {returns.index[0].date()} .. {returns.index[-1].date()} "
          f"({len(returns)} months)\n")

    Sigma_s = sample_cov(returns)
    Sigma_lw, alpha, mu_s = ledoit_wolf(returns)

    rho_s = Sigma_s[0, 1] / np.sqrt(Sigma_s[0, 0] * Sigma_s[1, 1])
    rho_lw = Sigma_lw[0, 1] / np.sqrt(Sigma_lw[0, 0] * Sigma_lw[1, 1])

    v_std = view_variance(p, Sigma_s)
    v_lw = view_variance(p, Sigma_lw)

    # omega calibrated once on the sample estimator, then frozen
    omega = idzorek_omega(C, TAU, v_std)

    ce_s = c_eff(v_std, omega, TAU)
    ce_lw = c_eff(v_lw, omega, TAU)

    pi = DELTA * Sigma_s @ W_EQ
    mu_s_vec = mu_bl_single_view(Sigma_s, pi, p, Q, omega, TAU)
    mu_lw_vec = mu_bl_single_view(Sigma_lw, pi, p, Q, omega, TAU)

    # return-channel isolation: same Sigma_s^-1 in both weight steps
    w_s = np.linalg.solve(Sigma_s, mu_s_vec) / DELTA
    w_lw = np.linalg.solve(Sigma_s, mu_lw_vec) / DELTA

    print(f"alpha (Ledoit-Wolf) = {alpha:.3f}   mu_s = {mu_s:.5f}")
    print(f"rho_s = {rho_s:.3f}   rho_LW = {rho_lw:.3f}")
    print(f"v_std = {v_std:.5f}   v_lw = {v_lw:.5f}   omega = {omega:.3e}\n")

    print(f"c_eff(Sample) = {ce_s*100:.2f}%   c_eff(LW) = {ce_lw*100:.2f}%   "
          f"Delta c_eff = {(ce_lw-ce_s)*100:+.2f} pp\n")

    print(f"Delta mu (pp) = {(mu_lw_vec-mu_s_vec)*100}")
    print(f"Delta w  (pp) = {(w_lw-w_s)*100}\n")

    # closed-form cross-check, Eq. (deltas)
    qmp = Q - p @ pi
    dmu_cf = 0.5 * (ce_lw - ce_s) * qmp
    dw_cf = (ce_lw - ce_s) * qmp / (DELTA * v_std)
    print("Closed-form check (Eq. deltas):")
    print(f"  Delta mu_IWF ~ {dmu_cf*100:+.3f} pp   (numeric {(mu_lw_vec-mu_s_vec)[0]*100:+.3f} pp)")
    print(f"  Delta w_IWF  ~ {dw_cf*100:+.1f} pp    (numeric {(w_lw-w_s)[0]*100:+.1f} pp)")


if __name__ == "__main__":
    main()

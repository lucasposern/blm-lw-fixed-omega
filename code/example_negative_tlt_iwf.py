"""
example_negative_tlt_iwf.py
===========================
Negative-drift example of Section "Real Data".

TLT (iShares 20+ Year Treasury) vs IWF (US Large-Cap Growth), a flight-to-safety
pair with rho ~ -0.35 over 2008-2018. Here the view variance is large
(v_std > 2*mu_s), so freezing omega and switching to Ledoit-Wolf LOWERS the
effective absorption: the model trusts the view slightly less than specified.

Downloads live monthly data from Yahoo Finance and reports the sign and size of
Delta c_eff, contrasting it with the positive case.
"""

import numpy as np
from blm_core import (
    download_monthly_returns, sample_cov, ledoit_wolf,
    view_variance, idzorek_omega, c_eff,
)

TICKERS = ["TLT", "IWF"]
START, END = "2008-01-01", "2018-12-31"
TAU = 0.05
C = 0.50
p = np.array([1.0, -1.0])     # view  TLT > IWF


def main():
    returns = download_monthly_returns(TICKERS, START, END)
    print(f"Window: {returns.index[0].date()} .. {returns.index[-1].date()} "
          f"({len(returns)} months)\n")

    Sigma_s = sample_cov(returns)
    Sigma_lw, alpha, mu_s = ledoit_wolf(returns)

    rho = Sigma_s[0, 1] / np.sqrt(Sigma_s[0, 0] * Sigma_s[1, 1])
    sig_tlt = np.sqrt(Sigma_s[0, 0])
    sig_iwf = np.sqrt(Sigma_s[1, 1])

    v_std = view_variance(p, Sigma_s)
    v_lw = view_variance(p, Sigma_lw)
    omega = idzorek_omega(C, TAU, v_std)   # frozen on the sample estimator

    ce_s = c_eff(v_std, omega, TAU)
    ce_lw = c_eff(v_lw, omega, TAU)

    print(f"rho(TLT, IWF) = {rho:+.3f}")
    print(f"sigma(TLT) = {sig_tlt*100:.1f}%   sigma(IWF) = {sig_iwf*100:.1f}%")
    print(f"alpha (LW) = {alpha:.3f}   mu_s = {mu_s:.4f}\n")

    print(f"v_std = {v_std:.4f}   2*mu_s = {2*mu_s:.4f}   "
          f"2*mu_s - v_std = {2*mu_s - v_std:+.4f}  (<0  =>  Delta c_eff < 0)")
    print(f"omega (x1e-4) = {omega*1e4:.2f}   v_lw = {v_lw:.4f}\n")

    print(f"c_eff(Sample) = {ce_s*100:.1f}%   c_eff(LW) = {ce_lw*100:.1f}%   "
          f"Delta c_eff = {(ce_lw-ce_s)*100:+.2f} pp")
    print("\nNegative correlation => view absorbed LESS than specified, "
          "and the effect is mild (large v damps the shift).")


if __name__ == "__main__":
    main()

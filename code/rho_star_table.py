"""
rho_star_table.py
=================
Reproduces Table (rho-star) of the report: the critical correlation

    rho*(alpha, eps) = 1 - alpha*mu_s*[c(1-c) - eps]
                           / ( sigma^2 * [c(1-c)*alpha + eps*(2 - alpha)] )

A symmetric relative view triggers |Delta c_eff| > eps once rho_ij > rho*.
Pure analytic, no data download required.
"""

import numpy as np

# Report parameters (Section "Example Values")
MU_S = 0.004      # average asset variance (monthly order of magnitude)
SIGMA2 = 0.0035   # per-asset variance
C = 0.50          # stated confidence

ALPHAS = [0.04, 0.10, 0.166, 0.30]
EPSILONS = [0.01, 0.02, 0.05, 0.10]


def rho_star(alpha, eps, mu_s=MU_S, sigma2=SIGMA2, c=C):
    num = alpha * mu_s * (c * (1 - c) - eps)
    den = sigma2 * (c * (1 - c) * alpha + eps * (2 - alpha))
    return 1.0 - num / den


def main():
    print(f"rho*(alpha, eps)   mu_s={MU_S}, sigma^2={SIGMA2}, c={C}\n")
    header = "  eps \\ alpha |" + "".join(f"  {a:>6}" for a in ALPHAS)
    print(header)
    print("  " + "-" * (len(header) - 2))
    for eps in EPSILONS:
        row = "".join(f"  {rho_star(a, eps):6.3f}" for a in ALPHAS)
        print(f"  {eps*100:6.0f}%    |{row}")

    print(
        "\nReading: at alpha=0.04, eps=5%, any symmetric pair with rho > 0.915 "
        "\ndrifts c_eff by more than 5 percentage points."
    )


if __name__ == "__main__":
    main()

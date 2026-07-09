"""
rho_star_table.py
=================
Reproduces the critical-correlation table of the paper (main.tex):

    rho*(alpha, eps) = 1 - c*alpha*mu_s*[(1-c) - eps]
                           / ( sigma^2 * [eps*(1 - c*alpha) + c*(1-c)*alpha] )

with eps in (0, 1-c). A symmetric relative view (sigma_i = sigma_j) triggers
|Delta c_eff| > eps once rho_ij > rho*.

NOTE (2026-07): this replaces an earlier, incorrect version of the formula
(a factor-2 error shifted every row of the table by one eps level). Each
entry below is verified by numerically checking the exact Delta c_eff.

Pure analytic, no data required. For every number in the paper, including
the exact per-pair evaluation without the symmetry assumption, run
paper_numbers.py instead.
"""

# Illustrative parameters (Section "Values and Economic Statement");
# only the ratio mu_s/sigma^2 and c matter.
MU_S = 0.004      # average asset variance (monthly order of magnitude)
SIGMA2 = 0.0035   # per-asset variance
C = 0.50          # stated confidence

ALPHAS = [0.04, 0.10, 0.166, 0.30]
EPSILONS = [0.01, 0.02, 0.05, 0.10]


def rho_star(alpha, eps, mu_s=MU_S, sigma2=SIGMA2, c=C):
    """Corrected closed form; eps must lie in (0, 1-c)."""
    num = c * alpha * mu_s * ((1 - c) - eps)
    den = sigma2 * (eps * (1 - c * alpha) + c * (1 - c) * alpha)
    return 1.0 - num / den


def dceff_symmetric(rho, alpha, mu_s=MU_S, sigma2=SIGMA2, c=C):
    """Exact Delta c_eff in the symmetric case, omega frozen on sample."""
    v_std = 2 * sigma2 * (1 - rho)
    v_lw = (1 - alpha) * v_std + 2 * alpha * mu_s
    return c * (1 - c) * (v_lw - v_std) / (c * v_lw + (1 - c) * v_std)


def main():
    print(f"rho*(alpha, eps)   mu_s={MU_S}, sigma^2={SIGMA2}, c={C}\n")
    header = "  eps \\ alpha |" + "".join(f"  {a:>6}" for a in ALPHAS)
    print(header)
    print("  " + "-" * (len(header) - 2))
    for eps in EPSILONS:
        vals = []
        for a in ALPHAS:
            rs = rho_star(a, eps)
            if -0.999 < rs < 1.0:  # verify against the exact expression
                assert abs(dceff_symmetric(rs, a) - eps) < 1e-12, (a, eps)
            vals.append(rs)
        row = "".join(f"  {v:6.3f}" for v in vals)
        print(f"  {eps*100:6.0f}%    |{row}")

    print(
        "\nReading: at alpha=0.04, eps=5%, any symmetric pair with rho > 0.826"
        "\ndrifts c_eff by more than 5 percentage points. Every entry is"
        "\nverified against the exact Delta c_eff (assertion above)."
    )


if __name__ == "__main__":
    main()

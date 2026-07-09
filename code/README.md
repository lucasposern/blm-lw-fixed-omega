# Code for "Black-Litterman with Ledoit-Wolf Shrinkage: Fixed-Ω Leakage and the Critical Correlation Threshold"

Reproduces every number in `main.tex` (the current paper draft) from the
frozen offline data in `../data/` -- no network access required.

**Main entry point:**

```bash
pip install -r requirements.txt   # plus: pandas scikit-learn cvxpy
python paper_numbers.py           # reproduces every number in main.tex
```

`paper_numbers.py` uses the exact thesis pipeline (excess returns, sample
covariance with T−1 normalisation, analytic sklearn Ledoit-Wolf, implied
risk aversion from S&P 500 excess returns) via the modules in `core/`
(frozen copies of the thesis code). Sections: estimators and condition
numbers, prior shift Δπ, corrected ρ*(α, ε) table, exact per-pair Δc_eff
for all 28 pairs, the 2x2 IWF/IWD example (both shrinkage designs), and
the out-of-sample backtest 2019-2024.

## Files

| File | Reproduces | Data |
|------|------------|------|
| `paper_numbers.py` | **every number in `main.tex`** | `../data/` (offline) |
| `core/` | frozen thesis pipeline (covariance, views, BLM, metrics, QP) | — |
| `rho_star_table.py` | standalone ρ\*(α, ε) table, corrected formula | offline |
| `blm_core.py` | legacy helpers (old report) | — |
| `example_positive_iwf_iwd.py` | legacy: positive drift example (old report) | yfinance |
| `example_negative_tlt_iwf.py` | legacy: negative drift TLT/IWF (old report) | yfinance |

## History / caveats

- **2026-07 revision.** `main.tex` replaces `tau_omega_independence_BLM.tex`
  (old draft, kept for reference). The old draft and the pre-revision
  `rho_star_table.py` contained an incorrect ρ* formula (factor-2 error:
  every table row shifted by one ε level). Corrected and verified against
  the exact Δc_eff; do not quote the old table.
- The legacy yfinance examples use total returns and live downloads; the
  paper pipeline uses frozen excess returns. Numbers can differ in the last
  digit. The TLT/IWF negative-drift example is not part of the current
  paper (TLT is outside the Idzorek universe) but is kept as a candidate
  for a sign-reversal illustration in a later revision.

## Method notes

- **Ledoit-Wolf.** α from the analytic sklearn estimator
  (`LedoitWolf().fit(X).shrinkage_`), applied on monthly excess returns.
- **Frozen Ω.** ω = (1−c)/c · τ · pᵀΣ_s p is calibrated once on the sample
  estimator and kept fixed when Σ switches to Ledoit-Wolf; that frozen ω is
  the source of the c_eff drift.
- **Channel isolation.** The 2x2 example uses Σ_s⁻¹ in the weight step for
  both scenarios, so the reported Δw reflects only the view-absorption
  channel, not the change in Σ⁻¹ (treated separately in the paper).

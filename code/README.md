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

- **2026-07 extension pass (`main_v2.tex`).** Edge cases closed: general
  view vectors (v_LW = (1−α)v + α·μ_s·‖p‖², criterion v ≶ μ_s‖p‖²;
  absolute view AGG +12.9pp = largest leakage, EEM −0.5pp = realised
  deflation, basket +5.8pp), strict monotonicity in v with range
  (−c(1−c)α/(1−cα), 1−c) incl. deflation floor (answers ρ=−1),
  PD-boundary notes (ρ=±1 as continuous limits), worst-case confidence
  c* = 1/(1+√(v_LW/v)). Numbers in section [7] of `paper_numbers.py`.
  `main.tex` = referee-fix state without extensions.
- **2026-07 revision 2 (referee fixes).** Two conventions pinned for
  exactness: (i) Σ_LW = (1−α)Σ_s + α·μ_s·I with sklearn's analytic α
  applied to the ddof=1 sample matrix (sklearn's own covariance output
  normalises by T and would rescale every closed form by (T−1)/T);
  (ii) weights via w = (δΣ)⁻¹μ_BL exactly as in the paper
  (`add_M_to_Sigma=False`), not He-Litterman's Σ_post = Σ + M.
  Backtest gains a third BLM column (LW with frozen Ω = the leakage case),
  multi-view absorption diagnostics (paper Remark 1) and a Memmel (2003)
  Sharpe-difference test. Numbers changed at the 0.1–0.5pp level
  (e.g. IWF/IWD pair shift +10.6 → +10.7pp; backtest 17.78 → 18.06%).
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
  (`LedoitWolf().fit(X).shrinkage_`) on monthly excess returns; the
  shrinkage is then applied to the ddof=1 sample matrix so that
  Σ_LW = (1−α)Σ_s + α·μ_s·I holds exactly (see History above).
- **Weights.** w = (δΣ)⁻¹μ_BL throughout; relative views keep the
  portfolios exactly fully invested (net weight 1.0000, printed).
- **Frozen Ω.** ω = (1−c)/c · τ · pᵀΣ_s p is calibrated once on the sample
  estimator and kept fixed when Σ switches to Ledoit-Wolf; that frozen ω is
  the source of the c_eff drift. The backtest now runs this workflow
  out of sample as a third BLM column.
- **Channel isolation.** The 2x2 example uses Σ_s⁻¹ in the weight step for
  both scenarios, so the reported Δw reflects only the view-absorption
  channel, not the change in Σ⁻¹ (treated separately in the paper).
- **Performance conventions.** Arithmetic annualisation (mean ×12,
  vol ×√12) on monthly excess returns; max drawdown on the compounded
  excess-return path; no transaction or shorting costs.

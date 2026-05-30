# Code for "When the Independence of τ and Ω Backfires"

Reproduces every number in `tau_omega_independence_BLM.tex`. Pure-analytic
parts run offline; the two real-data examples pull live monthly prices from
Yahoo Finance.

## Files

| File | Reproduces | Data |
|------|------------|------|
| `blm_core.py` | shared helpers (sample cov, Ledoit-Wolf, c_eff, single-view posterior) | — |
| `rho_star_table.py` | the ρ\*(α, ε) table (Section "Example Values") | offline |
| `example_positive_iwf_iwd.py` | positive drift, Δc_eff ≈ +9.4 pp (Section "Real Data") | yfinance |
| `example_negative_tlt_iwf.py` | negative drift, Δc_eff ≈ −1.4 pp (Section "Real Data") | yfinance |

## Setup

```bash
pip install -r requirements.txt
```

## Run

```bash
python rho_star_table.py            # offline, deterministic
python example_positive_iwf_iwd.py  # downloads IWF, IWD  (2008-2018)
python example_negative_tlt_iwf.py  # downloads TLT, IWF  (2008-2018)
```

## Method notes

- **Ledoit-Wolf.** The shrinkage intensity α comes from the analytic sklearn
  estimator (`LedoitWolf().fit(X).shrinkage_`); the shrunk matrix is rebuilt as
  the explicit convex combination `Σ_LW = (1−α)Σ_s + α·μ_s·I` with
  `μ_s = tr(Σ_s)/n`, matching Eq. (lw) in the paper.
- **Frozen Ω.** Each example calibrates `ω = (1−c)/c · τ · pᵀΣ_s p` once on the
  sample estimator and keeps it fixed when Σ switches to Ledoit-Wolf. That frozen
  ω is the source of the c_eff drift.
- **Channel isolation.** The positive example uses `Σ_s⁻¹` in the weight step for
  both scenarios, so the reported Δw reflects only the view-absorption channel,
  not the change in Σ⁻¹.
- **Excess vs. total returns.** The thesis computes the covariance on excess
  returns. These scripts use total returns; subtracting a near-constant monthly
  risk-free rate moves the covariance only marginally, so the figures may differ
  in the last reported digit depending on the data vintage.

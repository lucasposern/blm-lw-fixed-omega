# Black-Litterman with Ledoit-Wolf Shrinkage: Fixed-Ω Leakage and the Critical Correlation Threshold

Lucas Posern (TU Dresden) — standalone paper, grown out of my bachelor thesis.

**Claim in one paragraph.** Whether switching the covariance estimator inside
Black-Litterman from the sample estimator to Ledoit-Wolf is harmless depends on
how the view-uncertainty matrix Ω is handled. If Ω follows Idzorek's confidence
calibration, the switch is exactly inert in the view channel for a single view.
If Ω is calibrated once on the sample covariance and then frozen, the effective
view confidence shifts by an exact, τ-free closed form; the sign criterion is
`v < 2·μ_s` (view variance vs. twice the average asset variance), which catches
highly correlated style pairs *and* low-variance bond pairs. On the eight-asset
Idzorek universe (α = 0.041), five of 28 pairs shift by more than 5 pp, IWF/IWD
by +10.7 pp. Out of sample, the frozen-Ω workflow forfeits close to half of the
de-risking the estimator switch was meant to deliver.

## Repository layout

| Path | Content |
|------|---------|
| `fixed_omega.tex` / `fixed_omega.pdf` | The paper (compiled on Overleaf) |
| `code/paper_numbers.py` | Reproduces **every number** in the paper from `data/` |
| `code/core/` | Frozen thesis pipeline (covariance, views, BLM, metrics, QP) |
| `code/README.md` | Conventions, method notes, changelog |
| `data/` | Frozen monthly excess-return inputs (Yahoo Finance vintage) |

## Reproduce

```bash
pip install -r code/requirements.txt
python code/paper_numbers.py
```

Runs offline from the frozen CSVs. Two conventions are pinned so that every
closed form in the paper holds exactly in the pipeline: Σ_LW = (1−α)Σ_s + α·μ_s·I
with sklearn's analytic α applied to the ddof=1 sample matrix, and weights via
w = (δΣ)⁻¹μ_BL. Details in `code/README.md`.

## License

Code: MIT (see `LICENSE`). Paper text: © Lucas Posern.

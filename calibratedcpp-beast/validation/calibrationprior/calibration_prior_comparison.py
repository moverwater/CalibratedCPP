"""
Compare BEAST MCMC samples vs LPhy direct simulation for the calibration prior.

Usage (from the validation/calibrationprior directory):
    python calibration_prior_comparison.py

Requires: numpy, scipy, pandas, matplotlib
"""

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend

BEAST_LOG = "test_calibration_prior.log"
LPHY_TSV  = "lphy_wsim.tsv"
BURNIN    = 0.10
ALPHA     = 0.05
OUT_PDF   = "mcmc_vs_lphy_comparison.pdf"

CLADES = [f"mrca.clade{i}" for i in range(1, 11)]

# Bounds in clade1…clade10 order (clade1 = root, matching calibrationprior_simulation.r)
LOWER = [10.0, 9.4, 4.8, 4.0, 8.0, 6.8, 1.8, 2.0, 1.7, 1.6]
UPPER = [10.5, 9.6, 5.0, 5.0, 9.5, 8.0, 2.5, 3.2, 2.5, 2.5]

# ── Load data ─────────────────────────────────────────────────────────────────
beast_raw = pd.read_csv(BEAST_LOG, sep="\t", comment="#")
n_burn    = int(len(beast_raw) * BURNIN)
beast     = beast_raw.iloc[n_burn:].reset_index(drop=True)

lphy = pd.read_csv(LPHY_TSV, sep="\t")

# ── Coverage + t-tests ───────────────────────────────────────────────────────
print(f"\n{'Clade':<16} {'BEAST-cov':>9} {'LPhy-cov':>9} {'t p-val':>9} {'t':>5}")
for col, lo, hi in zip(CLADES, LOWER, UPPER):
    b = beast[col].values
    l = lphy[col].values
    b_cov = np.mean((b >= lo) & (b <= hi))
    l_cov = np.mean((l >= lo) & (l <= hi))
    t_p   = stats.ttest_ind(b, l).pvalue
    print(f"{col:<16} {b_cov:>9.3f} {l_cov:>9.3f} "
          f"{t_p:>9.4f} {'PASS' if t_p > ALPHA else 'FAIL':>5}")

# ── Overlaid density plots ────────────────────────────────────────────────────
fig, axes = plt.subplots(5, 2, figsize=(6.5, 9))
axes = axes.flatten()

for k, (col, lo, hi) in enumerate(zip(CLADES, LOWER, UPPER)):
    ax  = axes[k]
    b   = beast[col].values
    l   = lphy[col].values
    b_cov = np.mean((b >= lo) & (b <= hi))
    l_cov = np.mean((l >= lo) & (l <= hi))
    clade_num = col.replace("mrca.clade", "")

    # KDE
    x_min = min(b.min(), l.min())
    x_max = max(b.max(), l.max())
    xs    = np.linspace(x_min, x_max, 500)

    for values, label, color in [
        (b, "BEAST MCMC", "#E69F00"),
        (l, "LPhy sim",   "#56B4E9"),
    ]:
        kde = stats.gaussian_kde(values)
        ax.fill_between(xs, kde(xs), alpha=0.35, color=color, label=label)
        ax.plot(xs, kde(xs), color=color, linewidth=0.9)

    ax.axvline(lo, color="black", linestyle="--", linewidth=0.7)
    ax.axvline(hi, color="black", linestyle="--", linewidth=0.7)
    ax.set_title(f"Clade {clade_num}  BEAST={b_cov:.3f}  LPhy={l_cov:.3f}",
                 fontsize=8, fontweight="bold")
    ax.set_xlabel(f"Clade {clade_num}", fontsize=8)
    ax.set_ylabel("Density", fontsize=8)
    ax.tick_params(labelsize=7)
    if k == 0:
        ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig(OUT_PDF)
print(f"\nSaved: {OUT_PDF}")
plt.show()

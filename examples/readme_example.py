"""Reproduce the README example and figures."""
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from badlm import BADLM, readme_lag_curve, simulate_dlm

here = os.path.dirname(os.path.abspath(__file__))
figdir = os.path.join(here, "..", "figures")
os.makedirs(figdir, exist_ok=True)

truth = readme_lag_curve(50)
x, y = simulate_dlm(truth, n=500, standardize_x=True, ar_eps=0.0, sd_eps=0.01, seed=2018)

model = BADLM(nlag=50, k=30)
t0 = time.time()
model.fit(x, y, samples=1000, warmup=1000, chains=4, seed=1)
print(f"fitted in {time.time() - t0:.1f}s")
print(model.summary())

fig, ax = plt.subplots(figsize=(7, 4))
model.plot(ax=ax, truth=truth)
fig.tight_layout(); fig.savefig(os.path.join(figdir, "readme_fit.png"), dpi=120)

# compare the adaptive model with the non-adaptive one (M1 in the paper)
fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for ax, tm, name in zip(axes, ("constant", "rw"), ("non-adaptive (M1)", "adaptive (M3)")):
    m = BADLM(nlag=50, k=30, tau_model=tm).fit(x, y, samples=1000, warmup=1000, chains=2, seed=1)
    m.plot(ax=ax, truth=truth, title=f"{name}: ED = {m.effective_df()['mean']:.1f}")
fig.tight_layout(); fig.savefig(os.path.join(figdir, "readme_adaptive_vs_not.png"), dpi=120)

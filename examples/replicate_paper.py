"""Reduced replication of Section 3.1 of the paper: the adaptive model (M3)
against the non-adaptive model (M1) on five lag-curve shapes.

Usage: python examples/replicate_paper.py [n_reps]

The published curves are not available in closed form, so
``badlm.paper_lag_curves`` provides approximations with the same features.
RMSE and squared bias are scaled by 1e3 as in Table 1 of the paper.
"""
import sys
import time

import numpy as np

from badlm import BADLM, paper_lag_curves, simulate_dlm

n_reps = int(sys.argv[1]) if len(sys.argv) > 1 else 10
curves = paper_lag_curves(50)
models = {"M1": dict(tau_model="constant"), "M3": dict(tau_model="rw")}
results = {}
t0 = time.time()
for cname, truth in curves.items():
    for mname, cfg in models.items():
        rmse, bias, ed = [], [], []
        for r in range(n_reps):
            x, y = simulate_dlm(truth, n=500, seed=1000 * r + 7)
            m = BADLM(nlag=50, k=34, **cfg).fit(x, y, samples=500, warmup=500, chains=2, seed=r)
            est = m.summary_.mean
            rmse.append(np.sqrt(np.mean((est - truth) ** 2)))
            bias.append(np.mean(est - truth) ** 2)
            ed.append(m.effective_df()["mean"])
        results[(cname, mname)] = (1e3 * np.mean(rmse), 1e3 * np.mean(bias), np.mean(ed))
        print(f"{cname:15s} {mname}: RMSE {1e3*np.mean(rmse):6.2f}  Bias2 {1e3*np.mean(bias):6.2f}  "
              f"ED {np.mean(ed):5.2f}   [{time.time()-t0:.0f}s]", flush=True)

print("\n| Lag curve | Summary | M1 | M3 |")
print("|---|---|---|---|")
for cname in curves:
    for i, lab in enumerate(("RMSE", "Bias2", "ED")):
        print(f"| {cname} | {lab} | {results[(cname, 'M1')][i]:.2f} | {results[(cname, 'M3')][i]:.2f} |")

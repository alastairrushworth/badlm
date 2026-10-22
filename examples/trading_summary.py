"""Figure summarising the trading experiments from the saved JSON results:
the volatility-dependent lag surface and the Sharpe ratio of each signal.
Run after trading_pooled.py and trading_modulated.py."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
resdir = os.path.join(here, "..", "paper", "results")
figdir = os.path.join(here, "..", "figures")


def load(name):
    with open(os.path.join(resdir, name)) as f:
        return json.load(f)


mod = load("trading_modulated_fx.json")
pooled = load("trading_pooled_fx.json")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
ax = axes[0]
surf = np.array(mod["volatility"]["surface"])
lags = np.arange(1, surf.shape[1] + 1)
labels = ["very low", "low", "normal", "high", "very high"]
for i, (s, lab) in enumerate(zip(surf, labels)):
    ax.plot(lags, s, color=plt.cm.coolwarm(i / 4), lw=1.6, label=f"{lab} volatility")
ax.plot(lags, mod["static"]["curve"], color="k", ls="--", lw=1, label="static curve")
ax.axhline(0, color="grey", lw=0.5)
ax.set_xscale("log")
ax.set_xlabel("lag in days (log scale)")
ax.set_ylabel("coefficient (scaled return per scaled return)")
ax.set_title("(a) currency lag curve by volatility state (final refit)")
ax.legend(frameon=False, fontsize=8)

ax = axes[1]
rows = []
for r in pooled["results"]:
    if r["signal"] in ("badlm curve", "TSMOM linear(252d)", "sum21 + sum252") and r["rule"] == "gp 0.1":
        rows.append((r["signal"] if r["signal"] != "badlm curve" else "static curve", r))
rows.append(("volatility-modulated", [r for r in mod["volatility"]["results"] if r["rule"] == "gp 0.1"][0]))
names = [n for n, _ in rows]
y = np.arange(len(rows))
for j, (cost, col) in enumerate(((0.0, "C0"), (0.0005, "C1"))):
    key = "sharpe_gross" if cost == 0 else "net_5bp"
    ax.barh(y + (j - 0.5) * 0.38, [r[key] for _, r in rows], height=0.36, color=col,
            label="gross" if cost == 0 else "net of 5 bp")
ax.set_yticks(y)
ax.set_yticklabels(names, fontsize=9)
ax.invert_yaxis()
ax.set_xlabel("annualised Sharpe ratio, 1996-2026, slow trading (gp 0.1)")
ax.set_title("(b) currencies: out-of-sample Sharpe ratio by signal")
ax.legend(frameon=False, fontsize=8, loc="lower right")

fig.tight_layout()
os.makedirs(figdir, exist_ok=True)
fig.savefig(os.path.join(figdir, "trading_summary.png"), dpi=130)
print("figure written")

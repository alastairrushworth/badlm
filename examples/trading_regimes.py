"""Experiment 4: a latent regime on the amplitude of the pooled curve.

The walk-forward predictions ``m_t`` of the static curve are the signal; a
two-state hidden Markov model lets the amplitude of the signal switch between
regimes, ``y_t = a_s m_t + e_t`` with ``e_t ~ N(0, sigma^2)``, the state being
common to every asset on a day. The transition matrix carries a sticky prior
(pseudo-counts on the diagonal) so that the states are persistent regimes
rather than day-to-day dispersion. Parameters are re-estimated by EM at each
refit on training days only, and the filtered state probabilities scale the
aim out of sample: ``aim_t = E[a_s | data to t] m_t``.

Run: python examples/trading_regimes.py [--universe fx] [options]
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _markets import load_panel  # noqa: E402
from _trading import evaluate, ic, print_table, stacked_rows, strip  # noqa: E402

from badlm import BADLM  # noqa: E402
from badlm.backtest import walk_forward  # noqa: E402
from badlm.priors import Priors  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
resdir = os.path.join(here, "..", "paper", "results")


def day_stats(rows):
    """Per-day sufficient statistics of y on m: sum m^2, sum m y, sum y^2, n."""
    g = rows.assign(mm=rows["pred"] ** 2, my=rows["pred"] * rows["y"], yy=rows["y"] ** 2).groupby("time")
    s = g[["mm", "my", "yy"]].sum()
    s["n"] = g.size()
    return s


def emission(stats, a, s2):
    """Log-likelihood of each day under amplitude ``a`` and variance ``s2``."""
    rss = stats["yy"].to_numpy() - 2 * a * stats["my"].to_numpy() + a**2 * stats["mm"].to_numpy()
    return -0.5 * stats["n"].to_numpy() * np.log(2 * np.pi * s2) - 0.5 * rss / s2


def _scaled_emissions(stats, a, s2):
    """Emission densities scaled by each day's maximum over states."""
    logB = np.stack([emission(stats, a[k], s2[k]) for k in range(a.size)], axis=1)
    return np.exp(logB - logB.max(axis=1, keepdims=True))


def fit_hmm(stats, n_iter=20, K=2, sticky=250.0):
    """EM for a K-state HMM with regression emissions, a common variance and a
    sticky prior (``sticky`` pseudo-counts on the diagonal of the transition
    matrix), run with the scaled forward-backward recursions."""
    a = np.linspace(0.0, 1.5, K)
    s2 = np.ones(K)
    P = np.full((K, K), 0.005) + np.eye(K) * (1 - 0.005 * K)
    pi = np.full(K, 1.0 / K)
    T = len(stats)
    mm, my, yy, n = (stats[c].to_numpy() for c in ("mm", "my", "yy", "n"))
    for _ in range(n_iter):
        B = _scaled_emissions(stats, a, s2)
        alpha = np.zeros((T, K)); beta = np.ones((T, K)); c = np.zeros(T)
        alpha[0] = pi * B[0]; c[0] = alpha[0].sum(); alpha[0] /= c[0]
        for t in range(1, T):
            v = (alpha[t - 1] @ P) * B[t]
            c[t] = v.sum()
            alpha[t] = v / c[t]
        for t in range(T - 2, -1, -1):
            beta[t] = (P @ (B[t + 1] * beta[t + 1])) / c[t + 1]
        g = alpha * beta
        g /= g.sum(axis=1, keepdims=True)
        xi = np.einsum("ti,tj->ij", alpha[:-1], (B[1:] * beta[1:]) / c[1:, None]) * P
        xi += sticky * np.eye(K)
        pi = g[0]
        P = xi / xi.sum(axis=1, keepdims=True)
        rss = 0.0
        for k in range(K):
            w = g[:, k]
            a[k] = (w * my).sum() / max((w * mm).sum(), 1e-12)
            rss += (w * (yy - 2 * a[k] * my + a[k] ** 2 * mm)).sum()
        s2[:] = rss / n.sum()
    return {"a": a, "s2": s2, "P": P, "pi": pi, "gamma": g}


def filter_forward(stats, par, start_prob=None):
    """Filtered state probabilities P(s_t | days <= t) for every day in stats."""
    B = _scaled_emissions(stats, par["a"], par["s2"])
    out = np.zeros((len(stats), par["a"].size))
    prev = par["pi"] if start_prob is None else start_prob
    for t in range(len(stats)):
        v = (prev @ par["P"]) * B[t]
        prev = v / v.sum()
        out[t] = prev
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", default="fx")
    ap.add_argument("--start", type=int, default=1996)
    ap.add_argument("--step", type=int, default=63)
    ap.add_argument("--nlag", type=int, default=250)
    ap.add_argument("--k", type=int, default=25)
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--temper", type=float, default=1.0)
    ap.add_argument("--states", type=int, default=2)
    ap.add_argument("--sticky", type=float, default=250.0, help="pseudo-counts on the diagonal of the transition matrix")
    args = ap.parse_args()
    panel = load_panel()
    close, kinds = panel["close"], panel["kind"]
    use = ("equity", "fx", "commodity") if args.universe == "all" else (args.universe,)
    df = stacked_rows(close, kinds=kinds, use=use)
    start = int(close.index.get_indexer([pd.Timestamp(f"{args.start}-01-01")], method="bfill")[0])
    pri = Priors(tau_level_scale=args.scale, temper=args.temper)
    wf = walk_forward(BADLM(nlag=args.nlag, k=args.k, priors=pri), df["x"].to_numpy(), df["y"].to_numpy(),
                      df["time"].to_numpy(), groups=df["asset"].to_numpy(), start=start, step=args.step)
    rows = df.iloc[wf.row].reset_index(drop=True)
    rows["pred"], rows["pvar"] = wf.mean, wf.var
    print(f"{len(wf.refit_times)} curve refits; static IC {ic(rows['pred'], rows['y']):.4f}", flush=True)

    # the HMM is refitted every `hmm_step` days on the out-of-sample predictions available so far,
    # after a burn-in of `burn` days, and filtered forward between refits
    stats = day_stats(rows)
    days = stats.index.to_numpy()
    burn, hmm_step = 252 * 4, args.step * 2
    amp = pd.Series(1.0, index=stats.index)
    s2f = pd.Series(1.0, index=stats.index)
    par = None
    fits = []
    t0 = days[0] + burn
    while t0 < days[-1]:
        train = stats[stats.index <= t0]
        par = fit_hmm(train, K=args.states, sticky=args.sticky)
        order = np.argsort(par["a"])
        fits.append({"time": int(t0), "a": par["a"][order].round(3).tolist(), "s2": par["s2"][order].round(3).tolist(),
                     "stay": np.diag(par["P"])[order].round(4).tolist()})
        block = stats[(stats.index > t0) & (stats.index <= t0 + hmm_step)]
        if len(block):
            # start from the filtered state at the end of training
            start_prob = filter_forward(train, par)[-1]
            probs = filter_forward(block, par, start_prob=start_prob)
            # the state used on day t + 1 is the one filtered through day t
            probs = np.vstack([start_prob, probs[:-1]])
            amp.loc[block.index] = probs @ par["a"]
            s2f.loc[block.index] = probs @ par["s2"]
        t0 += hmm_step
    print("HMM parameters over time (amplitudes a, variances s2, stay probabilities):")
    for f in fits[::max(1, len(fits) // 8)]:
        print(f"  {close.index[f['time']].date()}: a {f['a']} s2 {f['s2']} stay {f['stay']}")
    oos = rows["time"] > days[0] + burn
    sub = rows[oos].reset_index(drop=True)
    a_t = amp.reindex(sub["time"]).to_numpy()
    v_t = s2f.reindex(sub["time"]).to_numpy()
    results = evaluate("static curve", (sub["pred"] / sub["pvar"]).to_numpy(), sub, rules=(("aim", {}), ("gp", {"speed": 0.1})))
    results += evaluate("regime-scaled curve", (a_t * sub["pred"] / sub["pvar"]).to_numpy(), sub, rules=(("aim", {}), ("gp", {"speed": 0.1})))
    print(f"\nregime-scaled IC {ic(a_t * sub['pred'], sub['y']):.4f} vs static {ic(sub['pred'], sub['y']):.4f} "
          f"on {close.index[sub['time'].min()].year}-{close.index[sub['time'].max()].year}; "
          f"amplitude: mean {a_t.mean():.2f}, 10th/90th pct {np.quantile(a_t, 0.1):.2f}/{np.quantile(a_t, 0.9):.2f}")
    print_table(results, f"{args.universe}: latent amplitude regime, out of sample")
    os.makedirs(resdir, exist_ok=True)
    suffix = "" if args.sticky == 250.0 else f"_sticky{int(args.sticky)}"
    with open(os.path.join(resdir, f"trading_regimes_{args.universe}{suffix}.json"), "w") as f:
        json.dump({"args": vars(args), "fits": fits, "results": strip(results),
                   "ic_static": ic(sub["pred"], sub["y"]), "ic_regime": ic(a_t * sub["pred"], sub["y"])}, f, indent=1)
    print("written")

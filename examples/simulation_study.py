"""Simulation study of the paper, re-run with the collapsed NUTS sampler.

Five models are compared on simulated series, as in Section 4 of the paper:

  M1  non-adaptive first-difference penalty (``tau_model='constant'``)
  M2  ridge penalty increasing linearly with lag (Muggeo, 2008), Bayesian
  M3  adaptive first-difference penalty, the proposed model (``tau_model='rw'``)
  M4  first-difference penalty plus a ridge (Obermeier et al., 2015), Bayesian
  M5  B-splines with logarithmically spaced knots, their number chosen by AIC,
      unpenalised (Gasparrini et al., 2010)

M1 to M4 share one basis of K = ceil(2(p + 1) / 3) cubic B-splines with the
curve constrained to zero at lag p; M2 and M4 are fitted with the same collapsed
sampler as M1 and M3, with penalised-complexity priors on their scales. One
chain of 300 retained draws after 300 warm-up iterations, started at the Laplace
mode, is used for each fit; only posterior means are needed. Data:
x_t standardised AR(1) with coefficient 0.5; y_t = sum_j beta_j x_{t-j} + e_t
with AR(1) errors (coefficient 0.2, innovation sd 0.1), t = p + 1, ..., 500.

Study 1 compares the models on five lag-curve shapes with p = 50. Study 2
fits the displacement curve with p = 50, 75, 100 and 125.

Usage:
  python examples/simulation_study.py --time            # one replicate of each fit, timed
  python examples/simulation_study.py --reps 50 --reps2 20
Outputs go to paper/results (JSON), paper/tables (LaTeX) and paper/figures.
"""
import argparse
import json
import os
import sys
import time

import jax
import jax.numpy as jnp
import numpy as np
import numpyro
import numpyro.distributions as dist
from scipy.interpolate import BSpline

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from badlm import BADLM, paper_lag_curves, simulate_dlm  # noqa: E402
from badlm.collapsed import collapsed_loglik_Q, conditional_coefficients_Q  # noqa: E402
from badlm.design import build_design, build_structure, lag_matrix  # noqa: E402
from badlm.numpyro_model import run_nuts  # noqa: E402
from badlm.priors import Priors  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.join(here, "..")
MCMC = dict(samples=300, warmup=300, chains=1)
SIM = dict(n=500, standardize_x=True, ar_x=0.5, ar_eps=0.2, sd_eps=0.1)
PRI = Priors()


# ------------------------------------------------------------------ M2, M4
def ridge_model(stats, struct, q, kind):
    """Collapsed Gaussian model with a fixed-shape penalty and one or two scales."""
    K_r = struct.K_r
    Z = jnp.asarray(struct.Z)
    sd = numpyro.sample("sd", dist.Exponential(PRI.tau_level_rate))
    lam = sd**-2
    if kind == "M2":
        w = jnp.arange(1, struct.K + 1, dtype=float)
        Q_r = lam * (Z.T * w) @ Z
    else:
        sd_r = numpyro.sample("sd_r", dist.Exponential(PRI.tau_level_rate))
        D = jnp.asarray(np.diff(np.eye(struct.K), axis=0))
        Q_r = Z.T @ (lam * D.T @ D + sd_r**-2 * jnp.eye(struct.K)) @ Z
    numpyro.deterministic("Q_r", Q_r)
    logdet = 2.0 * jnp.sum(jnp.log(jnp.diag(jnp.linalg.cholesky(Q_r))))
    sigma2 = numpyro.sample("sigma2", dist.InverseGamma(PRI.sigma2_a, PRI.sigma2_b))
    numpyro.factor("loglik", collapsed_loglik_Q(Q_r, logdet, sigma2, 0.0, stats, q, False))


def fit_ridge(x, y, struct, kind, seed):
    data, std = build_design(x, y, struct)
    draws, extra, _ = run_nuts(ridge_model, (data.stats, struct, data.q, kind), seed=seed, **MCMC)
    Q = jnp.asarray(draws["Q_r"])
    s2 = jnp.asarray(draws["sigma2"])
    z = jax.random.normal(jax.random.PRNGKey(seed + 1), s2.shape + (struct.K_r + data.q,))
    f = jax.jit(jax.vmap(jax.vmap(lambda Qr, ss, zz: conditional_coefficients_Q(Qr, ss, 0.0, data.stats, data.q, zz))))
    theta, ed = f(Q, s2, z)
    theta = np.asarray(theta)
    curve = theta[..., :struct.K_r] @ struct.B_r.T * (std.y_sd / std.x_sd)
    jax.clear_caches()
    return curve.reshape(-1, struct.p + 1).mean(axis=0), float(np.asarray(ed).mean() - data.q)


# ------------------------------------------------------------------ M5
def logknot_basis(p, n_int):
    lags = np.arange(p + 1, dtype=float)
    interior = np.exp(np.linspace(np.log(1.0), np.log(p + 1.0), n_int + 2))[1:-1] - 1.0
    t = np.r_[np.zeros(4), interior, np.full(4, float(p))]
    return np.asarray(BSpline(t, np.eye(n_int + 4), 3, extrapolate=True)(lags))


def fit_logknots(x, y, p, max_int=12):
    L = lag_matrix(x, p)
    y = y[-L.shape[0]:] if y.size != L.shape[0] else y
    best = None
    for n_int in range(0, max_int + 1):
        B = logknot_basis(p, n_int)
        X = np.hstack([np.ones((L.shape[0], 1)), L @ B])
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        rss = np.sum((y - X @ coef) ** 2)
        aic = L.shape[0] * np.log(rss / L.shape[0]) + 2 * X.shape[1]
        if best is None or aic < best[0]:
            best = (aic, B @ coef[1:], B.shape[1])
    return best[1], float(best[2])


# ------------------------------------------------------------------ M1, M3
def fit_badlm(x, y, p, K, tau_model, seed):
    m = BADLM(nlag=p, k=K, tau_model=tau_model).fit(x, y, seed=seed, **MCMC)
    return m.summary_.mean, m.effective_df()["mean"]


def fit_all(x, y, p, K, seed, models=("M1", "M2", "M3", "M4", "M5")):
    y = y[y.size - (x.size - p):]                 # rows t = p + 1, ..., n
    struct = build_structure(p, K, end_constraint="zero")
    out = {}
    for name in models:
        t0 = time.time()
        if name == "M1":
            est, ed = fit_badlm(x, y, p, K, "constant", seed)
        elif name == "M3":
            est, ed = fit_badlm(x, y, p, K, "rw", seed)
        elif name in ("M2", "M4"):
            est, ed = fit_ridge(x, y, struct, name, seed)
        else:
            est, ed = fit_logknots(x, y, p)
        out[name] = (np.asarray(est), ed, time.time() - t0)
    return out


def metrics(est, truth):
    return {"rmse": float(np.sqrt(np.mean((est - truth) ** 2))), "bias2": float(np.mean(est - truth) ** 2)}


def k_for(p):
    return int(np.ceil(2.0 * (p + 1) / 3.0))


# ------------------------------------------------------------------ driver
def load_partial(out, keys):
    """Resume from a previous run: results saved after every replicate."""
    results = {k: {m: [] for m in ("M1", "M2", "M3", "M4", "M5")} for k in keys}
    if os.path.exists(out):
        saved = json.load(open(out))
        for k in keys:
            for m in results[k]:
                results[k][m] = saved.get(k, {}).get(m, [])
    return results


def run_study1(reps, out):
    curves = paper_lag_curves(50)
    results = load_partial(out, list(curves))
    t0 = time.time()
    for cname, truth in curves.items():
        for r in range(reps):
            if len(results[cname]["M3"]) > r:
                continue
            x, y = simulate_dlm(truth, seed=10_000 + r, **SIM)
            fits = fit_all(x, y, 50, k_for(50), seed=r)
            for m, (est, ed, dt) in fits.items():
                results[cname][m].append({**metrics(est, truth), "ed": ed, "seconds": dt})
            if (r + 1) % 5 == 0:
                print(f"study 1 {cname} rep {r + 1}/{reps} [{time.time() - t0:.0f}s]", flush=True)
            json.dump(results, open(out, "w"))
    return results


def run_study2(reps, out):
    truth = paper_lag_curves(50)["displacement"]
    ps = (50, 75, 100, 125)
    results = load_partial(out, [str(p) for p in ps])
    t0 = time.time()
    for r in range(reps):
        if len(results[str(ps[-1])]["M3"]) > r:
            continue
        x, y = simulate_dlm(truth, seed=20_000 + r, **SIM)
        for p in ps:
            if len(results[str(p)]["M3"]) > r:
                continue
            full = np.r_[truth, np.zeros(p - 50)]
            fits = fit_all(x, y, p, k_for(p), seed=r)
            for m, (est, ed, dt) in fits.items():
                results[str(p)][m].append({**metrics(est, full), "ed": ed, "seconds": dt})
        print(f"study 2 rep {r + 1}/{reps} [{time.time() - t0:.0f}s]", flush=True)
        json.dump(results, open(out, "w"))
    return results


def summarise(results):
    return {k: {m: {s: float(np.mean([r[s] for r in v])) for s in ("rmse", "bias2", "ed", "seconds")}
                for m, v in d.items() if v} for k, d in results.items()}


def write_outputs(res1, res2):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(os.path.join(root, "paper", "tables"), exist_ok=True)
    s1, s2 = summarise(res1), summarise(res2)
    labels = {"delayed_peak": "Delayed peak", "decay_curve": "Decay curve", "flat_response": "Flat response",
              "displacement": "Displacement", "sharp_peak": "Sharp peak"}
    lines = ["\\begin{tabular}{llrrrrr}", "\\hline",
             "Lag curve shape & Summary & $\\M{1}$ & $\\M{2}$ & $\\M{3}$ & $\\M{4}$ & $\\M{5}$ \\\\", "\\hline"]
    for c, lab in labels.items():
        if c not in s1:
            continue
        for i, (stat, name) in enumerate((("rmse", "RMSE"), ("bias2", "Bias$^2$"), ("ed", "ED"))):
            scale = {"rmse": 1e3, "bias2": 1e6, "ed": 1.0}[stat]
            cells = " & ".join(f"{scale * s1[c][m][stat]:.2f}" for m in ("M1", "M2", "M3", "M4", "M5"))
            first = f"\\scen{{{lab}}}" if i == 1 else ""
            lines.append(f"{first} & {name} & {cells} \\\\")
        lines.append("\\hline")
    lines.append("\\end{tabular}")
    open(os.path.join(root, "paper", "tables", "simulation.tex"), "w").write("\n".join(lines) + "\n")

    # Figure 1: the five lag curves
    curves = paper_lag_curves(50)
    fig, ax = plt.subplots(figsize=(7, 3.6))
    for (c, truth), col in zip(curves.items(), ("k", "C3", "C2", "C0", "C1")):
        ax.plot(np.arange(51), truth, color=col, lw=2, label=labels[c])
    ax.set_xlabel("lag")
    ax.set_ylabel("lag function")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(root, "paper", "figures", "fig1_lag_curves.pdf"))
    fig.savefig(os.path.join(root, "figures", "simulation_curves.png"), dpi=130)

    # Figure 2: ED against RMSE for each model and p
    fig, ax = plt.subplots(figsize=(7, 5))
    cols = {"M1": "C2", "M2": "C0", "M3": "k", "M4": "C3", "M5": "C1"}
    for m in ("M1", "M2", "M3", "M4", "M5"):
        pts = [(1e3 * s2[str(p)][m]["rmse"], s2[str(p)][m]["ed"]) for p in (50, 75, 100, 125) if str(p) in s2]
        if not pts:
            continue
        ax.plot([a for a, _ in pts], [b for _, b in pts], color=cols[m], lw=2, zorder=1)
        for (a, b), size in zip(pts, (20, 45, 80, 125)):
            ax.scatter(a, b, s=size, color=cols[m], zorder=2)
        ax.annotate(f"$\\mathcal{{M}}_{m[1]}$", pts[-1], textcoords="offset points", xytext=(8, 4), fontsize=12)
    for size, p in zip((20, 45, 80, 125), (50, 75, 100, 125)):
        ax.scatter([], [], s=size, color="grey", label=f"p = {p}")
    ax.legend(title="Maximum lag", frameon=False, loc="upper left")
    ax.set_xlabel("RMSE (x 1000)")
    ax.set_ylabel("effective degrees of freedom (ED)")
    fig.tight_layout()
    fig.savefig(os.path.join(root, "paper", "figures", "fig2_ed_rmse.pdf"))
    fig.savefig(os.path.join(root, "figures", "simulation_ed_rmse.png"), dpi=130)
    json.dump({"study1": s1, "study2": s2}, open(os.path.join(root, "paper", "results", "summary.json"), "w"), indent=1)
    print(json.dumps({"study1": s1, "study2": s2}, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--reps2", type=int, default=20)
    ap.add_argument("--time", action="store_true")
    ap.add_argument("--outputs-only", action="store_true", help="rebuild tables and figures from saved results")
    args = ap.parse_args()
    resdir = os.path.join(root, "paper", "results")
    os.makedirs(resdir, exist_ok=True)
    if args.time:
        for p in (50, 125):
            truth = paper_lag_curves(50)["delayed_peak"]
            x, y = simulate_dlm(truth, seed=1, **SIM)
            L = lag_matrix(x, 50)
            sig = L @ truth
            print(f"p={p}: signal sd {sig.std():.3f}, noise sd {(y - sig).std():.3f}, R2 {np.var(sig) / np.var(y):.2f}")
            fits = fit_all(x, y, p, k_for(p), seed=1)
            full = np.r_[truth, np.zeros(p - 50)]
            for m, (est, ed, dt) in fits.items():
                print(f"  {m}: {dt:6.1f}s  rmse {1e3*metrics(est, full)['rmse']:.2f}  ED {ed:.1f}", flush=True)
    elif args.outputs_only:
        write_outputs(json.load(open(os.path.join(resdir, "study1.json"))), json.load(open(os.path.join(resdir, "study2.json"))))
    else:
        res1 = run_study1(args.reps, os.path.join(resdir, "study1.json"))
        res2 = run_study2(args.reps2, os.path.join(resdir, "study2.json"))
        write_outputs(res1, res2)

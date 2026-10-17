"""Gibbs sampler for the Gaussian model with slice updates of the log precisions.

Kept as an independent cross-check of the collapsed NUTS sampler and as the
sampler validated by the "getting it right" test. Supports ``tau_model='rw'``,
any end constraint, and unpenalised columns; AR(1) errors are not supported.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular

from .priors import Priors, pc_logprior_logzeta, rw_logprior
from .slice import slice_sample


class GibbsSampler:
    def __init__(self, data, struct, pri: Priors, rng, interweave=True):
        self.struct = struct
        self.pri = pri
        self.rng = rng
        self.q = data.q
        self.K_r = struct.K_r
        self.N = data.N
        self.X = data.X
        self.S = data.stats.S_a
        self.set_y(data.y - data.offset)
        self.interweave = interweave and struct.m <= 1
        kt = struct.K_tau
        # initial state: heavily smoothed curve
        self.tau = np.full(kt, pri.tau_level_mean)
        self.zeta = 0.5
        self.sigma2 = 1.0
        self.theta = np.zeros(self.K_r + self.q)
        self._cache_prior()

    # ------------------------------------------------------------------ data
    def set_y(self, y):
        self.y = np.asarray(y, float)
        self.s = self.X.T @ self.y
        self.yy = float(self.y @ self.y)

    # ---------------------------------------------------------------- helpers
    def _cache_prior(self):
        self.lam = np.exp(self.tau)

    def _rss(self, theta):
        return self.yy - 2.0 * self.s @ theta + theta @ self.S @ theta

    def loglik(self, theta, sigma2=None):
        sigma2 = self.sigma2 if sigma2 is None else sigma2
        return -0.5 * self.N * np.log(2 * np.pi * sigma2) - 0.5 * self._rss(theta) / sigma2

    def _Q_full(self, tau=None):
        tau = self.tau if tau is None else tau
        Q = np.zeros((self.K_r + self.q, self.K_r + self.q))
        Q[: self.K_r, : self.K_r] = self.struct.lam_matrix(tau)
        return Q

    # ----------------------------------------------------------------- steps
    def update_theta(self):
        A = self._Q_full() + self.S / self.sigma2
        c, low = cho_factor(A, lower=True)
        mu = cho_solve((c, low), self.s / self.sigma2)
        z = self.rng.normal(size=mu.size)
        self.theta = mu + solve_triangular(c, z, lower=True, trans="T")

    def update_sigma2(self):
        a = self.pri.sigma2_a + 0.5 * self.N
        b = self.pri.sigma2_b + 0.5 * self._rss(self.theta)
        self.sigma2 = 1.0 / self.rng.gamma(a, 1.0 / b)

    def update_tau_centred(self):
        st = self.struct
        b_r = self.theta[: self.K_r]
        d2 = (st.D_r @ b_r) ** 2
        tau = self.tau.copy()
        for k in range(tau.size):
            def logf(t, k=k):
                tt = tau.copy()
                tt[k] = t
                ld = (0.5 * t) if st.m <= 1 else 0.5 * st.logdet_prior(tt)
                return ld - 0.5 * np.exp(t) * d2[k] + rw_logprior(tt, self.zeta, self.pri)
            tau[k], _ = slice_sample(logf, tau[k], self.rng)
        self.tau = tau
        self._cache_prior()

    def _theta_from_whitened(self, tau, w, level, gamma):
        st = self.struct
        d = np.exp(-0.5 * tau) * w
        b = st.coefficients_from_increments(d, level)
        b_r = st.Z.T @ b
        return np.concatenate([b_r, gamma])

    def update_tau_interweaved(self):
        """Non-centred update: hold the whitened increments ``w = exp(tau/2) * d``
        fixed and update each ``tau_k`` against the likelihood."""
        st = self.struct
        b = st.Z @ self.theta[: self.K_r]
        gamma = self.theta[self.K_r:]
        level = b[0] if st.m == 0 else None
        d = np.diff(b)
        w = np.exp(0.5 * self.tau) * d
        tau = self.tau.copy()
        for k in range(tau.size):
            def logf(t, k=k):
                tt = tau.copy()
                tt[k] = t
                theta = self._theta_from_whitened(tt, w, level, gamma)
                return self.loglik(theta) + rw_logprior(tt, self.zeta, self.pri)
            tau[k], _ = slice_sample(logf, tau[k], self.rng)
        self.tau = tau
        self._cache_prior()
        self.theta = self._theta_from_whitened(tau, w, level, gamma)

    def update_tau_level(self):
        """Shift all log precisions together: a cheap move along the direction
        in which one-at-a-time updates mix worst."""
        st = self.struct
        b_r = self.theta[: self.K_r]
        d2 = (st.D_r @ b_r) ** 2
        tau0 = self.tau.copy()

        def logf(delta):
            tt = tau0 + delta
            ld = 0.5 * st.logdet_prior(tt)
            return ld - 0.5 * np.sum(np.exp(tt) * d2) + rw_logprior(tt, self.zeta, self.pri)

        delta, _ = slice_sample(logf, 0.0, self.rng)
        self.tau = tau0 + delta
        self._cache_prior()

    def update_zeta(self):
        def logf(lz):
            return rw_logprior(self.tau, np.exp(lz), self.pri) + pc_logprior_logzeta(lz, self.pri)
        lz, _ = slice_sample(logf, np.log(self.zeta), self.rng)
        self.zeta = float(np.exp(lz))

    def effective_df(self):
        A = self.S + self.sigma2 * self._Q_full()
        return float(np.trace(np.linalg.solve(A, self.S)))

    def step(self):
        self.update_theta()
        self.update_sigma2()
        self.update_tau_centred()
        self.update_tau_level()
        if self.interweave:
            self.update_tau_interweaved()
        self.update_zeta()

    def state(self):
        return {"theta": self.theta.copy(), "sigma2": self.sigma2, "tau": self.tau.copy(),
                "zeta": self.zeta}

    def run(self, samples, warmup, thin=1, ed=True):
        n_save = (samples) // thin
        out = {"theta": np.empty((n_save, self.theta.size)), "sigma2": np.empty(n_save),
               "tau": np.empty((n_save, self.tau.size)), "zeta": np.empty(n_save),
               "ed": np.empty(n_save)}
        for i in range(warmup + samples):
            self.step()
            j = i - warmup
            if j >= 0 and j % thin == 0 and j // thin < n_save:
                k = j // thin
                out["theta"][k] = self.theta
                out["sigma2"][k] = self.sigma2
                out["tau"][k] = self.tau
                out["zeta"][k] = self.zeta
                out["ed"][k] = self.effective_df() if ed else np.nan
        return out


def run_gibbs(data, struct, pri, *, samples=2000, warmup=1000, chains=2, seed=0, thin=1,
              interweave=True):
    """Run independent chains and stack the draws as ``(chains, draws, ...)``."""
    runs = []
    for c in range(chains):
        rng = np.random.default_rng([seed, c])
        sampler = GibbsSampler(data, struct, pri, rng, interweave=interweave)
        runs.append(sampler.run(samples, warmup, thin=thin))
    return {k: np.stack([r[k] for r in runs]) for k in runs[0]}

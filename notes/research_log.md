# Research log: lag curves as trading signals

Short record of what was tried, what it showed and why it was dropped, so
that the repository only carries the ideas that survived. Numbers are
out of sample (expanding window, quarterly refits, 1996-2026 unless stated),
Sharpe ratios annualised, "slow" is the Garleanu-Pedersen rule moving a tenth
of the way to the aim each day, costs proportional to notional traded.

## Kept (see README, "Trading experiments")

- Pooled cross-asset curve with calibrated evidence (25 basis functions,
  `tau_level_scale=0.001`, `temper=10`): currencies 0.70 gross / 0.64 net of
  5 bp, the same as one-parameter time-series momentum. A baseline, not an edge.
- Volatility-modulated surface (curve varies with the asset's own volatility
  state): currencies 0.88 / 0.78, all 32 series 0.84 / 0.77. Beats the
  two-horizon hand rule (0.81 / 0.72). The one extension that paid.

## Dropped

- **Flexible basis with default priors** (84 functions, scale 1, no temper):
  in-sample IC 0.050, out of sample 0.014, about 38 effective degrees of
  freedom. The penalised-complexity priors are uninformative at coefficient
  sizes of a few thousandths and correlated panel rows overstate the evidence.
  Fixed by the prior scale and the tempered likelihood; a coarse basis of 10
  functions gives the same out-of-sample quality (IC 0.0215). Constant weights
  cannot temper: the free error variance absorbs them exactly.
- **Log-spaced knots** (`lag_scale='log'`): no better than linear knots once
  the evidence is calibrated (IC 0.016-0.018 against 0.019). Kept as an option
  because it is cheap, not used in the experiments.
- **Gains and losses as two exposures**: losses are followed by continued
  losses (short-lag coefficient 0.005) and gains barely by anything (0.001),
  but out of sample the gain is within noise (IC 0.0212 against 0.0204, net
  0.66 against 0.64). Interesting shape, no edge.
- **Return x relative volume** (equities only): the volume curve is zero, IC
  stays slightly negative. Eleven US names at a daily horizon show nothing,
  with or without volume; momentum rules are equally weak there (0.3).
- **Commodities** (two spot oil series and gold): 0.15, nothing to pool.
- **Time-drifting surface** (modifier = calendar time): tracks the fading of
  the momentum plateau from 0.006 in 1999 to below zero now, but prediction
  does not improve (0.68 / 0.63): the drift is only known after the fact.
- **Latent amplitude regime** (two-state hidden Markov model on the sign and
  size of the pooled signal, state shared across assets, EM refit yearly,
  filtered state scaling the position): the states flip almost daily (stay
  probabilities 0.6 with 250 pseudo-counts, still only 0.92-0.97 with 5000),
  the amplitudes drift to about -7 and +8, the filtered amplitude is a noisy
  verdict on whether the signal worked yesterday, IC falls (0.015 / 0.012
  against 0.016 on 1999-2026) and turnover rises to about 200-270 a year, so
  it loses after costs. Observed-state modulation beats latent regimes.
- **Design-effect weights** (row weight 1 / (1 + (n_t - 1) rho)): no effect,
  for the same reason constant weights have none.

## Data notes

- Lean's QQQ file is empty from 2004 to March 2011 (forward-filled zeros
  collapsed its volatility and made it dominate the portfolio); it starts
  after the gap. GOOG before April 2014 is the pre-split class and is dropped.
  WM has a 112% annual volatility from ticker reuse and is dropped.
- Henry Hub spot is stale on a third of days and moves over 100% on others:
  excluded. The emerging-market currencies are managed before 1995.
- The frictionlessdata exchange-rate package quotes every currency per
  dollar, including AUD, EUR, GBP and NZD.

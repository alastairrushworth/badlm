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

## Dropped, second batch (currencies, same set-up)

- **Exposure-level surface** `beta(j, x_{t-j})`, the cross-basis of
  Gasparrini et al. (`modifier='exposure'`, kept in the package and tested):
  the response is close to linear in the size of the past move, slightly
  larger after losses than after gains, IC 0.020 against 0.020 for the linear
  curve, 0.68 / 0.62 net. No non-linearity in the exposure worth having.
- **Momentum-crash state** (exposures split by the sign of the past-year
  trend, each with the volatility surface): 0.79 / 0.73 against 0.88 / 0.78
  for the plain volatility surface. Splitting halves the data behind each
  surface and the sign of the trend adds nothing the state did not have.
- **Dollar-neutral (cross-sectional) returns**, x and y demeaned across the
  currencies each day: static 0.51 / 0.45, surface 0.48 / 0.34. The common
  dollar component carries much of the predictability.
- **Market-wide volatility state** (cross-sectional mean of the states):
  0.60 / 0.53; own minus market: 0.77 / 0.70. The asset's own state is the
  right modifier, not a global stress indicator.

## Third batch: variants of the volatility surface (currencies)

- **Modifier resolution and definition**: 4, 5 or 7 state splines give
  0.71, 0.78 and 0.77 net of 5 bp; the absolute log volatility level
  instead of the level relative to the asset's expanding median gives 0.81
  (IC 0.029 against 0.027); a fast 5-day-over-60-day realised-volatility
  ratio gives 0.61. The state that matters is the slow level of volatility,
  not a recent burst. The relative definition is kept as the default because
  it is scale-free across assets; the differences are within noise.
- **Trading rule on the surface**: moving a twentieth of the way to the aim
  each day (gp 0.05) is best net of costs, 0.81 at 5 bp and 0.75 at 10 bp
  with a break-even cost of 73 bp; the full-aim rule earns 0.90 gross but only
  0.44 net at 5 bp. A no-trade band of twice the mean aim gives 0.76 at 5 bp.
- **Mixing with the two-horizon hand rule**: a half-and-half mix gives
  0.79 net against 0.78 for the surface alone, so the surface subsumes it.
- **Weekly target** (sum of the next five scaled returns): IC 0.035 on the
  weekly target, 0.67 / 0.58, the same as the daily target after costs. With
  `ar1=True` the walk-forward predictions came out strongly negative (IC
  -0.045) because the AR(1) model is right and the prediction is incomplete: with
  overlapping targets the whitened regression finds phi near one and a
  coefficient of about -phi at lag one (the target of day t - 1 contains the
  exposure of day t), and `predict`/`walk_forward` return only the regression
  mean `X beta`, without the carried-over residual `phi * e_{t-1}` that the
  AR(1) model needs for a forecast. Documented in the package; do not use
  `ar1=True` for prediction until that term is added.

## Data notes

- Lean's QQQ file is empty from 2004 to March 2011 (forward-filled zeros
  collapsed its volatility and made it dominate the portfolio); it starts
  after the gap. GOOG before April 2014 is the pre-split class and is dropped.
  WM has a 112% annual volatility from ticker reuse and is dropped.
- Henry Hub spot is stale on a third of days and moves over 100% on others:
  excluded. The emerging-market currencies are managed before 1995.
- The frictionlessdata exchange-rate package quotes every currency per
  dollar, including AUD, EUR, GBP and NZD.

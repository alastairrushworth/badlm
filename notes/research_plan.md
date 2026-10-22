# Research plan: from lag curves to robust multi-asset trading

*Status after the first pass (see the README section "The multivariate
plan: first pass" and `notes/research_log.md`): pillars A to D were built
and run on the wide panel; none beat the asset-level volatility surface,
and that surface earned nothing in the held-out years. The plan's
principles stand; its next iteration needs the data listed in pillar F
before any further modelling.*

## Where we stand

Two sessions of experiments established three facts. A single adaptive lag
curve of returns on their own past, pooled across assets, rediscovers
time-series momentum and nothing else, and that signal has faded since 2009.
The adaptive smoother overfits returns unless the evidence of a dependent,
heavy-tailed panel is tempered and the prior is set at the scale of return
predictability. And the only extension that beat the hand-made rules out of
sample gave the curve a second dimension along an observed state, the
asset's own volatility: the shape of predictability changes with conditions,
and that structure is what a flexible-but-shrunk Bayesian model is good at
finding.

The univariate question "what does an asset's past say about its future" is
close to exhausted. The plan below moves the method to multivariate
questions, where the statistical machinery (adaptive smoothing over lag,
shrinkage, calibrated evidence, honest evaluation) is reused but the model
is rebuilt around them.

## Principles

1. **Borrow strength deliberately.** Every parameter that is specific to one
   asset must be shrunk towards something shared: a group, a factor, a
   common shape. Thirty assets with thirty free curves is thirty univariate
   problems; thirty assets with one shape and thirty loadings is one problem
   with thirty times the data.
2. **Model the panel as a panel.** Cross-sectional correlation is not a
   nuisance to be tempered away for ever; it is the factor structure, and
   factors are where predictability concentrates (dollar-neutral returns
   were less predictable than raw ones).
3. **Prefer durable targets.** Return predictability is small and decays as
   it is traded. Volatility, correlation and the state-dependence of returns
   are large and have not decayed. Durability comes from targets that are
   structural rather than from cleverer return forecasts.
4. **Estimate for the decision, evaluate for the decision.** The quantity we
   care about is portfolio utility after costs, not squared error on a
   scaled return; where the two disagree the estimator should follow the
   first.
5. **Every idea gets a kill criterion before it is run**, a held-out period
   that is never touched during development, and its result logged either
   way.

## Pillar A: the lag-curve matrix (lead-lag across assets)

Replace the univariate curve beta(j) with a matrix of curves beta_ij(j):
asset i's next return on asset j's past. This is a distributed-lag vector
autoregression. Unrestricted it has N^2 curves and is hopeless; the point is
the structure that makes it estimable, which is also what makes it
multivariate thinking rather than N^2 technical rules.

- **Reduced rank over pairs, adaptive over lag.** beta_ij(l) =
  sum_r u_ir v_jr g_r(l): a few shared lag shapes g_r with adaptive
  smoothing (the existing machinery), and sparse loadings u, v (horseshoe)
  saying which assets lead and which follow. The Minnesota prior of Bayesian
  VARs is a lag-curve prior with precision rising in lag; this is its
  adaptive, reduced-rank version. Inference: g_r integrated out as now, u, v
  and the precisions sampled; the design for fixed u, v is a stack of
  weighted lag matrices.
- **Known lead-lag channels to test first**, because economics says they
  are durable: commodity prices leading commodity-exporter currencies
  (Chen, Rogoff and Rossi); large and liquid leading small and illiquid
  (Lo and MacKinlay; Hou's industry lead-lag); US close leading Asian open
  through time zones; the dollar factor leading everything priced in
  dollars; volatility (VIX) leading risk assets. Each is a row or column of
  the matrix with a hypothesised sign.
- **Kill criterion.** The reduced-rank matrix must beat the diagonal (own
  past only) out of sample on the same rows, net of costs, in two of three
  sub-periods. If the off-diagonal loadings are all shrunk to zero, the
  answer is "no lead-lag at a daily horizon in this universe" and the pillar
  closes.

## Pillar B: factors first, lags second

Decompose the panel into a few common factors and idiosyncratic residuals,
and put the lag curves on the factors.

- **Factor momentum.** Factor returns are more autocorrelated than asset
  returns (Gupta and Kelly; Ehsani and Linnainmaa), and a factor is by
  construction a basket, so trading it borrows strength from every asset
  that loads on it. The model is a small system: each factor's next return
  on its own past and the other factors' past, with adaptive curves and the
  volatility-state surface that already works. Positions in assets follow
  from loadings. Factors: statistical (first few principal components of
  scaled returns, estimated on a rolling window) and economic (dollar,
  carry, global equity, energy, metals).
- **State-dependent surfaces with a vector state.** The volatility surface
  used one state. Learn a single index w'z_t over several states (own
  volatility, dispersion across assets, average correlation, the factor's
  own volatility, carry) and let the curve vary with the index: a
  varying-coefficient model with a learned direction. One extra vector of
  hyperparameters, so the collapsed sampler still applies, and the learned
  direction is itself a finding.
- **Hierarchical amplitudes.** Asset-specific curves as a shared shape times
  an asset loading, loadings shrunk towards one within a group and groups
  towards a global value. Three levels of pooling in one fit; the loadings
  are hyperparameters, the shape is integrated out.
- **Kill criterion.** Factor-level prediction must give a higher net Sharpe
  ratio than the asset-level surface on the same universe and period; the
  learned state index must improve on own volatility alone in the held-out
  years.

## Pillar C: durable targets, volatility and correlation across assets

The most predictable quantities in the panel are second moments, and the
package's best-fitting example was already a volatility lag curve.

- **A lag-curve VAR on log realised volatility across assets** is an
  adaptive version of the Diebold-Yilmaz spillover model: who transmits
  volatility to whom, over what horizon, with the reduced-rank structure of
  pillar A. Volatility clustering has not faded in a century of data.
- **Monetisation.** First, volatility targeting of everything in pillars A
  and B with the multi-asset forecast rather than a per-asset EWMA; the
  literature (Moreira and Muir) and our own vol-state result say the sizing
  is where much of the value is. Second, with options or VIX-futures data,
  relative-value variance trades where the forecast disagrees with the
  implied level. Third, forecast correlation for the portfolio construction
  of pillar B.
- **Kill criterion.** The multi-asset volatility forecast must beat the
  per-asset EWMA and HAR on log score out of sample, and the strategy sized
  with it must have a higher net Sharpe ratio than the same strategy sized
  with the EWMA.

## Pillar D: estimate the policy, not the forecast

The lag curve is a linear policy: today's position is a weighted sum of
past returns. Brandt, Santa-Clara and Valkanov's parametric portfolio
policies estimate such weights by maximising realised portfolio utility
pooled across assets, which is a different objective from squared forecast
error and a better one for trading.

- **Adaptive-smooth policy functions.** Keep the B-spline representation
  and the adaptive smoothness prior over lag, but replace the Gaussian
  likelihood by a pseudo-likelihood proportional to exp(portfolio utility),
  with the utility computed on the whole cross-section (so covariance
  between assets is accounted for, not assumed away) and with a turnover
  penalty inside the objective, so the slow-trading rule is estimated
  rather than imposed. This is a Bayesian M-estimator; the tempering
  question becomes the choice of the utility's scale, which has a decision-
  theoretic meaning.
- **Why this is multivariate.** A forecast-based rule sizes each asset on
  its own prediction; a policy maximising portfolio utility learns to down-
  weight signals that are redundant across correlated assets and to use
  cross-asset signals only where they diversify.
- **Kill criterion.** The policy estimate must beat the plug-in rule from
  the forecast model net of costs on held-out years, and its curve must be
  recognisably stable across sub-periods; a policy that changes shape every
  refit is noise.

## Pillar E: durability and honesty

- **Held-out years** (2019 onwards, say) frozen before any further tuning,
  used once at the end of each pillar.
- **Sub-period consistency** as a first-class statistic: a signal that works
  only in 1996-2005 is not durable; report the three-decade split for every
  result and prefer the smaller edge that holds in all three.
- **Combinatorial purged cross-validation** and the deflated Sharpe ratio
  over the full count of configurations tried, kept in the research log.
- **Posterior predictive position sizing**: positions shrink with parameter
  uncertainty (the predictive variance already includes it), and ensembles
  over basis size and temperature replace point choices.
- **Decay monitoring**: fit the time-drifting surface for diagnosis only, as
  a running estimate of how much of each edge remains.

## Pillar F: data, the binding constraint

Everything above needs a wider and better panel than nineteen currencies
and spot oil. In rough order of value per effort:

1. Kenneth French daily industry and characteristic portfolios (1926 to
   now, dozens of series, free): the natural universe for the lag-curve
   matrix and factor momentum, with Hou's lead-lag as a known target.
2. Futures continuous contracts with roll adjustment (equity indices,
   bonds, currencies, commodities): the universe where the costs assumed
   here are real.
3. Short rates and forward points for carry, both as an exposure and as a
   state.
4. Daily realised volatility from intraday data, and VIX futures or
   variance swap levels, for pillar C.
5. Crypto (24-hour, many coins, strong cross-coin lead-lag, less crowded):
   a cheap testbed for pillar A at hourly resolution.

## Sequencing

1. **Data and evaluation (weeks 1-2).** Assemble the French and futures
   panels, freeze the held-out years, extend `badlm.backtest` with
   purged cross-validation, sub-period reporting and ensemble sizing.
2. **Pillar B first (weeks 3-5).** Factor construction, factor-level
   curves with the volatility surface, the learned state index. Lowest
   cost, builds directly on what works.
3. **Pillar A (weeks 6-9).** Reduced-rank lag-curve matrix; test the
   economic channels one at a time before the full matrix.
4. **Pillar C (weeks 10-12).** Volatility VAR; re-size pillars A and B with
   it.
5. **Pillar D (weeks 13-16).** Policy estimation on whichever of A and B
   survived.
6. **Held-out evaluation and write-up.**

## What changes in the package

- A design layer for several *series* per observation (asset-by-asset lag
  matrices with shared shapes and asset loadings), not just several
  exposures of one series: the reduced-rank panel component.
- Loadings and state-index directions as hyperparameters of the collapsed
  model, with horseshoe and hierarchical priors.
- A pseudo-likelihood option (utility-based) alongside the Gaussian one.
- Log-score and portfolio-utility evaluation, purged cross-validation and
  sub-period reports in `badlm.backtest`.
- The AR(1) carry-over in `predict`, so error dynamics can be used for
  forecasting.

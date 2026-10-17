# Changes relative to arXiv:1801.06670v1

`badlm.tex` is a reconstruction of the preprint from its PDF (the original
LaTeX source was not available), with the edits listed below. Figure 1 is the
raster image embedded in the PDF; Figure 2 is the vector plot cropped from
page 12 of the PDF.

## Changes to the model description (Section 2)

These bring the paper into line with the `badlm` package.

- The end-of-curve penalty `rho`, added to the last diagonal element of `Q`
  in v1, is replaced by the pseudo-observation the package actually uses: the
  basis row at lag `p` with a zero response, which shrinks `beta_p` towards
  zero at the scale of the error variance. The displayed `Q` no longer carries
  `rho`, and the text notes that the prior is intrinsic (rank `K - 1`).
- A paragraph on the MCMC scheme now states the closed form of the
  generalised determinant, `|Q|_+ = K * prod(lambda_k)`, which the sampler uses
  in the Metropolis update of the log precisions. The earlier implementation
  evaluated `determinant(Q)` numerically on a singular matrix.
- `tau = log(lambda)` is now defined before the random walk prior on `tau` is
  stated; v1 used `tau` without defining it.
- The prior hyperparameters are stated as `sigma^2 ~ IG(0.001, 0.001)` and
  `zeta^2 ~ IG(1, 1/2)`, the values used in the package. v1 stated
  `IG(1, 1/2)` for both. **Please check which values produced the simulation
  results in Section 3 and adjust the sentence if necessary.**
- The B-spline basis is described as uniformly spaced over the lag range
  `0, ..., p` rather than "over the range of x".
- The sentence crediting a mixture of R and C++ (Rcpp) code is replaced by a
  description of the inference scheme used by the Python package `badlm`:
  the spline coefficients are integrated out analytically and the
  low-dimensional marginal posterior of the log precisions and error variance
  is sampled with the No-U-Turn Sampler. The paragraph also notes the Gibbs
  sampler, Laplace approximation, exact end constraint, AR(1) errors and count
  responses that the software provides.
- A sentence notes that the end constraint can be imposed exactly through a
  null-space restriction, the software's default.

## Corrections

- Conditional distribution of `b_k`: the last case read `k = p` and the first
  case had `lambda_i`; now `k = K` and `lambda_1`.
- Introduction: the unresolved citation placeholders
  `(zanobetti, peng, welty, Gasparrini et al. (2010))` are replaced by
  Zanobetti et al. (2000), Welty et al. (2009) and Gasparrini et al. (2010).
  No standalone Peng reference appeared in the bibliography, so none is cited;
  add one if intended.
- Table 1 caption said "four different types of lag curve shape"; there are
  five. The table is now referenced as Table 1 rather than "Table 3.1".
- The null scenario is called **Flat response** throughout, matching the
  figure legend and the table (v1 also used "Null curve").
- Discussion: the sentence fragment "and that some existing approaches." is
  removed.
- Minor wording: "a time-varying explanatory variables" -> "variable";
  "various form of" -> "various forms of"; "in the prior in in Equation" ->
  "in the prior in Equation"; "could grew" -> "could grow"; "even although"
  -> "even though"; "the adaptive models M3" -> "the adaptive model M3";
  "sensitivity of different model" -> "models"; "(i, j)th element" ->
  "(i, k)th element" for `X_{i,k}`; "the number of logarithmically spaced
  knots prior that minimises AIC" -> "that minimises AIC".

## Not changed

- All simulation results, figures and the structure of Sections 3 and 4.
- The bibliography entries were reconstructed from the reference list; volume
  and page numbers for a few entries (Stan, Koyck, Nerlove) were completed
  from the published versions and should be checked.

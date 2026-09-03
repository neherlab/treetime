use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

/// Which topology-cleanup steps run during the optimize loop.
///
/// Each step is applied once per iteration unless disabled here. All default to
/// on, so an optimize run with no topology flags behaves as it did before the
/// switches existed. The steps map to the loop's cleanup passes:
///
/// - `collapse_short_branches`: contract internal edges the per-edge optimizer
///   drove to exactly zero that carry no mutations (v0's `prune_short_branches`).
/// - `merge_siblings`: group polytomy siblings that share substitutions under a
///   new internal node. Requires discrete per-edge mutations, so sparse only.
/// - `flip_parent_child`: hoist a child that carries the exact reversion of a
///   parent-edge substitution under a new node, removing one mutation per
///   reverted position. Requires discrete per-edge mutations, so sparse only.
#[derive(Copy, Clone, Debug, PartialEq, Eq, SmartDefault)]
pub struct TopologyOps {
  #[default = true]
  pub collapse_short_branches: bool,
  #[default = true]
  pub merge_siblings: bool,
  #[default = true]
  pub flip_parent_child: bool,
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
pub enum BranchLengthMode {
  Input,
  #[default]
  Marginal,
}

/// Per-edge branch length optimization method.
///
/// Controls how `run_optimize_mixed()` finds the maximum-likelihood branch
/// length for each edge. Two orthogonal axes: algorithm (Newton-Raphson
/// vs Brent's method) and parameterization ($t$, $\sqrt{t}$, $\ln(t)$).
#[derive(Copy, Clone, Debug, PartialEq, Eq, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
pub enum BranchOptMethod {
  /// Brent's method in $t$ space (derivative-free, bracket-based).
  ///
  /// Finds the maximum within a bracket derived from the grid search bounds.
  /// Convergence is independent of Hessian conditioning. Uses `argmin::BrentOpt`.
  /// Included for completeness; `brent-sqrt` dominates for convergence speed.
  Brent,

  /// Brent's method in $\sqrt{t}$ space.
  ///
  /// Matches v0 exactly (same algorithm, same parameterization). The $\sqrt{t}$
  /// reparameterization smooths the objective, giving parabolic interpolation
  /// a better fit. Default method for golden master comparison against v0.
  #[default]
  BrentSqrt,

  /// Brent's method in $\ln(t)$ space.
  ///
  /// Smoothest objective of all parameterizations, giving the best parabolic
  /// interpolation. Requires a finite lower bound in log-space.
  BrentLog,

  /// Newton-Raphson in $t$ space.
  ///
  /// Baseline Newton method matching RAxML-NG/IQ-TREE. The Poisson indel
  /// Hessian ($-k/t^2$) can dominate the substitution Hessian on short
  /// branches, causing the step-size convergence criterion to fire before
  /// the combined gradient reaches zero.
  Newton,

  /// Newton-Raphson in $\sqrt{t}$ space.
  ///
  /// Reparameterizes the optimization variable as $s = \sqrt{t}$ and applies
  /// the chain rule to transform derivatives. Reduces the indel Hessian
  /// singularity from $O(1/t^2)$ to $O(1/t)$. Residual dominance on extreme
  /// cases ($t < 0.001$, $k > 10$).
  NewtonSqrt,

  /// Newton-Raphson in $\ln(t)$ space.
  ///
  /// Eliminates the indel singularity entirely ($\ell''_{\text{indel}} = -\mu t$,
  /// bounded). Natural relative tolerance. Best conditioning of all Newton
  /// variants.
  NewtonLog,
}

/// Controls the initial branch length estimate that runs before Newton
/// optimization.
///
/// The estimate computes `#substitutions / effective_alignment_length` per
/// edge from the marginal reconstruction. When input trees already carry
/// well-calibrated branch lengths (e.g. from RAxML, IQ-TREE, or a previous
/// TreeTime run), preserving those values lets Newton converge from a
/// better starting position.
#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
pub enum InitialGuessMode {
  /// Estimate only edges with missing or invalid branch lengths, preserve
  /// valid input values. No-op when all edges have finite branch lengths.
  #[default]
  Auto,
  /// Estimate all edges, overwriting input branch lengths.
  Always,
  /// Use input branch lengths as-is. Fails if any edge has a missing or
  /// invalid branch length.
  Never,
}

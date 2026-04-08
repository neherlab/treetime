use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

/// Per-edge branch length optimization method.
///
/// Controls how `run_optimize_mixed()` finds the maximum-likelihood branch
/// length for each edge.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum, SmartDefault, Serialize, Deserialize)]
#[value(rename_all = "kebab-case")]
pub enum BranchOptMethod {
  /// Newton-Raphson in $\sqrt{t}$ space.
  ///
  /// Reparameterizes the optimization variable as $s = \sqrt{t}$ and applies
  /// the chain rule to transform derivatives:
  ///
  ///   $d\ell/ds = 2s \cdot d\ell/dt$
  ///   $d^2\ell/ds^2 = 4s^2 \cdot d^2\ell/dt^2 + 2 \cdot d\ell/dt$
  ///
  /// The reparameterization reduces the indel Hessian singularity from
  /// $O(1/t^2)$ to $O(1/t)$, improving conditioning of the combined
  /// (substitution + indel) objective on short branches with indels.
  #[default]
  NewtonSqrt,

  /// Newton-Raphson in $t$ space.
  ///
  /// The Poisson indel Hessian ($-k/t^2$) can dominate the substitution
  /// Hessian on short branches, causing the step-size convergence criterion
  /// to fire before the combined gradient reaches zero.
  Newton,

  /// Brent's method (derivative-free, bracket-based).
  ///
  /// Finds the maximum within a bracket derived from the grid search bounds.
  /// Convergence is independent of Hessian conditioning. Uses `argmin::BrentOpt`.
  Brent,
}

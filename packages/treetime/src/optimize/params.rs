use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

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
#[serde(rename_all = "kebab-case")]
pub enum BranchLengthMode {
  Input,
  #[default]
  Marginal,
}

/// Per-edge maximum-likelihood branch-length optimizer. Variants combine an
/// algorithm with a parameterization ($t$, $\sqrt{t}$, or $\ln(t)$).
#[derive(Copy, Clone, Debug, PartialEq, Eq, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum BranchOptMethod {
  /// Brent's derivative-free method in $t$ space, with a bracket from the grid
  /// bounds. Its convergence does not depend on Hessian conditioning.
  Brent,

  /// Brent's method in $\sqrt{t}$ space. This default matches v0's algorithm
  /// and parameterization.
  #[default]
  BrentSqrt,

  /// Brent's method in $\ln(t)$ space.
  ///
  /// Smoothest objective of all parameterizations, giving the best parabolic
  /// interpolation. Requires a finite lower bound in log-space.
  BrentLog,

  /// Newton-Raphson in $t$ space. On short branches, the Poisson indel Hessian
  /// can make step-size convergence precede a zero combined gradient.
  Newton,

  /// Newton-Raphson in $\sqrt{t}$ space. The chain-rule transformation reduces
  /// the indel Hessian singularity from $O(1/t^2)$ to $O(1/t)$.
  NewtonSqrt,

  /// Newton-Raphson in $\ln(t)$ space. This removes the indel singularity and
  /// gives a natural relative tolerance.
  NewtonLog,
}

/// Controls whether marginal reconstruction estimates initial branch lengths
/// from substitutions divided by effective alignment length. Preserving valid
/// input lengths can provide a better Newton starting point.
#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, Serialize, Deserialize, JsonSchema)]
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

#[derive(Copy, Debug, Clone, PartialEq, Eq)]
pub enum ExistingBranchLengths {
  Keep,
  Overwrite,
}

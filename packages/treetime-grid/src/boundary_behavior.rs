use crate::hard_approach_law::HardApproachLaw;
use crate::soft_tail_law::SoftTailLaw;
use serde::{Deserialize, Serialize};

/// Behavior of a [`GridFn`](crate::GridFn) when evaluated outside its grid support.
///
/// A bare grid function is a generic interpolant with no probabilistic meaning, so the
/// default is [`BoundaryBehavior::Error`]: a query outside the grid is a programming error
/// unless the caller has declared how the tail should behave. The two non-default variants
/// are explicit opt-ins used by finite-support distributions and by the timetree message
/// passes, which assign a per-side tail policy to each message.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub enum BoundaryBehavior {
  /// Out-of-support evaluation is an error.
  #[default]
  Error,
  /// Hard boundary: the domain terminates and probability is zero beyond the grid edge.
  /// An optional [`HardApproachLaw`] provides edge-relative interpolation between the hard boundary
  /// and the nearest grid point: a power-law term for divergent branches plus a linear term that
  /// carries a finite mode sitting on the boundary. The raw `GridFn` returns `T::zero()` outside support;
  /// the distribution layer maps this to the policy-correct zero-probability value
  /// (e.g. `+inf` under negative-log representation).
  Hard(Option<HardApproachLaw>),
  /// Return the nearest boundary value (`y[0]` to the left, `y[n-1]` to the right),
  /// i.e. a flat tail. Use when the function is genuinely uninformative beyond the edge.
  Constant,
  /// Soft boundary: the density continues past the grid edge along a log-linear
  /// [`SoftTailLaw`] (an exponential probability tail). Unlike `Constant` it has finite mass,
  /// so it does not corrupt the quantile and HPD integrals. `Linear(None)` declares a soft edge
  /// with no fitted law yet and falls back to the flat boundary value.
  Linear(Option<SoftTailLaw>),
}

impl BoundaryBehavior {
  /// Whether this boundary is *soft*: the distribution continues past the grid edge, so the
  /// grid boundary is a representation choice (where interpolation stopped), not a fact about
  /// the distribution. A soft boundary is evaluable beyond the grid via its declared tail law
  /// and therefore *extends* the evaluable domain under multiplication.
  ///
  /// The complement is a *hard* boundary: the grid edge terminates the domain (`Hard`: zero
  /// probability beyond; `Error`: out-of-support evaluation is undefined). A hard boundary
  /// *restricts* the result domain and is never silently extended.
  ///
  /// `Constant` (flat) and `Linear` (log-linear) are the soft tails. This predicate is what the
  /// multiplication rule keys off, so both extend the domain through the same rule.
  pub fn is_soft(self) -> bool {
    matches!(self, BoundaryBehavior::Constant | BoundaryBehavior::Linear(_))
  }

  /// The fitted power-law approach for a `Hard` boundary, if any.
  pub fn approach_law(self) -> Option<HardApproachLaw> {
    match self {
      BoundaryBehavior::Hard(law) => law,
      _ => None,
    }
  }

  /// The fitted log-linear tail for a `Linear` soft boundary, if any.
  pub fn soft_law(self) -> Option<SoftTailLaw> {
    match self {
      BoundaryBehavior::Linear(law) => law,
      _ => None,
    }
  }
}

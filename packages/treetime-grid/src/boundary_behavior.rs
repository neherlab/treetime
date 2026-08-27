use crate::hard_approach_law::HardApproachLaw;
use crate::soft_tail_law::SoftTailLaw;
use serde::{Deserialize, Serialize};

/// Number of edge grid points a boundary law is fit from by default.
///
/// Both [`HardApproachLaw::fit`] and [`SoftTailLaw::fit`] read the innermost/outermost points nearest
/// an edge; this is the shared count so the branch-length approach fit and the regrid refit use one
/// value rather than drifting copies.
pub const DEFAULT_TAIL_FIT_POINTS: usize = 5;

/// Behavior of a [`GridFn`](crate::GridFn) when evaluated outside its grid support.
///
/// A bare grid function is a generic interpolant with no probabilistic meaning, so the
/// default is [`BoundaryBehavior::Error`]: a query outside the grid is a programming error
/// unless the caller has declared how the tail should behave. The non-default variants are
/// explicit opt-ins used by finite-support distributions and by the timetree message passes,
/// which assign a per-side tail policy to each message.
///
/// Every variant is a complete value: the two law-carrying variants ([`Self::HardApproach`],
/// [`Self::Linear`]) always hold a fitted law, and the two nullary variants ([`Self::Error`],
/// [`Self::Hard`]) declare a boundary that needs no law. A grid too small to fit a required law is
/// an error at the fitting site, never a silent flat fallback, so no variant stands for "a law was
/// wanted here but is missing".
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum BoundaryBehavior {
  /// Out-of-support evaluation is an error.
  #[default]
  Error,
  /// Hard boundary with no sub-grid gap: the domain terminates at the grid edge and probability is
  /// zero beyond it. The raw `GridFn` returns `T::zero()` outside support; the distribution layer
  /// maps this to the policy-correct zero-probability value (e.g. `+inf` under negative-log
  /// representation). This is the correct value when the grid edge *is* the hard boundary, so there
  /// is no gap between them to interpolate.
  Hard,
  /// Hard boundary with a fitted approach law across the sub-grid gap. The domain still terminates
  /// (probability is zero past [`HardApproachLaw::t_hard`]), but between the hard boundary and the
  /// nearest grid point the density follows the edge-relative [`HardApproachLaw`]: the single-exponent
  /// power law of a divergent branch. A finite mode sitting exactly on the boundary is instead stored
  /// as an exact grid endpoint ([`Self::Hard`]), so this variant carries only the divergent case.
  HardApproach(HardApproachLaw),
  /// Soft boundary: the density continues past the grid edge along the log-linear [`SoftTailLaw`]
  /// (an exponential probability tail). The tail decays and has finite mass, so it keeps the
  /// quantile and HPD integrals well-defined.
  Linear(SoftTailLaw),
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
  /// `Linear` (log-linear) is the only soft tail. This predicate is what the multiplication rule
  /// keys off, so a soft edge extends the domain through that rule.
  pub fn is_soft(self) -> bool {
    matches!(self, BoundaryBehavior::Linear(_))
  }

  /// Whether this boundary terminates the evaluable domain on its side: a `Hard`/`HardApproach`
  /// boundary is zero probability beyond the edge and an `Error` boundary is undefined beyond it, so
  /// all three restrict the result domain under multiplication. The complement is [`Self::is_soft`].
  pub fn is_hard(self) -> bool {
    matches!(
      self,
      BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_) | BoundaryBehavior::Error
    )
  }

  /// The fitted power-law approach for a hard boundary, if the side carries one.
  pub fn approach_law(self) -> Option<HardApproachLaw> {
    match self {
      BoundaryBehavior::HardApproach(law) => Some(law),
      _ => None,
    }
  }

  /// The fitted log-linear tail for a `Linear` soft boundary, if the side carries one.
  pub fn soft_law(self) -> Option<SoftTailLaw> {
    match self {
      BoundaryBehavior::Linear(law) => Some(law),
      _ => None,
    }
  }
}

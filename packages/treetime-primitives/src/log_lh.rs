use serde::{Deserialize, Serialize};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Neg, Sub};

/// Aggregate log-likelihood in natural-log space.
///
/// The wrapper distinguishes complete and component log-likelihoods from raw
/// probabilities, branch lengths, dates, and optimizer coordinates while
/// retaining the layout of `f64`.
///
/// Raw floats cannot be mixed into log-likelihood arithmetic implicitly:
///
/// ```compile_fail
/// use treetime_primitives::LogLh;
///
/// let _ = LogLh::new(-2.0) + 1.0;
/// ```
///
/// ```compile_fail
/// use treetime_primitives::LogLh;
///
/// let _ = 1.0 + LogLh::new(-2.0);
/// ```
///
/// ```compile_fail
/// use treetime_primitives::LogLh;
///
/// let _: LogLh = -2.0;
/// ```
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
#[repr(transparent)]
#[serde(transparent)]
pub struct LogLh(f64);

impl LogLh {
  /// Additive identity for independent log-likelihood components.
  pub const ZERO: Self = Self(0.0);

  /// Log-likelihood of an impossible event.
  pub const IMPOSSIBLE: Self = Self(f64::NEG_INFINITY);

  pub const fn new(value: f64) -> Self {
    Self(value)
  }

  pub const fn value(self) -> f64 {
    self.0
  }
}

impl Default for LogLh {
  fn default() -> Self {
    Self::ZERO
  }
}

impl Add for LogLh {
  type Output = Self;

  fn add(self, rhs: Self) -> Self::Output {
    Self(self.0 + rhs.0)
  }
}

impl AddAssign for LogLh {
  fn add_assign(&mut self, rhs: Self) {
    self.0 += rhs.0;
  }
}

impl Sum for LogLh {
  fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
    iter.fold(Self::ZERO, Add::add)
  }
}

impl<'a> Sum<&'a LogLh> for LogLh {
  fn sum<I: Iterator<Item = &'a LogLh>>(iter: I) -> Self {
    iter.copied().sum()
  }
}

impl Sub for LogLh {
  type Output = f64;

  fn sub(self, rhs: Self) -> Self::Output {
    self.0 - rhs.0
  }
}

impl Neg for LogLh {
  type Output = f64;

  fn neg(self) -> Self::Output {
    -self.0
  }
}

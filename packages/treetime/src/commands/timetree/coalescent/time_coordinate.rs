use std::ops::{Add, Sub};

/// Calendar time as decimal year fraction (e.g., 2013.5 = July 2013).
///
/// Primary time coordinate in v1. Date constraints, node times, and distribution
/// domains all use calendar time. The coalescent module converts to [`Tbp`] for
/// internal merger rate computations and back to `CalendarTime` for distribution
/// domains consumed by the backward pass.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct CalendarTime(f64);

impl CalendarTime {
  pub fn new(value: f64) -> Self {
    Self(value)
  }

  /// Convert to time-before-present relative to `present`.
  pub fn to_tbp(self, present: CalendarTime) -> Tbp {
    Tbp(present.0 - self.0)
  }

  /// Raw f64 for coordinate-agnostic APIs (distributions, arrays).
  pub fn value(self) -> f64 {
    self.0
  }

  pub fn max(self, other: Self) -> Self {
    if self.0 >= other.0 { self } else { other }
  }

  pub fn is_finite(self) -> bool {
    self.0.is_finite()
  }
}

/// Time before present (TBP). Zero at the most recent sample, increases into the past.
///
/// Used internally by the coalescent module for lineage counts, merger rates,
/// and integral computations. All piecewise functions operate in TBP.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct Tbp(f64);

impl Tbp {
  pub fn new(value: f64) -> Self {
    Self(value)
  }

  /// Convert to calendar time relative to `present`.
  pub fn to_calendar(self, present: CalendarTime) -> CalendarTime {
    CalendarTime(present.0 - self.0)
  }

  /// Raw f64 for coordinate-agnostic APIs (piecewise functions, arrays).
  pub fn value(self) -> f64 {
    self.0
  }
}

/// `Tbp - Tbp` yields a duration (plain f64).
impl Sub for Tbp {
  type Output = f64;

  fn sub(self, rhs: Tbp) -> f64 {
    self.0 - rhs.0
  }
}

/// `Tbp + duration` yields a shifted `Tbp`.
impl Add<f64> for Tbp {
  type Output = Tbp;

  fn add(self, duration: f64) -> Tbp {
    Tbp(self.0 + duration)
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;

  #[test]
  fn test_calendar_to_tbp_roundtrip() {
    let present = CalendarTime::new(2013.0);
    let cal = CalendarTime::new(2005.0);
    let tbp = cal.to_tbp(present);
    assert_ulps_eq!(tbp.value(), 8.0, max_ulps = 4);
    let back = tbp.to_calendar(present);
    assert_ulps_eq!(back.value(), 2005.0, max_ulps = 4);
  }

  #[test]
  fn test_tbp_at_present_is_zero() {
    let present = CalendarTime::new(2013.0);
    let tbp = present.to_tbp(present);
    assert_ulps_eq!(tbp.value(), 0.0, max_ulps = 4);
  }

  #[test]
  fn test_tbp_arithmetic() {
    let a = Tbp::new(8.0);
    let b = Tbp::new(3.0);
    assert_ulps_eq!(a - b, 5.0, max_ulps = 4);
    assert_ulps_eq!((b + 2.0).value(), 5.0, max_ulps = 4);
  }

  #[test]
  fn test_calendar_time_max() {
    let a = CalendarTime::new(2005.0);
    let b = CalendarTime::new(2013.0);
    assert_eq!(a.max(b), b);
    assert_eq!(b.max(a), b);
  }
}

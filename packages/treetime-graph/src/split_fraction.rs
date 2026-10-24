use eyre::Report;
use treetime_utils::make_internal_report;

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct SplitFraction(f64);

impl SplitFraction {
  pub const PARENT_END: Self = Self(0.0);

  pub const MIDDLE: Self = Self(0.5);

  pub const CHILD_END: Self = Self(1.0);

  pub fn new(value: f64) -> Result<Self, Report> {
    if (0.0..=1.0).contains(&value) {
      Ok(Self(value))
    } else {
      Err(make_internal_report!(
        "A split fraction must be a finite number in [0, 1], got {value}"
      ))
    }
  }

  #[must_use]
  pub fn value(self) -> f64 {
    self.0
  }

  #[must_use]
  pub fn lengths(self, length: f64) -> SplitLengths {
    SplitLengths {
      parent_side: self.0 * length,
      child_side: (1.0 - self.0) * length,
    }
  }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SplitLengths {
  pub parent_side: f64,
  pub child_side: f64,
}

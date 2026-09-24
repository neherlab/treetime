#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
#[repr(transparent)]
pub(crate) struct CalendarTime(f64);

impl CalendarTime {
  pub(crate) fn new(value: f64) -> Self {
    Self(value)
  }

  pub(crate) fn value(self) -> f64 {
    self.0
  }

  pub(crate) fn max(self, other: Self) -> Self {
    if self.0 >= other.0 { self } else { other }
  }

  pub(crate) fn is_finite(self) -> bool {
    self.0.is_finite()
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_time_coordinate_calendar_time_max() {
    let earlier = CalendarTime::new(2005.0);
    let later = CalendarTime::new(2013.0);
    assert_eq!(earlier.max(later), later);
    assert_eq!(later.max(earlier), later);
  }
}

use treetime_utils::cast::f64_to_i32_saturating;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearScale {
  domain: (f64, f64),
  range: (f64, f64),
}

impl LinearScale {
  pub fn new(domain: (f64, f64), range: (f64, f64)) -> Self {
    Self { domain, range }
  }

  pub fn map(&self, value: f64) -> f64 {
    let (d0, d1) = self.domain;
    let (r0, r1) = self.range;
    let span = d1 - d0;
    if span == 0.0 {
      return f64::midpoint(r0, r1);
    }
    r0 + (value - d0) / span * (r1 - r0)
  }

  pub fn pixel(&self, value: f64) -> i32 {
    f64_to_i32_saturating(self.map(value).round())
  }
}

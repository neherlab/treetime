use itertools::Itertools;
use std::iter::repeat_n;
use std::ops::Add;

#[derive(Clone, Debug)]
pub(crate) struct CellSums {
  nodes: Vec<ScaledMass>,
}

impl CellSums {
  pub(crate) fn new(cells: impl ExactSizeIterator<Item = ScaledMass>) -> Self {
    let n = cells.len();
    let mut nodes = repeat_n(ScaledMass::ZERO, n).chain(cells).collect_vec();
    for i in (1..n).rev() {
      nodes[i] = nodes[2 * i] + nodes[2 * i + 1];
    }
    Self { nodes }
  }

  pub(crate) fn sum(&self, from: usize, to: usize) -> ScaledMass {
    let n = self.nodes.len() >> 1;
    let (mut lo, mut hi) = (from + n, to + n);
    let mut sum = ScaledMass::ZERO;
    while lo < hi {
      if lo & 1 == 1 {
        sum = sum + self.nodes[lo];
        lo += 1;
      }
      if hi & 1 == 1 {
        hi -= 1;
        sum = sum + self.nodes[hi];
      }
      lo >>= 1;
      hi >>= 1;
    }
    sum
  }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct ScaledMass {
  r: f64,
  factor: f64,
}

impl ScaledMass {
  pub(crate) const ZERO: Self = Self {
    r: f64::INFINITY,
    factor: 0.0,
  };

  pub(crate) const INFINITE: Self = Self {
    r: f64::NEG_INFINITY,
    factor: 1.0,
  };

  pub(crate) fn new(r: f64, factor: f64) -> Self {
    if factor > 0.0 && r < f64::INFINITY {
      Self { r, factor }
    } else {
      Self::ZERO
    }
  }

  pub(crate) fn linear(self) -> f64 {
    if self.factor == 0.0 {
      return 0.0;
    }
    let scale = (-self.r).exp();
    if scale.is_normal() {
      self.factor * scale
    } else {
      (self.factor.ln() - self.r).exp()
    }
  }

  pub(crate) fn neg_log(self) -> f64 {
    self.r - self.factor.ln()
  }

  fn ln_mass(self) -> f64 {
    -self.neg_log()
  }

  fn is_infinite(self) -> bool {
    self.factor > 0.0 && (self.r == f64::NEG_INFINITY || (self.factor == f64::INFINITY && self.r != f64::INFINITY))
  }
}

impl Add for ScaledMass {
  type Output = Self;

  fn add(self, other: Self) -> Self {
    if self.factor == 0.0 || other.is_infinite() {
      return other;
    }
    if other.factor == 0.0 || self.is_infinite() {
      return self;
    }
    let r = self.r.min(other.r);
    let (scale, other_scale) = (scale_below(r, self.r), scale_below(r, other.r));
    if scale.is_normal() && other_scale.is_normal() {
      return Self {
        r,
        factor: self.factor * scale + other.factor * other_scale,
      };
    }
    let (large, small) = if self.ln_mass() >= other.ln_mass() {
      (self, other)
    } else {
      (other, self)
    };
    Self {
      r: large.r,
      factor: large.factor * (1.0 + (small.ln_mass() - large.ln_mass()).exp()),
    }
  }
}

#[expect(
  clippy::float_cmp,
  reason = "r_min is the minimum of r and another ordinate, so exact equality selects the operand whose scale is exactly 1"
)]
pub(crate) fn scale_below(r_min: f64, r: f64) -> f64 {
  if r == r_min && r.is_finite() {
    1.0
  } else {
    (r_min - r).exp()
  }
}

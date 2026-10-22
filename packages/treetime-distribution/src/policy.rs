use std::fmt::Debug;
use std::marker::PhantomData;

const MAX_VALUE_ULPS: f64 = 10.0;

pub trait YAxisPolicy: Clone + Copy + Debug + Default + PartialEq + Send + Sync + 'static {
  fn from_plain(p: f64) -> f64;
  fn to_plain(y: f64) -> f64;
  fn multiply(a: f64, b: f64) -> f64;
  fn divide(a: f64, b: f64) -> f64;
  fn is_defined(val: f64) -> bool;
  fn safe_divisor(val: f64) -> f64;

  fn probability_zero() -> f64;

  fn likely_is_maximum() -> bool;

  fn to_neg_log(y: f64) -> f64;

  fn from_neg_log(nl: f64) -> f64;

  fn same_value(a: f64, b: f64) -> bool;
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Plain;

impl YAxisPolicy for Plain {
  fn from_plain(p: f64) -> f64 {
    p
  }

  fn to_plain(y: f64) -> f64 {
    y
  }

  fn multiply(a: f64, b: f64) -> f64 {
    a * b
  }

  fn divide(a: f64, b: f64) -> f64 {
    a / b
  }

  fn is_defined(val: f64) -> bool {
    val > 0.0
  }

  fn safe_divisor(val: f64) -> f64 {
    const TINY_NUMBER: f64 = 1e-10;
    val.max(TINY_NUMBER)
  }

  fn probability_zero() -> f64 {
    0.0
  }

  fn likely_is_maximum() -> bool {
    true
  }

  fn to_neg_log(y: f64) -> f64 {
    if y > 0.0 { -y.ln() } else { f64::INFINITY }
  }

  fn from_neg_log(nl: f64) -> f64 {
    (-nl).exp()
  }

  fn same_value(a: f64, b: f64) -> bool {
    values_coincide(a, b, a.abs().max(b.abs()))
  }
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NegLog;

impl YAxisPolicy for NegLog {
  fn from_plain(p: f64) -> f64 {
    -p.ln()
  }

  fn to_plain(y: f64) -> f64 {
    (-y).exp()
  }

  fn multiply(a: f64, b: f64) -> f64 {
    a + b
  }

  fn divide(a: f64, b: f64) -> f64 {
    a - b
  }

  fn is_defined(val: f64) -> bool {
    val.is_finite()
  }

  fn safe_divisor(val: f64) -> f64 {
    val
  }

  fn probability_zero() -> f64 {
    f64::INFINITY
  }

  fn likely_is_maximum() -> bool {
    false
  }

  fn to_neg_log(y: f64) -> f64 {
    y
  }

  fn from_neg_log(nl: f64) -> f64 {
    nl
  }

  fn same_value(a: f64, b: f64) -> bool {
    values_coincide(a, b, a.abs().max(b.abs()).max(1.0))
  }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PolicyMarker<Y: YAxisPolicy>(PhantomData<Y>);

impl<Y: YAxisPolicy> PolicyMarker<Y> {
  pub(crate) fn new() -> Self {
    Self(PhantomData)
  }
}

fn values_coincide(a: f64, b: f64, scale: f64) -> bool {
  a == b || (a.is_finite() && b.is_finite() && (a - b).abs() <= MAX_VALUE_ULPS * f64::EPSILON * scale)
}

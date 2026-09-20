use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::marker::PhantomData;

pub trait YAxisPolicy: Clone + Copy + Debug + Default + PartialEq + Send + Sync + 'static {
  fn from_plain(p: f64) -> f64;
  fn to_plain(y: f64) -> f64;
  fn multiplicative_identity() -> f64;
  fn multiply(a: f64, b: f64) -> f64;
  fn divide(a: f64, b: f64) -> f64;
  fn is_defined(val: f64) -> bool;
  fn safe_divisor(val: f64) -> f64;

  fn probability_zero() -> f64;

  fn likely_is_maximum() -> bool;

  fn to_neg_log(y: f64) -> f64;

  fn from_neg_log(nl: f64) -> f64;
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Plain;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct NegLog;

pub trait SupportsConvolution: YAxisPolicy {}

impl SupportsConvolution for Plain {}
impl SupportsConvolution for NegLog {}

pub trait SupportsSubtraction: YAxisPolicy {}

impl SupportsSubtraction for Plain {}

impl YAxisPolicy for Plain {
  fn from_plain(p: f64) -> f64 {
    p
  }

  fn to_plain(y: f64) -> f64 {
    y
  }

  fn multiplicative_identity() -> f64 {
    1.0
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
}

impl YAxisPolicy for NegLog {
  fn from_plain(p: f64) -> f64 {
    -p.ln()
  }

  fn to_plain(y: f64) -> f64 {
    (-y).exp()
  }

  fn multiplicative_identity() -> f64 {
    0.0
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
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct PolicyMarker<Y: YAxisPolicy>(#[serde(skip)] PhantomData<Y>);

impl<Y: YAxisPolicy> PolicyMarker<Y> {
  pub fn new() -> Self {
    Self(PhantomData)
  }
}

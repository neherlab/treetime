use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::marker::PhantomData;

pub trait YAxisPolicy: Clone + Copy + Debug + Default + PartialEq + Send + Sync + 'static {
  fn from_plain(p: f64) -> f64;
  fn to_plain(y: f64) -> f64;
  fn multiplicative_identity() -> f64;
  fn multiply(a: f64, b: f64) -> f64;
  fn divide(a: f64, b: f64) -> f64;
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Plain;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct NegLog;

/// Marker trait for policies that support convolution operations.
/// Convolution requires integral sums which are not representable in log space.
pub trait SupportsConvolution: YAxisPolicy {}

impl SupportsConvolution for Plain {}

/// Marker trait for policies that support subtraction operations.
/// Subtraction uses direct y value subtraction which is only valid for plain values.
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
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct PolicyMarker<Y: YAxisPolicy>(#[serde(skip)] PhantomData<Y>);

impl<Y: YAxisPolicy> PolicyMarker<Y> {
  pub fn new() -> Self {
    Self(PhantomData)
  }
}

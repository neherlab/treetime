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

  /// Whether a `Hard` boundary tail is representable under this policy.
  ///
  /// A `Hard` boundary tail writes the literal `0.0` outside support. That is zero
  /// probability under [`Plain`], but the multiplicative identity (probability one) under
  /// [`NegLog`], where zero probability is `+inf`. The distribution layer rejects a `Hard`
  /// tail on a policy that returns `false` here.
  fn supports_hard_boundary() -> bool;

  /// Whether the most likely time sits at the maximum stored ordinate.
  ///
  /// Under [`Plain`] the ordinate is probability, so the mode is the maximum ordinate. Under
  /// [`NegLog`] the ordinate is `-ln(probability)`, so the mode is the minimum ordinate.
  /// Most-likely-time selection dispatches on this instead of assuming a maximum.
  fn likely_is_maximum() -> bool;

  /// Map a stored ordinate to negative-log probability (`-ln p`).
  ///
  /// This is the common space where the convolution tail law is linear (an exponential tail in
  /// probability is a straight line in `-ln p`). [`Plain`] takes `-ln`, treating a non-positive
  /// probability as `+inf` (zero probability); [`NegLog`] is already there and is the identity.
  fn to_neg_log(y: f64) -> f64;

  /// Inverse of [`Self::to_neg_log`]: map a negative-log value back to a stored ordinate.
  ///
  /// [`Plain`] exponentiates (`exp(-nl)`), so values below the plain underflow floor collapse to
  /// zero; [`NegLog`] is the identity and preserves the full dynamic range. This asymmetry is why
  /// log-space storage can represent tail values the plain axis cannot.
  fn from_neg_log(nl: f64) -> f64;
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Plain;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct NegLog;

/// Marker trait for policies whose convolution is defined.
///
/// The convolution integral runs in plain probability space regardless of storage policy: each
/// operand is converted through [`YAxisPolicy::to_neg_log`] and peak-normalized to plain, convolved,
/// then converted back through [`YAxisPolicy::from_neg_log`]. Both [`Plain`] and [`NegLog`] support
/// this round-trip.
pub trait SupportsConvolution: YAxisPolicy {}

impl SupportsConvolution for Plain {}
impl SupportsConvolution for NegLog {}

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

  fn is_defined(val: f64) -> bool {
    val > 0.0
  }

  fn safe_divisor(val: f64) -> f64 {
    const TINY_NUMBER: f64 = 1e-10;
    val.max(TINY_NUMBER)
  }

  fn supports_hard_boundary() -> bool {
    true
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

  fn supports_hard_boundary() -> bool {
    false
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

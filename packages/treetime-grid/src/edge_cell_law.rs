use crate::coordinates_coincide;
use eyre::Report;
use treetime_utils::make_internal_error;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeCellLaw {
  pub t_edge: f64,
  pub offset: f64,
  pub b: f64,
}

impl EdgeCellLaw {
  #[must_use]
  pub fn scale_arg(&self, factor: f64) -> Self {
    Self {
      t_edge: self.t_edge * factor,
      offset: self.offset * factor.abs(),
      b: self.b,
    }
  }

  #[must_use]
  pub fn shift_arg(&self, shift: f64) -> Self {
    Self {
      t_edge: self.t_edge + shift,
      ..*self
    }
  }

  pub fn at_edge(&self, t: f64) -> bool {
    coordinates_coincide(self.t_edge, t, self.t_edge.abs().max(self.offset))
  }

  #[must_use]
  pub fn compose_multiply(&self, other: &Self) -> Option<Self> {
    if !self.at_edge(other.t_edge) {
      return None;
    }
    let b = self.b + other.b;
    if b <= 0.0 || self.same_offset(other) {
      return Some(Self { b, ..*self });
    }
    let ln_offset = (self.b * self.offset.ln() + other.b * other.offset.ln()) / b;
    Some(Self {
      offset: ln_offset.exp(),
      b,
      ..*self
    })
  }

  pub fn compose_divide(&self, divisor: &Self) -> Result<Option<Self>, Report> {
    if !self.at_edge(divisor.t_edge) {
      return Ok(None);
    }
    let b = self.b - divisor.b;
    if coordinates_coincide(self.b, divisor.b, self.b) {
      return Ok(Some(Self { b: 0.0, ..*self }));
    }
    if b < 0.0 {
      return make_internal_error!(
        "Cannot divide an edge cell law of exponent {} by one of exponent {}: the quotient has a negative exponent",
        self.b,
        divisor.b
      );
    }
    if self.same_offset(divisor) {
      return Ok(Some(Self { b, ..*self }));
    }
    let ln_offset = (self.b * self.offset.ln() - divisor.b * divisor.offset.ln()) / b;
    Ok(Some(Self {
      offset: ln_offset.exp(),
      b,
      ..*self
    }))
  }

  fn same_offset(&self, other: &Self) -> bool {
    coordinates_coincide(self.offset, other.offset, self.offset)
  }
}

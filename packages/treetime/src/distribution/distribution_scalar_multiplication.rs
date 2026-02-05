use crate::distribution::distribution::Distribution;
use crate::distribution::y_axis_policy::YAxisPolicy;
use eyre::Report;

/// Multiply a distribution by a scalar value.
///
/// For all distribution types, multiply the amplitude values by the scalar.
/// This scales the probability density or likelihood represented by the distribution.
pub fn distribution_scalar_multiplication<Y: YAxisPolicy>(
  dist: &Distribution<Y>,
  scalar: f64,
) -> Result<Distribution<Y>, Report> {
  match dist {
    Distribution::Function(f) => {
      let new_y = f.y().mapv(|y| Y::multiply(y, Y::from_plain(scalar)));
      Distribution::function(f.t(), new_y) //
    },
    Distribution::Point(p) => {
      let amplitude = Y::multiply(p.amplitude(), Y::from_plain(scalar));
      Ok(Distribution::point(p.t(), amplitude)) //
    },
    Distribution::Range(r) => {
      let amplitude = Y::multiply(r.amplitude(), Y::from_plain(scalar));
      Ok(Distribution::range((r.start(), r.end()), amplitude)) //
    },
    Distribution::Empty => {
      Ok(Distribution::empty()) //
    },
    Distribution::Formula(_) => {
      panic!("Scalar multiplication not implemented for Formula distributions") //
    },
  }
}

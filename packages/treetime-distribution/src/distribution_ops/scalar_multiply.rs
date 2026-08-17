use crate::Distribution;
use crate::distribution_core::function::DistributionFunction;
use crate::policy::YAxisPolicy;
use eyre::Report;
use treetime_utils::make_error;

/// Multiply a distribution by a scalar value.
///
/// For all distribution types, multiply the amplitude values by the scalar.
/// This scales the probability density or likelihood represented by the distribution.
///
/// Tail policy survives (`kb/decisions/distribution-tails-and-arithmetic.md`). Scaling every
/// probability by a constant is a constant shift of the negative-log ordinate, so both fitted laws
/// carry through unchanged: a `SoftTailLaw` stores only its neg-log slope and a `HardApproachLaw`
/// only its shape, and neither depends on the constant scale; each law re-reads the now-scaled edge
/// ordinate on evaluation. This mirrors how `GridFn::shift_y`/`scale_y` preserve laws, whereas
/// rebuilding the grid alone would erase them to `Error`.
pub fn distribution_scalar_multiplication<Y: YAxisPolicy>(
  dist: &Distribution<Y>,
  scalar: f64,
) -> Result<Distribution<Y>, Report> {
  match dist {
    Distribution::Function(f) => {
      let new_y = f.y().mapv(|y| Y::multiply(y, Y::from_plain(scalar)));
      DistributionFunction::from_start_dx_values(f.x_min(), f.dx(), new_y)?
        .with_left_extrap(f.left_extrap())?
        .with_right_extrap(f.right_extrap())
        .map(Distribution::Function)
    },
    Distribution::Point(p) => {
      let amplitude = Y::multiply(p.amplitude(), Y::from_plain(scalar));
      Ok(Distribution::point(p.t(), amplitude))
    },
    Distribution::Range(r) => {
      let amplitude = Y::multiply(r.amplitude(), Y::from_plain(scalar));
      Ok(Distribution::range((r.start(), r.end()), amplitude))
    },
    Distribution::Empty => Ok(Distribution::empty()),
    Distribution::Formula(_) => make_error!("Cannot multiply Formula by scalar: operation not implemented"),
  }
}

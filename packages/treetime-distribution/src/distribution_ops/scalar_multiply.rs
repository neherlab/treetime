use crate::Distribution;
use crate::distribution_core::function::DistributionFunction;
use crate::policy::YAxisPolicy;
use eyre::Report;
use treetime_utils::make_error;

pub(crate) fn distribution_scalar_multiplication<Y: YAxisPolicy>(
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

use crate::Distribution;
use crate::distribution_core::function::DistributionFunction;
use crate::policy::SupportsSubtraction;
use eyre::Report;
use treetime_utils::make_error;

pub fn distribution_subtraction<Y: SupportsSubtraction>(
  a: &Distribution<Y>,
  b: &Distribution<Y>,
) -> Result<Distribution<Y>, Report> {
  match (a, b) {
    (Distribution::Function(af), Distribution::Function(bf)) => {
      if af.grid() != bf.grid() {
        return make_error!("Cannot subtract distributions with different grids");
      }
      DistributionFunction::from_start_dx_values(af.x_min(), af.dx(), af.y() - bf.y()).map(Distribution::Function)
    },
    (Distribution::Empty | Distribution::Point(_) | Distribution::Range(_) | Distribution::Formula(_), _)
    | (_, Distribution::Empty | Distribution::Point(_) | Distribution::Range(_) | Distribution::Formula(_)) => {
      make_error!("Subtraction only supported for Function distributions with matching grids")
    },
  }
}

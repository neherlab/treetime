use crate::distribution::distribution::Distribution;
use crate::distribution::y_axis_policy::SupportsSubtraction;
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
      Distribution::function(af.t(), af.y() - bf.y())
    },
    (Distribution::Empty | Distribution::Point(_) | Distribution::Range(_) | Distribution::Formula(_), _)
    | (_, Distribution::Empty | Distribution::Point(_) | Distribution::Range(_) | Distribution::Formula(_)) => {
      make_error!("Subtraction only supported for Function distributions with matching grids")
    },
  }
}

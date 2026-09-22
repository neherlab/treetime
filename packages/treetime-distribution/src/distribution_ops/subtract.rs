use crate::Distribution;
use crate::distribution_core::function::DistributionFunction;
use crate::policy::SupportsSubtraction;
use eyre::Report;
use treetime_grid::BoundaryBehavior;
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
      let left = subtraction_result_tail(af.left_extrap(), bf.left_extrap());
      let right = subtraction_result_tail(af.right_extrap(), bf.right_extrap());
      DistributionFunction::from_start_dx_values(af.x_min(), af.dx(), af.y() - bf.y())?
        .with_left_extrap(left)?
        .with_right_extrap(right)
        .map(Distribution::Function)
    },
    (Distribution::Empty | Distribution::Point(_) | Distribution::Range(_) | Distribution::Formula(_), _)
    | (_, Distribution::Empty | Distribution::Point(_) | Distribution::Range(_) | Distribution::Formula(_)) => {
      make_error!("Subtraction only supported for Function distributions with matching grids")
    },
  }
}

pub fn subtraction_result_tail(a: BoundaryBehavior, b: BoundaryBehavior) -> BoundaryBehavior {
  let zero_beyond = |tail| matches!(tail, BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_));
  if zero_beyond(a) && zero_beyond(b) {
    BoundaryBehavior::Hard
  } else {
    BoundaryBehavior::Error
  }
}

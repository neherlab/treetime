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
      // Tail policy survives (kb/decisions/distribution-tails-and-arithmetic.md). A pointwise
      // difference of two densities is not log-linear, so no soft tail law can represent it: a side
      // becomes `Hard` (zero beyond) only when both operands vanish beyond the edge (`Hard` or
      // `HardApproach`), otherwise `Error`. This preserves the genuine zero-beyond fact without
      // fabricating a slope the difference cannot carry, rather than erasing every side to `Error`.
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

/// Per-side tail of a Function difference. Zero-beyond-hard on both operands (`Hard`/`HardApproach`)
/// gives a `Hard` result (`0 - 0 = 0`); any soft or `Error` operand leaves the difference without a
/// representable tail law, so the side is `Error`.
fn subtraction_result_tail(a: BoundaryBehavior, b: BoundaryBehavior) -> BoundaryBehavior {
  let zero_beyond = |tail| matches!(tail, BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_));
  if zero_beyond(a) && zero_beyond(b) {
    BoundaryBehavior::Hard
  } else {
    BoundaryBehavior::Error
  }
}

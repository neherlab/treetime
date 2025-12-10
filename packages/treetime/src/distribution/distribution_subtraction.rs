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
      if af.t().len() != bf.t().len() {
        return make_error!("Cannot subtract distributions with different time points");
      }
      Distribution::function(af.t(), af.y() - bf.y())
    },
    _ => make_error!("Subtraction only supported for Function distributions with matching time points"),
  }
}

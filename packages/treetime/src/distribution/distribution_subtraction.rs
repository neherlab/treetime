use crate::distribution::distribution::Distribution;
use eyre::Report;
use treetime_utils::make_error;

pub fn distribution_subtraction(a: &Distribution, b: &Distribution) -> Result<Distribution, Report> {
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

#[cfg(test)]
pub(super) mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::optimize;
  use crate::partition::optimize::contribution::OptimizationContribution;
  use ndarray::Array2;

  pub(crate) fn make_dense_contribution(coefficients: Array2<f64>) -> OptimizationContribution {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");
    OptimizationContribution::Dense(optimize::dense::PartitionContribution::new(coefficients, gtr))
  }
}

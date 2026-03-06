#[cfg(test)]
pub mod tests {
  use crate::commands::optimize::optimize_dense;
  use crate::commands::optimize::optimize_unified::OptimizationContribution;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use ndarray::Array2;

  pub fn make_dense_contribution(coefficients: Array2<f64>) -> OptimizationContribution {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");
    OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr))
  }
}

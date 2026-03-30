#[cfg(test)]
pub mod tests {
  use crate::commands::optimize::optimize_dense;
  use crate::commands::optimize::optimize_unified::{OptimizationContribution, evaluate_mixed};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use ndarray::{Array1, Array2};
  use ordered_float::OrderedFloat;

  pub fn make_dense_contribution(coefficients: Array2<f64>) -> OptimizationContribution {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");
    OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr))
  }

  pub fn grid_search(contributions: &[OptimizationContribution], branch_length: f64, one_mutation: f64) -> f64 {
    let branch_lengths = Array1::linspace(0.0, 1.5 * branch_length + one_mutation, 100);

    branch_lengths
      .iter()
      .max_by_key(|&&bl| {
        let metrics = evaluate_mixed(contributions, bl);
        OrderedFloat(metrics.log_lh)
      })
      .copied()
      .unwrap_or(branch_length)
  }
}

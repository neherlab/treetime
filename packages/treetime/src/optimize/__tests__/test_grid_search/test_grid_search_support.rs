#[cfg(test)]
pub mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::optimize_unified::{evaluate_mixed, grid_search_branch_lengths};
  use crate::representation::partition::optimization_contribution::OptimizationContribution;
  use crate::representation::partition::optimize_dense;
  use ndarray::Array2;
  use ordered_float::OrderedFloat;

  pub fn make_dense_contribution(coefficients: Array2<f64>) -> OptimizationContribution {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");
    OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr))
  }

  pub fn grid_search(contributions: &[OptimizationContribution], branch_length: f64, one_mutation: f64) -> f64 {
    let branch_lengths = grid_search_branch_lengths(branch_length, one_mutation).unwrap();

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

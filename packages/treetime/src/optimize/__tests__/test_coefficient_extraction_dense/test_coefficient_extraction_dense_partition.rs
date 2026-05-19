#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::optimize_dense::PartitionContribution;
  use ndarray::array;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_partition_contribution_new() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");
    let coefficients = array![[0.1, 0.2, 0.3, 0.4], [0.4, 0.3, 0.2, 0.1]];

    let contribution = PartitionContribution::new(coefficients.clone(), gtr);

    pretty_assert_ulps_eq!(contribution.coefficients, coefficients, max_ulps = 10);
  }
}

#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_dense::PartitionContribution;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use approx::assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_partition_contribution_new() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");
    let coefficients = array![[0.1, 0.2, 0.3, 0.4], [0.4, 0.3, 0.2, 0.1]];

    let contribution = PartitionContribution::new(coefficients.clone(), gtr);

    assert_eq!(contribution.coefficients.shape(), coefficients.shape());
    for i in 0..2 {
      for j in 0..4 {
        assert_ulps_eq!(contribution.coefficients[[i, j]], coefficients[[i, j]], max_ulps = 10);
      }
    }
  }
}

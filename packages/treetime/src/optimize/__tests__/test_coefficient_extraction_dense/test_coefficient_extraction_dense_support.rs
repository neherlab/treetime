#[cfg(test)]
pub mod tests {
  use crate::partition::storage::dense::DenseSeqDistribution;
  use ndarray::Array2;
  use treetime_primitives::LogLh;

  pub fn make_dense_seq_dis(dis: Array2<f64>) -> DenseSeqDistribution {
    DenseSeqDistribution::new(dis, LogLh::ZERO)
  }
}

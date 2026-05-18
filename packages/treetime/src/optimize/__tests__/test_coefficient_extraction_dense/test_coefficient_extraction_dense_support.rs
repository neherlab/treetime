#[cfg(test)]
pub mod tests {
  use crate::partition::dense::DenseSeqDistribution;
  use ndarray::Array2;

  pub fn make_dense_seq_dis(dis: Array2<f64>) -> DenseSeqDistribution {
    DenseSeqDistribution::new(dis, 0.0)
  }
}

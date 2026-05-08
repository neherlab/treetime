#[cfg(test)]
pub mod tests {
  use crate::representation::payload::dense::DenseSeqDis;
  use ndarray::Array2;

  pub fn make_dense_seq_dis(dis: Array2<f64>) -> DenseSeqDis {
    DenseSeqDis::new(dis, 0.0)
  }
}

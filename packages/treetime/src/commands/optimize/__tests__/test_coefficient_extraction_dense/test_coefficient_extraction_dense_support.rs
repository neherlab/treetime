use crate::representation::payload::dense::DenseSeqDis;
use ndarray::Array2;

pub(super) fn make_dense_seq_dis(dis: Array2<f64>) -> DenseSeqDis {
  DenseSeqDis::new(dis)
}

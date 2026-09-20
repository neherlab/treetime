use crate::gtr::gtr::GTR;
use crate::partition::storage::dense::DenseSeqDistribution;
use ndarray::Array2;

pub struct PartitionContribution {
  pub coefficients: Array2<f64>,
  pub gtr: GTR,
}

impl PartitionContribution {
  pub fn new(coefficients: Array2<f64>, gtr: GTR) -> Self {
    Self { coefficients, gtr }
  }
}

pub fn get_coefficients(
  msg_to_parent: &DenseSeqDistribution,
  msg_to_child: &DenseSeqDistribution,
  gtr: &GTR,
) -> PartitionContribution {
  let coefficients = msg_to_child.dis.dot(&gtr.v) * msg_to_parent.dis.dot(&gtr.v_inv.t());
  PartitionContribution::new(coefficients, gtr.clone())
}

//! The likelihood of an edge length is the product of the likelihoods of all positions of all partitions
//!
//!   Lh = prod_i prod_j \sum_{ab} s^{ij}_a exp(Q_i t)_{ab} r^{ij}_b
//!
//! The log likelihood is the sum of many terms:
//!
//!   logLh = sum_i sum_j \log(\sum_{ab} s^{ij}_a exp(Q_i t)_{ab} r^{ij}_b)
//!
//! To effectively calculate this, we need to reformulate the likelihood in terms of the eigenvectors of the GTR matrix. Dropping the {ij} superscripts for brevity, we can write the likelihood as:
//!
//!   s_a exp(Qt)_{ab} r_b = s_a \sum_c v_{ac} exp(\lambda_c t) vinv_{cb} r_b = g_c exp(\lambda_c t) h_c = k_c exp(\lambda_c t)
//!
//! The `k_c` can be reused for different iterations of the branch length optimization
//!
//!   logLh = sum_i sum_j \log(\sum_c k_c exp(\lambda^i_c t))
//!
//! The derivative is simply:
//!
//!   dlogLh/dt = sum_i sum_j \sum_c k_c \lambda_c exp(\lambda^i_c t) / \sum_c k_c exp(\lambda^i_c t)
//!
//!   d^2logLh/dt^2 = sum_i sum_j \sum_c k_c \lambda_c*\lambda^i_c exp(\lambda^i_c t) / \sum_c k_c exp(\lambda^i_c t) - k_c \lambda_c*\exp(\lambda^i_c t) / \sum_c k_c exp(\lambda^i_c t)
//!
use crate::gtr::gtr::GTR;
use crate::partition::payload::dense::DenseSeqDis;
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

pub fn get_coefficients(msg_to_parent: &DenseSeqDis, msg_to_child: &DenseSeqDis, gtr: &GTR) -> PartitionContribution {
  // Multiply the messages by the eigenvectors of the GTR matrix, multiply elementwise, and sum over the rows:
  //    s_a eQt_{ab} r_b =  \sum_{abc} s_a v_{ac} e^{\lambda_c t} vinv_{cb} r_b
  let coefficients = msg_to_child.dis.dot(&gtr.v) * msg_to_parent.dis.dot(&gtr.v_inv.t());
  PartitionContribution::new(
    coefficients,
    gtr.clone(), // TODO: avoid clone
  )
}

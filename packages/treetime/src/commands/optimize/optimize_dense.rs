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
use crate::commands::optimize::optimize_unified::OptimizationMetrics;
use crate::graph::edge::HasBranchLength;
use crate::gtr::gtr::GTR;
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::graph_dense::DenseSeqDis;
use crate::representation::partition_marginal_dense::PartitionMarginalDense;
use eyre::Report;
use ndarray::{Array2, Axis};
use num::clamp;
use parking_lot::RwLock;
use std::sync::Arc;

pub struct PartitionContribution {
  pub coefficients: Array2<f64>,
  pub gtr: GTR,
}

impl PartitionContribution {
  pub fn new(coefficients: Array2<f64>, gtr: GTR) -> Self {
    Self { coefficients, gtr }
  }
}

pub fn evaluate(contributions: &[PartitionContribution], branch_length: f64) -> OptimizationMetrics {
  let mut log_likelihood = 0.0;
  let mut derivative = 0.0;
  let mut second_derivative = 0.0;
  for contribution in contributions {
    let gtr = &contribution.gtr;
    let exp_ev = gtr.eigvals.mapv(|ev| (ev * branch_length).exp());
    let ev_exp_ev = &gtr.eigvals * &exp_ev;
    let ev2_exp_ev = &gtr.eigvals * &ev_exp_ev;
    // This loop could be coded more efficiently
    for coeff in contribution.coefficients.outer_iter() {
      let val = (&coeff * &exp_ev).sum();
      log_likelihood += val.ln();
      derivative += (&coeff * &ev_exp_ev).sum() / val;
      second_derivative += (&coeff * &ev2_exp_ev).sum() / val - ((&coeff * &ev_exp_ev).sum() / val).powi(2);
    }
  }
  OptimizationMetrics::new(log_likelihood, derivative, second_derivative)
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

pub fn run_optimize_dense(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalDense>>],
) -> Result<(), Report> {
  let total_length: usize = partitions
    .iter()
    .map(|part| part.read_arc().get_sequence_length())
    .sum();
  let one_mutation = 1.0 / total_length as f64;
  let n_partitions = partitions.len();
  graph.get_edges().iter().for_each(|edge_ref| {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let mut branch_length = edge.branch_length().unwrap_or(0.0);
    let mut new_branch_length;
    let contributions = (0..n_partitions)
      .map(|pi| {
        let partition = partitions[pi].read_arc();
        let edge_partition = &partition.edges[&edge_key];
        get_coefficients(
          &edge_partition.msg_to_parent,
          &edge_partition.msg_to_child,
          &partition.gtr,
        )
      })
      .collect::<Vec<_>>();

    // Cheap check whether the branch length is zero
    let zero_branch_length_lh: f64 = contributions
      .iter()
      .map(|contribution| contribution.coefficients.sum_axis(Axis(1)).product())
      .product();

    if zero_branch_length_lh > 0.01 {
      let zero_branch_length_derivative: f64 = contributions
        .iter()
        .map(|contribution| {
          let gtr = &contribution.gtr;
          ((&contribution.coefficients * &gtr.eigvals).sum_axis(Axis(1)) / &contribution.coefficients.sum_axis(Axis(1)))
            .sum()
        })
        .sum();
      if zero_branch_length_derivative < 0.0 {
        edge.set_branch_length(Some(0.0));
        return;
      }
    }

    // Otherwise, we need to optimize the branch length
    let metrics = evaluate(&contributions, branch_length);
    if metrics.log_lh.is_finite() && metrics.second_derivative < 0.0 {
      // Newton's method to find the optimal branch length
      new_branch_length = branch_length - clamp(metrics.derivative / metrics.second_derivative, -1.0, branch_length);
      let max_iter = 10;
      let mut n_iter = 0;
      while (new_branch_length - branch_length).abs() > 0.001 * branch_length && n_iter < max_iter {
        let new_metrics = evaluate(&contributions, new_branch_length);
        if new_metrics.second_derivative < 0.0 {
          branch_length = new_branch_length;
          new_branch_length = branch_length
            - clamp(
              new_metrics.derivative / new_metrics.second_derivative,
              -1.0,
              branch_length,
            );
        } else {
          break;
        }
        n_iter += 1;
      }
    } else {
      // Evaluate on a vector of branch lengths to find the maximum
      let branch_lengths = ndarray::Array1::linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 10);
      // This seems like a bit of a mess.
      let (best_branch_length, _) = branch_lengths
        .iter()
        .map(|&bl| {
          let metrics = evaluate(&contributions, bl);
          (bl, metrics.log_lh)
        })
        .max_by(|&(_, ll1), &(_, ll2)| ll1.partial_cmp(&ll2).unwrap())
        .unwrap();
      new_branch_length = best_branch_length;
    }
    edge.set_branch_length(Some(new_branch_length));
  });
  Ok(())
}

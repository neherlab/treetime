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
use crate::{
  gtr::gtr::GTR,
  representation::{
    graph_dense::{DenseGraph, DenseSeqDis},
    partitions_likelihood::PartitionLikelihood,
  },
};
use eyre::Report;
use ndarray::{Array2, Axis};
use ndarray_stats::QuantileExt;
use num::clamp;
use std::iter::zip;

pub fn evaluate(
  coefficients: &[Array2<f64>],
  partitions: &[PartitionLikelihood],
  branch_length: f64,
) -> (f64, f64, f64, f64) {
  let mut likelihood = 1.0;
  let mut log_likelihood = 0.0;
  let mut derivative = 0.0;
  let mut second_derivative = 0.0;
  for (pi, partition) in partitions.iter().enumerate() {
    let PartitionLikelihood { gtr, .. } = &partition;
    // let coefficients = &partition.msg_to_parent.dis * &partition.msg_to_child.dis;
    let exp_ev = gtr.eigvals.mapv(|ev| (ev * branch_length).exp());
    let ev_exp_ev = &gtr.eigvals * &exp_ev;
    let ev2_exp_ev = &gtr.eigvals * &ev_exp_ev;
    // this loop could be coded more efficiently
    for coeff in coefficients[pi].outer_iter() {
      let val = (&coeff * &exp_ev).sum();
      likelihood *= val;
      log_likelihood += val.ln();
      derivative += (&coeff * &ev_exp_ev).sum() / val;
      second_derivative += (&coeff * &ev2_exp_ev).sum() / val - ((&coeff * &ev_exp_ev).sum() / val).powi(2);
    }
  }
  (likelihood, log_likelihood, derivative, second_derivative)
}

pub fn get_coefficients(msg_to_parent: &DenseSeqDis, msg_to_child: &DenseSeqDis, gtr: &GTR) -> Array2<f64> {
  // multiply the messages by the eigenvectors of the GTR matrix, multiply elementwise, and sum over the rows
  // s_a eQt_{ab} r_b =  \sum_{abc} s_a v_{ac} e^{\lambda_c t} vinv_{cb} r_b
  msg_to_child.dis.dot(&gtr.v) * msg_to_parent.dis.dot(&gtr.v_inv.t())
}

pub fn initial_guess(graph: &DenseGraph, partitions: &[PartitionLikelihood]) -> () {
  //FIXME: this initial guess needs to be improved
  let total_length: usize = partitions.iter().map(|part| part.length).sum();
  let one_mutation = 1.0 / total_length as f64;
  for edge in graph.get_edges() {
    let mut edge = edge.write_arc().payload().write_arc();
    let mut differences: usize = 0;
    for (_, partition) in edge.dense_partitions.iter().enumerate() {
      for (row1, row2) in zip(partition.msg_to_parent.dis.rows(), partition.msg_to_child.dis.rows()) {
        if row1[row2.argmax().unwrap()] < 0.5 {
          differences += 1;
        };
      }
    }
    let new_branch_length = differences as f64 * one_mutation;
    edge.branch_length = Some(new_branch_length);
  }
}

pub fn run_optimize_dense(graph: &DenseGraph, partitions: &[PartitionLikelihood]) -> Result<(), Report> {
  let total_length: usize = partitions.iter().map(|part| part.length).sum();
  let one_mutation = 1.0 / total_length as f64;
  let n_partitions = partitions.len();
  graph.get_edges().iter_mut().for_each(|edge| {
    let mut edge = edge.write_arc().payload().write_arc();
    let mut branch_length = edge.branch_length.unwrap_or(0.0);
    let mut new_branch_length;
    let coefficients = (0..n_partitions)
      .map(|pi| {
        get_coefficients(
          &edge.dense_partitions[pi].msg_to_parent,
          &edge.dense_partitions[pi].msg_to_child,
          &partitions[pi].gtr,
        )
      })
      .collect::<Vec<_>>();

    // cheap check whether the branch length is zero
    let zero_branch_length_lh: f64 = coefficients
      .iter()
      .map(|coeff| coeff.sum_axis(Axis(1)).product())
      .product();

    if zero_branch_length_lh > 0.01 {
      let zero_branch_length_derivative: f64 = coefficients
        .iter()
        .enumerate()
        .map(|(pi, coeff)| ((coeff * &partitions[pi].gtr.eigvals).sum_axis(Axis(1)) / coeff.sum_axis(Axis(1))).sum())
        .sum();
      if zero_branch_length_derivative < 0.0 {
        edge.branch_length = Some(0.0);
        return;
      }
    }

    // otherwise, we need to optimize the branch length
    let (likelihood, _, derivative, second_derivative) = evaluate(&coefficients, partitions, branch_length);
    if likelihood > 0.0 && second_derivative < 0.0 {
      // newton's method to find the optimal branch length
      new_branch_length = branch_length - clamp(derivative / second_derivative, -1.0, branch_length);
      let max_iter = 10;
      let mut n_iter = 0;
      while (new_branch_length - branch_length).abs() > 0.001 * branch_length && n_iter < max_iter {
        let (_, _, derivative, second_derivative) = evaluate(&coefficients, partitions, new_branch_length);
        if second_derivative < 0.0 {
          branch_length = new_branch_length;
          new_branch_length = branch_length - clamp(derivative / second_derivative, -1.0, branch_length);
        } else {
          break;
        }
        n_iter += 1;
      }
    } else {
      // evaluate on a vector of branch lengths to find the maximum
      let branch_lengths = ndarray::Array1::linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 10);
      // this seems like a bit of a mess.
      let (best_branch_length, _) = branch_lengths
        .iter()
        .map(|&bl| {
          let (_, ll, _, _) = evaluate(&coefficients, partitions, bl);
          (bl, ll)
        })
        .max_by(|&(_, ll1), &(_, ll2)| ll1.partial_cmp(&ll2).unwrap())
        .unwrap();
      new_branch_length = best_branch_length;
    };
    edge.branch_length = Some(new_branch_length);
  });
  Ok(())
}

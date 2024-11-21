use crate::{
  gtr::gtr::GTR,
  representation::{
    graph_sparse::{SparseGraph, SparseSeqEdge},
    partitions_likelihood::PartitionLikelihood,
  },
};
use eyre::Report;
use itertools::Itertools;
use num::clamp;
use std::iter::zip;

// The likelihood of an edge length is the product of the likelihoods of all positions of all partitions
// Lh = prod_i prod_j \sum_{ab} s^{ij}_a exp(Q_i t)_{ab} r^{ij}_b
// The log likelihood is the sum of many terms
// logLh = sum_i sum_j \log(\sum_{ab} s^{ij}_a exp(Q_i t)_{ab} r^{ij}_b)
// to effectively calculate this, we need to reformulate the likelihood in terms of the eigenvectors of the GTR matrix
// Dropping the {ij} superscripts for brevity, we can write the likelihood as
// s_a exp(Qt)_{ab} r_b = s_a \sum_c v_{ac} exp(\lambda_c t) vinv_{cb} r_b = g_c exp(\lambda_c t) h_c = k_c exp(\lambda_c t)
// the k_c can be reused for different iterations of the branch length optimization
// logLh = sum_i sum_j \log(\sum_c k_c exp(\lambda^i_c t))
// the derivative is simply
// dlogLh/dt = sum_i sum_j \sum_c k_c \lambda_c exp(\lambda^i_c t) / \sum_c k_c exp(\lambda^i_c t)

struct SiteContribution {
  multiplicity: f64,
  coefficients: ndarray::Array1<f64>,
}

struct PartitionContribution {
  site_contributions: Vec<SiteContribution>,
  eigenvalues: ndarray::Array1<f64>,
}

fn get_coefficients(edge: &SparseSeqEdge, gtr: &GTR) -> PartitionContribution {
  let mut site_contributions: Vec<SiteContribution> = Vec::new();
  // collect variable positions from msg_to_child, msg_to_parent, and the substitutions along the edge
  let variable_positions: Vec<usize> = edge
    .msg_to_child
    .variable
    .keys()
    .copied()
    .chain(edge.msg_to_parent.variable.keys().copied())
    .chain(edge.subs.iter().map(|sub| sub.pos))
    .unique()
    .collect();

  let variable_states = variable_positions.iter().map(|pos| {
    // check whether the position is in substitutions
    let state_pair = if let Some(sub) = edge.subs.iter().find(|m| m.pos == *pos) {
      (sub.reff, sub.qry)
    } else {
      let parent = edge
        .msg_to_child
        .variable
        .get(pos)
        .or_else(|| edge.msg_to_parent.variable.get(pos))
        .unwrap()
        .state;
      let child = edge
        .msg_to_parent
        .variable
        .get(pos)
        .or_else(|| edge.msg_to_child.variable.get(pos))
        .unwrap()
        .state;
      (parent, child)
    };
    state_pair
  });

  for (&pos, (parent_state, child_state)) in zip(&variable_positions, variable_states) {
    let parent = if let Some(parent) = edge.msg_to_child.variable.get(&pos) {
      &parent.dis
    } else {
      &edge.msg_to_child.fixed[&parent_state]
    };

    let child = if let Some(child) = edge.msg_to_parent.variable.get(&pos) {
      &child.dis
    } else {
      &edge.msg_to_parent.fixed[&child_state]
    };
    site_contributions.push(SiteContribution {
      multiplicity: 1.0,
      coefficients: (parent.dot(&gtr.v) * child.dot(&gtr.v_inv.t())),
    });
  }
  for state in edge.msg_to_child.fixed.keys() {
    let parent = &edge.msg_to_child.fixed[state];
    let child = &edge.msg_to_parent.fixed[state];
    site_contributions.push(SiteContribution {
      multiplicity: edge.msg_to_child.fixed_counts.get(*state).unwrap() as f64,
      coefficients: (parent.dot(&gtr.v) * child.dot(&gtr.v_inv.t())),
    });
  }
  PartitionContribution {
    site_contributions,
    eigenvalues: gtr.eigvals.to_owned(),
  }
}

// function that takes two message projections, and gtr model, and the length of branch and returns the
// likelihood as well as its derivative with respect to the branch length
fn evaluate_sparse(coefficients: &Vec<PartitionContribution>, branch_length: f64) -> (f64, f64, f64, f64) {
  let likelihood = 0.0; //currently not calculate
  let mut log_likelihood = 0.0;
  let mut derivative = 0.0;
  let mut second_derivative = 0.0;
  for coeffs in coefficients {
    let exp_ev = coeffs.eigenvalues.mapv(|ev| (ev * branch_length).exp());
    let ev_exp_ev = &coeffs.eigenvalues * &exp_ev;
    let ev2_exp_ev = &coeffs.eigenvalues * &ev_exp_ev;
    for coeff in &coeffs.site_contributions {
      let val = (&coeff.coefficients * &exp_ev).sum();
      log_likelihood += coeff.multiplicity * val.ln();
      derivative += coeff.multiplicity * (&coeff.coefficients * &ev_exp_ev).sum() / val;
      second_derivative += coeff.multiplicity * (&coeff.coefficients * &ev2_exp_ev).sum() / val
        - (coeff.multiplicity * (&coeff.coefficients * &ev_exp_ev).sum() / val).powi(2);
    }
  }
  (likelihood, log_likelihood, derivative, second_derivative)
}

pub fn initial_guess_sparse(graph: &SparseGraph, partitions: &[PartitionLikelihood]) -> () {
  let total_length: usize = partitions.iter().map(|part| part.length).sum();
  let one_mutation = 1.0 / total_length as f64;
  for edge in graph.get_edges() {
    let mut edge = edge.write_arc().payload().write_arc();
    let mut differences: usize = 0;
    for partition in &edge.sparse_partitions {
      differences += partition.subs.len();
    }
    let new_branch_length = differences as f64 * one_mutation;
    edge.branch_length = Some(new_branch_length);
  }
}

pub fn run_optimize_sparse(graph: &SparseGraph, partitions: &[PartitionLikelihood]) -> Result<(), Report> {
  let total_length: usize = partitions.iter().map(|part| part.length).sum();
  let one_mutation = 1.0 / total_length as f64;
  let n_partitions = partitions.len();
  graph.get_edges().iter_mut().for_each(|edge| {
    let name = &graph
      .get_node(edge.read_arc().target())
      .unwrap()
      .read_arc()
      .payload()
      .read_arc()
      .name
      .clone();
    let mut edge = edge.write_arc().payload().write_arc();
    let coefficients = (0..n_partitions)
      .map(|pi| get_coefficients(&edge.sparse_partitions[pi], &partitions[pi].gtr))
      .collect_vec();
    let mut branch_length = edge.branch_length.unwrap_or(0.0);
    let mut new_branch_length;

    let zero_branch_length_lh: f64 = coefficients
      .iter()
      .map(|coeffs| {
        coeffs
          .site_contributions
          .iter()
          .map(|coeff| coeff.coefficients.sum().powf(coeff.multiplicity))
          .product::<f64>()
      })
      .product();

    if zero_branch_length_lh > 0.0001 {
      // TODO: could check that derivative is negative
      edge.branch_length = Some(0.0);
      return;
    }

    // otherwise, we need to optimize the branch length
    let (likelihood, log_likelihood, derivative, second_derivative) = evaluate_sparse(&coefficients, branch_length);
    if second_derivative < 0.0 {
      // newton's method to find the optimal branch length
      new_branch_length = branch_length - clamp(derivative / second_derivative, -1.0, branch_length);
      let max_iter = 10;
      let mut n_iter = 0;
      while (new_branch_length - branch_length).abs() > 0.001 * branch_length && n_iter < max_iter {
        let (likelihood, log_likelihood, derivative, second_derivative) =
          evaluate_sparse(&coefficients, new_branch_length);
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
      let branch_lengths = ndarray::Array1::linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 100);

      // this seems like a bit of a mess.
      let (best_branch_length, _) = branch_lengths
        .iter()
        .map(|&bl| {
          let (_, ll, _, _) = evaluate_sparse(&coefficients, bl);
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

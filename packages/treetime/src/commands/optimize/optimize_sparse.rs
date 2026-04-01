//! The likelihood of an edge length is the product of the likelihoods of all positions of all partitions:
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
//! The `k_c` can be reused for different iterations of the branch length optimization:
//!
//!   logLh = sum_i sum_j \log(\sum_c k_c exp(\lambda^i_c t))
//!
//! For compressed site patterns with multiplicity $m_i$, the derivatives are:
//!
//!   dlogLh/dt = sum_i m_i * (sum_c k_c \lambda_c exp(\lambda_c t)) / (sum_c k_c exp(\lambda_c t))
//!
//!   d^2logLh/dt^2 = sum_i m_i * [(sum_c k_c \lambda_c^2 exp(\lambda_c t)) / L_i
//!                                 - ((sum_c k_c \lambda_c exp(\lambda_c t)) / L_i)^2]
//!
//! where L_i = sum_c k_c exp(\lambda_c t). The multiplicity is a linear factor on
//! each site's contribution; the squared term applies only to the per-site ratio.
//!
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::seq::mutation::Sub;
use eyre::{OptionExt, Report};
use itertools::Itertools;
use std::iter::zip;
use treetime_graph::edge::GraphEdgeKey;

pub struct SiteContribution {
  pub multiplicity: f64,
  pub coefficients: ndarray::Array1<f64>,
}

pub struct PartitionContribution {
  pub site_contributions: Vec<SiteContribution>,
  pub eigenvalues: ndarray::Array1<f64>,
  pub unimodal_branch_likelihood: bool,
}

pub fn get_coefficients(
  edge_key: GraphEdgeKey,
  partition: &PartitionMarginalSparse,
) -> Result<PartitionContribution, Report> {
  let edge = &partition.edges[&edge_key];

  // Collect variable positions from msg_to_child, msg_to_parent, and the substitutions along the edge
  let variable_positions: Vec<usize> = edge
    .msg_to_child
    .variable
    .keys()
    .copied()
    .chain(edge.msg_to_parent.variable.keys().copied())
    .chain(edge.subs.iter().map(Sub::pos))
    .unique()
    .collect();

  let variable_states = variable_positions
    .iter()
    .map(|pos| -> Result<_, Report> {
      // Check whether the position is in substitutions
      if let Some(sub) = edge.subs.iter().find(|m| m.pos() == *pos) {
        Ok((sub.reff(), sub.qry()))
      } else {
        let parent = edge
          .msg_to_child
          .variable
          .get(pos)
          .or_else(|| edge.msg_to_parent.variable.get(pos))
          .ok_or_eyre("Unable to find msg_to_parent")?
          .state;
        let child = edge
          .msg_to_parent
          .variable
          .get(pos)
          .or_else(|| edge.msg_to_child.variable.get(pos))
          .ok_or_eyre("Unable to find msg_to_child")?
          .state;
        Ok((parent, child))
      }
    })
    .collect::<Result<Vec<_>, Report>>()?;

  let mut site_contributions: Vec<SiteContribution> = Vec::new();
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
      coefficients: parent.dot(&partition.gtr.v) * child.dot(&partition.gtr.v_inv.t()),
    });
  }
  for state in edge.msg_to_child.fixed.keys() {
    let parent = &edge.msg_to_child.fixed[state];
    let child = &edge.msg_to_parent.fixed[state];
    site_contributions.push(SiteContribution {
      multiplicity: edge.msg_to_child.fixed_counts.get(*state).unwrap() as f64,
      coefficients: parent.dot(&partition.gtr.v) * child.dot(&partition.gtr.v_inv.t()),
    });
  }
  Ok(PartitionContribution {
    site_contributions,
    eigenvalues: partition.gtr.eigvals.to_owned(),
    unimodal_branch_likelihood: partition.gtr.unimodal_branch_likelihood,
  })
}

use crate::gtr::gtr::GTR;
use crate::partition::storage::sparse::{SparseEdgeBackward, SparseEdgeForward, SparseEdgeObs};
use crate::seq::mutation::Sub;
use eyre::{OptionExt, Report};
use itertools::Itertools;
use std::iter::zip;

#[allow(
  clippy::as_conversions,
  clippy::unwrap_used,
  reason = "count/index numeric cast is exact for the domain range; unwrap on a value an upstream invariant guarantees is present"
)]
pub(crate) fn get_coefficients(
  gtr: &GTR,
  backward: &SparseEdgeBackward,
  forward: &SparseEdgeForward,
  edge_obs: &SparseEdgeObs,
) -> Result<PartitionContribution, Report> {
  let msg_to_child = &forward.msg_to_child;
  let msg_to_parent = &backward.msg_to_parent;

  let variable_positions: Vec<usize> = msg_to_child
    .variable
    .keys()
    .copied()
    .chain(msg_to_parent.variable.keys().copied())
    .chain(edge_obs.fitch_subs().iter().map(Sub::pos))
    .unique()
    .collect();

  let variable_states = variable_positions
    .iter()
    .map(|pos| -> Result<_, Report> {
      if let Some(sub) = edge_obs.fitch_subs().iter().find(|m| m.pos() == *pos) {
        Ok((sub.reff(), sub.qry()))
      } else {
        let parent = msg_to_child
          .variable
          .get(pos)
          .or_else(|| msg_to_parent.variable.get(pos))
          .ok_or_eyre("Unable to find msg_to_parent")?
          .state;
        let child = msg_to_parent
          .variable
          .get(pos)
          .or_else(|| msg_to_child.variable.get(pos))
          .ok_or_eyre("Unable to find msg_to_child")?
          .state;
        Ok((parent, child))
      }
    })
    .collect::<Result<Vec<_>, Report>>()?;

  let mut site_contributions: Vec<SiteContribution> = Vec::new();
  for (&pos, (parent_state, child_state)) in zip(&variable_positions, variable_states) {
    let parent = if let Some(parent) = msg_to_child.variable.get(&pos) {
      &parent.dis
    } else {
      &msg_to_child.fixed[&parent_state]
    };

    let child = if let Some(child) = msg_to_parent.variable.get(&pos) {
      &child.dis
    } else {
      &msg_to_parent.fixed[&child_state]
    };
    site_contributions.push(SiteContribution {
      multiplicity: 1.0,
      coefficients: parent.dot(&gtr.v) * child.dot(&gtr.v_inv.t()),
    });
  }
  for state in msg_to_child.fixed.keys() {
    let parent = &msg_to_child.fixed[state];
    let child = &msg_to_parent.fixed[state];
    site_contributions.push(SiteContribution {
      multiplicity: msg_to_child.fixed_counts.get(*state).unwrap() as f64,
      coefficients: parent.dot(&gtr.v) * child.dot(&gtr.v_inv.t()),
    });
  }
  Ok(PartitionContribution {
    site_contributions,
    gtr: gtr.clone(),
  })
}

pub struct PartitionContribution {
  pub(crate) site_contributions: Vec<SiteContribution>,
  pub(crate) gtr: GTR,
}

pub struct SiteContribution {
  pub(crate) multiplicity: f64,
  pub(crate) coefficients: ndarray::Array1<f64>,
}

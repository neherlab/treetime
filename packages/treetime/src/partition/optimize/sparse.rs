use crate::gtr::gtr::GTR;
use crate::make_internal_report;
use crate::partition::storage::sparse::{SparseEdgeBackward, SparseEdgeForward, SparseEdgeObs};
use crate::seq::mutation::Sub;
use eyre::{OptionExt, Report};
use itertools::Itertools;
use ndarray::{Array1, aview1};
use std::collections::BTreeMap;
use std::iter::zip;
use treetime_utils::cast::ToF64Exact;

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

  let mut fitch_states = BTreeMap::new();
  for sub in edge_obs.fitch_subs() {
    fitch_states.entry(sub.pos()).or_insert_with(|| (sub.reff(), sub.qry()));
  }

  let variable_states = variable_positions
    .iter()
    .map(|pos| -> Result<_, Report> {
      if let Some(&states) = fitch_states.get(pos) {
        Ok(states)
      } else {
        let parent = msg_to_child
          .variable
          .get(*pos)
          .or_else(|| msg_to_parent.variable.get(*pos))
          .ok_or_eyre("Unable to find msg_to_parent")?
          .state;
        let child = msg_to_parent
          .variable
          .get(*pos)
          .or_else(|| msg_to_child.variable.get(*pos))
          .ok_or_eyre("Unable to find msg_to_child")?
          .state;
        Ok((parent, child))
      }
    })
    .collect::<Result<Vec<_>, Report>>()?;

  let mut site_contributions = Vec::with_capacity(variable_positions.len() + msg_to_child.fixed.len());
  for (&pos, (parent_state, child_state)) in zip(&variable_positions, variable_states) {
    let parent = msg_to_child.profile_at(pos, parent_state).ok_or_else(|| {
      make_internal_report!("Outside message has no profile of state {parent_state} at position {pos}")
    })?;
    let child = msg_to_parent
      .profile_at(pos, child_state)
      .ok_or_else(|| make_internal_report!("Inside message has no profile of state {child_state} at position {pos}"))?;
    site_contributions.push(SiteContribution {
      multiplicity: 1.0,
      coefficients: aview1(child).dot(gtr.v()) * aview1(parent).dot(&gtr.v_inv().t()),
    });
  }
  for group in 0..msg_to_child.groups.len() {
    let child_group = msg_to_child
      .groups
      .first_position(group)
      .and_then(|pos| msg_to_parent.groups.group_at(pos));
    for parent_entry in msg_to_child.groups.profiles(group).iter() {
      let state = parent_entry.state;
      let child = child_group
        .and_then(|child_group| msg_to_parent.groups.profiles(child_group).get(state))
        .or_else(|| msg_to_parent.fixed.get(state))
        .ok_or_else(|| make_internal_report!("Inside message has no profile of grouped state {state}"))?;
      site_contributions.push(SiteContribution {
        multiplicity: parent_entry.count.to_f64_exact(),
        coefficients: aview1(&child.dis).dot(gtr.v()) * aview1(&parent_entry.dis).dot(&gtr.v_inv().t()),
      });
    }
  }
  for parent_fixed in msg_to_child.fixed.iter() {
    let parent = &parent_fixed.dis;
    let child = &msg_to_parent.fixed[parent_fixed.state].dis;
    site_contributions.push(SiteContribution {
      multiplicity: parent_fixed.count.to_f64_exact(),
      coefficients: aview1(child).dot(gtr.v()) * aview1(parent).dot(&gtr.v_inv().t()),
    });
  }
  Ok(PartitionContribution::new(site_contributions, gtr.clone()))
}

pub struct PartitionContribution {
  pub(crate) site_contributions: Vec<SiteContribution>,
  pub(crate) rate_eigvals: Array1<f64>,
  pub(crate) gtr: GTR,
}

impl PartitionContribution {
  pub(crate) fn new(site_contributions: Vec<SiteContribution>, gtr: GTR) -> Self {
    Self {
      site_contributions,
      rate_eigvals: gtr.rate_eigvals(),
      gtr,
    }
  }
}

pub struct SiteContribution {
  pub(crate) multiplicity: f64,
  pub(crate) coefficients: Array1<f64>,
}

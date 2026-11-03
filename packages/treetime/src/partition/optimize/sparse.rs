use crate::gtr::gtr::GTR;
use crate::partition::storage::sparse::INLINE_STATES;
use crate::partition::storage::sparse::{SparseEdgeBackward, SparseEdgeForward, SparseEdgeObs};
use crate::seq::mutation::Sub;
use eyre::{OptionExt, Report};
use indexmap::IndexMap;
use itertools::Itertools;
use ndarray::{Array1, aview1};
use smallvec::SmallVec;
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

  let mut site_contributions =
    GroupedSiteContributions::with_capacity(variable_positions.len() + msg_to_child.fixed.len());
  for (&pos, (parent_state, child_state)) in zip(&variable_positions, variable_states) {
    let parent = if let Some(parent) = msg_to_child.variable.get(pos) {
      parent.dis
    } else {
      msg_to_child.fixed[parent_state].dis.as_slice()
    };

    let child = if let Some(child) = msg_to_parent.variable.get(pos) {
      child.dis
    } else {
      msg_to_parent.fixed[child_state].dis.as_slice()
    };
    site_contributions.add(1.0, aview1(child).dot(gtr.v()) * aview1(parent).dot(&gtr.v_inv().t()));
  }
  for parent_fixed in msg_to_child.fixed.iter() {
    let parent = &parent_fixed.dis;
    let child = &msg_to_parent.fixed[parent_fixed.state].dis;
    site_contributions.add(
      parent_fixed.count.to_f64_exact(),
      aview1(child).dot(gtr.v()) * aview1(parent).dot(&gtr.v_inv().t()),
    );
  }
  Ok(PartitionContribution::new(site_contributions.into_vec(), gtr.clone()))
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

struct GroupedSiteContributions {
  sites: IndexMap<SmallVec<[u64; INLINE_STATES]>, SiteContribution>,
}

impl GroupedSiteContributions {
  fn with_capacity(capacity: usize) -> Self {
    Self {
      sites: IndexMap::with_capacity(capacity),
    }
  }

  fn add(&mut self, multiplicity: f64, coefficients: Array1<f64>) {
    let bits = coefficients.iter().map(|value| value.to_bits()).collect();
    self
      .sites
      .entry(bits)
      .and_modify(|site| site.multiplicity += multiplicity)
      .or_insert(SiteContribution {
        multiplicity,
        coefficients,
      });
  }

  fn into_vec(self) -> Vec<SiteContribution> {
    self.sites.into_values().collect()
  }
}

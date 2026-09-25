use crate::results::tree::{ResultTree, preorder};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use treetime_io::auspice_types::AuspiceTree;

pub const UNCERTAIN_STATE_PROBABILITY: f64 = 0.8;

/// Results of a `mugration` run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct MugrationResults {
  /// The reconstructed attribute.
  pub attribute: String,
  /// Number of distinct states among the samples.
  pub states: usize,
  /// Changes of state from parent to child, counted over branches, most frequent first.
  pub state_changes: Vec<StateChange>,
  /// Number of branches whose state differs from the parent's.
  pub changed_branches: usize,
  /// Probability below which an ancestor's most probable state counts as uncertain.
  pub uncertain_below: f64,
  /// Ancestors whose most probable state is uncertain, least certain first.
  pub uncertain_ancestors: Vec<AncestorState>,
  /// Most probable state of the root.
  pub root: Option<AncestorState>,
}

/// A change of state along branches.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct StateChange {
  /// State of the parent.
  pub from: String,
  /// State of the child.
  pub to: String,
  /// Number of branches with this change.
  pub branches: usize,
}

/// Most probable state of an ancestor.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct AncestorState {
  /// Name of the ancestor.
  pub name: String,
  /// Number of samples below the ancestor.
  pub tips: usize,
  /// Most probable state.
  pub state: String,
  /// Probability of the state.
  pub probability: f64,
}

pub fn mugration_results(
  auspice: Option<&AuspiceTree>,
  tree: Option<&ResultTree>,
  attribute: &str,
) -> Result<MugrationResults, Report> {
  let (Some(auspice), Some(tree)) = (auspice, tree) else {
    return Ok(MugrationResults {
      attribute: attribute.to_owned(),
      states: 0,
      state_changes: vec![],
      changed_branches: 0,
      uncertain_below: UNCERTAIN_STATE_PROBABILITY,
      uncertain_ancestors: vec![],
      root: None,
    });
  };
  let traits = preorder(&auspice.tree)
    .into_iter()
    .map(|(node, _)| node.node_attrs.attr(attribute))
    .collect::<Result<Vec<_>, Report>>()?;
  let values = traits
    .iter()
    .map(|attr| attr.as_ref().map(|attr| attr.value().to_owned()))
    .collect_vec();

  let mut counts: BTreeMap<(String, String), usize> = BTreeMap::new();
  for (index, node) in tree.nodes.iter().enumerate() {
    if let (Some(parent), Some(to)) = (node.parent, &values[index])
      && let Some(from) = &values[parent]
      && from != to
    {
      *counts.entry((from.clone(), to.clone())).or_default() += 1;
    }
  }
  let state_changes = counts
    .into_iter()
    .map(|((from, to), branches)| StateChange { from, to, branches })
    .sorted_by(|a, b| b.branches.cmp(&a.branches))
    .collect_vec();

  let mut ancestors = vec![];
  for (index, node) in tree.nodes.iter().enumerate() {
    if let Some(attr) = traits[index].as_ref().filter(|_| !node.is_tip()) {
      if let Some(probability) = attr
        .state_confidence()?
        .and_then(|confidence| confidence.get(attr.value()).copied())
      {
        ancestors.push((
          index,
          AncestorState {
            name: node.name.clone(),
            tips: node.tips,
            state: attr.value().to_owned(),
            probability,
          },
        ));
      }
    }
  }

  let states: BTreeSet<&String> = tree
    .nodes
    .iter()
    .zip(&values)
    .filter(|(node, _)| node.is_tip())
    .filter_map(|(_, value)| value.as_ref())
    .collect();

  Ok(MugrationResults {
    attribute: attribute.to_owned(),
    states: states.len(),
    changed_branches: state_changes.iter().map(|change| change.branches).sum(),
    state_changes,
    uncertain_below: UNCERTAIN_STATE_PROBABILITY,
    root: ancestors
      .iter()
      .find(|(index, _)| *index == 0)
      .map(|(_, state)| state.clone()),
    uncertain_ancestors: ancestors
      .into_iter()
      .map(|(_, state)| state)
      .filter(|state| state.probability < UNCERTAIN_STATE_PROBABILITY)
      .sorted_by(|a, b| a.probability.total_cmp(&b.probability))
      .collect(),
  })
}

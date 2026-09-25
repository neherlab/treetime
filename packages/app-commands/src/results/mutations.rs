use crate::results::tree::ResultTree;
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use treetime_utils::make_report;

/// Results of an `ancestral` run.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct AncestralResults {
  /// Number of nucleotide mutations on all branches.
  pub mutations: usize,
  /// Branches with at least one mutation, most mutations first.
  pub branches: Vec<BranchMutations>,
  /// Sequence positions that mutate on more than one branch, most branches first.
  pub recurrent_sites: Vec<RecurrentSite>,
}

/// Mutations on the branch above one node.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct BranchMutations {
  /// Name of the node below the branch.
  pub name: String,
  /// Number of samples below the branch.
  pub tips: usize,
  /// Mutations on the branch.
  pub mutations: Vec<String>,
}

/// A sequence position that mutates on several branches.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct RecurrentSite {
  /// Position in the sequence, from 1.
  pub position: usize,
  /// Number of branches with a mutation at the position.
  pub branches: usize,
}

pub fn ancestral_results(tree: Option<&ResultTree>) -> Result<AncestralResults, Report> {
  let Some(tree) = tree else {
    return Ok(AncestralResults {
      mutations: 0,
      branches: vec![],
      recurrent_sites: vec![],
    });
  };
  let branches = tree
    .nodes
    .iter()
    .filter(|node| !node.mutations.is_empty())
    .map(|node| BranchMutations {
      name: node.name.clone(),
      tips: node.tips,
      mutations: node.mutations.clone(),
    })
    .sorted_by(|a, b| {
      b.mutations
        .len()
        .cmp(&a.mutations.len())
        .then_with(|| a.name.cmp(&b.name))
    })
    .collect_vec();
  let mut branches_per_site: BTreeMap<usize, usize> = BTreeMap::new();
  for node in &tree.nodes {
    let positions = node
      .mutations
      .iter()
      .map(|mutation| mutation_position(mutation))
      .collect::<Result<BTreeSet<_>, Report>>()?;
    for position in positions {
      *branches_per_site.entry(position).or_default() += 1;
    }
  }
  let recurrent_sites = branches_per_site
    .into_iter()
    .filter(|&(_, branches)| branches > 1)
    .map(|(position, branches)| RecurrentSite { position, branches })
    .sorted_by(|a, b| b.branches.cmp(&a.branches).then_with(|| a.position.cmp(&b.position)))
    .collect();
  Ok(AncestralResults {
    mutations: branches.iter().map(|branch| branch.mutations.len()).sum(),
    branches,
    recurrent_sites,
  })
}

fn mutation_position(mutation: &str) -> Result<usize, Report> {
  mutation
    .split(|c: char| !c.is_ascii_digit())
    .filter(|digits| !digits.is_empty())
    .exactly_one()
    .map_err(|groups| {
      make_report!(
        "mutation `{mutation}` names {} sequence positions, expected one",
        groups.count()
      )
    })?
    .parse()
    .map_err(|err| make_report!("mutation `{mutation}` has an invalid position: {err}"))
}

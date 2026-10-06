use crate::results::tree::ResultTree;
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::str::FromStr;
use treetime::alphabet::alphabet::Alphabet;
use treetime::homoplasy::classify::MutationClass;
use treetime::homoplasy::site_branches::sites_by_branch_count;
use treetime::seq::mutation::{MutationEvent, Sub};

/// Results of an `ancestral` run.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct AncestralResults {
  /// Number of nucleotide mutations on all branches.
  pub mutations: usize,
  /// Branches with at least one mutation, most mutations first.
  pub branches: Vec<BranchMutations>,
  /// Sequence positions with substitutions between determined states on more than one branch, most
  /// branches first.
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

pub fn ancestral_results(tree: Option<&ResultTree>, alphabet: &Alphabet) -> Result<AncestralResults, Report> {
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
  let substitutions = tree
    .nodes
    .iter()
    .map(|node| {
      node
        .mutations
        .iter()
        .map(|mutation| substitution_event(mutation, alphabet))
        .flatten_ok()
        .collect::<Result<Vec<MutationEvent>, Report>>()
    })
    .collect::<Result<Vec<_>, Report>>()?;
  let recurrent_sites = sites_by_branch_count(&substitutions, alphabet, MutationClass::Substitution)
    .into_iter()
    .filter(|site| site.branches > 1)
    .map(|site| RecurrentSite {
      position: site.position + 1,
      branches: site.branches,
    })
    .collect();
  Ok(AncestralResults {
    mutations: branches.iter().map(|branch| branch.mutations.len()).sum(),
    branches,
    recurrent_sites,
  })
}

fn substitution_event(mutation: &str, alphabet: &Alphabet) -> Result<Option<MutationEvent>, Report> {
  if mutation.contains(&alphabet.gap().to_string()) {
    return Ok(None);
  }
  Ok(Some(MutationEvent::Substitution(Sub::from_str(mutation)?)))
}

use crate::alphabet::alphabet::Alphabet;
use crate::homoplasy::classify::{MutationClass, classify_mutation};
use crate::seq::mutation::MutationEvent;
use itertools::Itertools;
use std::borrow::Borrow;
use std::collections::{BTreeMap, BTreeSet};

pub fn sites_by_branch_count<Branch, Event>(
  branches: impl IntoIterator<Item = Branch>,
  alphabet: &Alphabet,
  class: MutationClass,
) -> Vec<SiteBranches>
where
  Branch: IntoIterator<Item = Event>,
  Event: Borrow<MutationEvent>,
{
  let mut branches_per_site: BTreeMap<usize, usize> = BTreeMap::new();
  for branch in branches {
    let positions: BTreeSet<usize> = branch
      .into_iter()
      .filter_map(|event| match event.borrow() {
        MutationEvent::Substitution(sub) if classify_mutation(event.borrow(), alphabet) == class => Some(sub.pos()),
        MutationEvent::Substitution(_) | MutationEvent::Insertion(_) | MutationEvent::Deletion(_) => None,
      })
      .collect();
    for position in positions {
      *branches_per_site.entry(position).or_default() += 1;
    }
  }
  branches_per_site
    .into_iter()
    .map(|(position, branches)| SiteBranches { position, branches })
    .sorted_by(|a, b| b.branches.cmp(&a.branches).then_with(|| a.position.cmp(&b.position)))
    .collect()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SiteBranches {
  pub position: usize,
  pub branches: usize,
}

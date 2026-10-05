use std::cmp::Ordering;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Recurrence<K> {
  pub mutation: K,
  pub branches: Vec<GraphNodeKey>,
}

impl<K> Recurrence<K> {
  pub fn multiplicity(&self) -> usize {
    self.branches.len()
  }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RecurrenceTable<K> {
  pub count: usize,
  pub histogram: BTreeMap<usize, usize>,
  pub ranked: Vec<Recurrence<K>>,
}

impl<K: Ord> RecurrenceTable<K> {
  pub(crate) fn new(
    events: impl IntoIterator<Item = (K, GraphNodeKey)>,
    tie_order: impl Fn(&K, &K) -> Ordering,
  ) -> Self {
    let mut branches_by_mutation: BTreeMap<K, Vec<GraphNodeKey>> = BTreeMap::new();
    for (mutation, branch) in events {
      branches_by_mutation.entry(mutation).or_default().push(branch);
    }
    let mut ranked: Vec<Recurrence<K>> = branches_by_mutation
      .into_iter()
      .map(|(mutation, branches)| Recurrence { mutation, branches })
      .collect();
    ranked.sort_by(|a, b| {
      b.multiplicity()
        .cmp(&a.multiplicity())
        .then_with(|| tie_order(&a.mutation, &b.mutation))
    });
    let mut histogram = BTreeMap::new();
    for recurrence in &ranked {
      *histogram.entry(recurrence.multiplicity()).or_default() += 1;
    }
    Self {
      count: ranked.iter().map(Recurrence::multiplicity).sum(),
      histogram,
      ranked,
    }
  }

  pub fn multiplicities(&self) -> BTreeMap<&K, usize> {
    self
      .ranked
      .iter()
      .map(|recurrence| (&recurrence.mutation, recurrence.multiplicity()))
      .collect()
  }
}

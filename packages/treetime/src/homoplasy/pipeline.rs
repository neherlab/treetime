use crate::alphabet::alphabet::Alphabet;
use crate::error::OperationError;
use crate::homoplasy::classify::{MutationClass, classify_mutation};
use crate::homoplasy::recurrence::RecurrenceTable;
use crate::homoplasy::site_hits::{SiteHistogram, site_histogram};
use crate::seq::indel::InDelKind;
use crate::seq::mutation::{AlignedMutation, Mutation, MutationEvent, Sub};
use std::cmp::Ordering;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

pub fn run(params: &HomoplasyParams, input: &HomoplasyInput<'_>) -> Result<HomoplasyOutput, OperationError> {
  let branches = branches(input);
  let genome_length = input.sequence_length + params.constant_sites;
  Ok(HomoplasyOutput {
    substitutions: substitution_stats(input, &branches, genome_length)?,
    ambiguous: ambiguous_stats(input, &branches),
    indels: indel_stats(input, &branches),
  })
}

#[derive(Clone, Copy, Debug)]
pub struct HomoplasyParams {
  pub constant_sites: usize,
}

pub struct HomoplasyInput<'a> {
  pub graph: &'a Graph,
  pub bridged_mutations: &'a BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  pub raw_mutations: &'a BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  pub branch_lengths: &'a BTreeMap<GraphEdgeKey, f64>,
  pub alphabet: &'a Alphabet,
  pub sequence_length: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct HomoplasyOutput {
  pub substitutions: SubstitutionStats,
  pub ambiguous: AmbiguousStats,
  pub indels: IndelStats,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SubstitutionStats {
  pub total_branch_length: f64,
  pub terminal_branch_length: f64,
  pub all: RecurrenceTable<Sub>,
  pub terminal: RecurrenceTable<Sub>,
  pub sites: SiteHistogram,
  pub leaves: BTreeMap<GraphNodeKey, Vec<Sub>>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AmbiguousStats {
  pub all: RecurrenceTable<Sub>,
  pub sites: Vec<SiteBranches>,
  pub leaves: BTreeMap<GraphNodeKey, usize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SiteBranches {
  pub position: usize,
  pub branches: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IndelStats {
  pub all: RecurrenceTable<IndelKey>,
  pub terminal: RecurrenceTable<IndelKey>,
  pub leaves: BTreeMap<GraphNodeKey, usize>,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct IndelKey {
  pub kind: InDelKind,
  pub range: (usize, usize),
  pub sequence: Seq,
}

impl IndelKey {
  fn from_event(event: &MutationEvent) -> Option<Self> {
    let (kind, AlignedMutation { range, sequence }) = match event {
      MutationEvent::Insertion(segment) => (InDelKind::Insertion, segment),
      MutationEvent::Deletion(segment) => (InDelKind::Deletion, segment),
      MutationEvent::Substitution(_) => return None,
    };
    Some(Self {
      kind,
      range: *range,
      sequence: sequence.clone(),
    })
  }
}

struct Branch {
  edge: GraphEdgeKey,
  node: GraphNodeKey,
  terminal: bool,
}

fn branches(input: &HomoplasyInput<'_>) -> Vec<Branch> {
  input
    .graph
    .get_edges()
    .map(|edge| Branch {
      edge: edge.key(),
      node: edge.target(),
      terminal: input.graph.is_leaf(edge.target()),
    })
    .collect()
}

fn events_of_class<'a>(
  input: &'a HomoplasyInput<'_>,
  mutations: &'a BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  branch: &'a Branch,
  class: MutationClass,
) -> impl Iterator<Item = &'a MutationEvent> + 'a {
  mutations[&branch.edge]
    .iter()
    .map(|mutation| &mutation.event)
    .filter(move |event| classify_mutation(event, input.alphabet) == class)
}

fn substitutions<'a>(
  input: &'a HomoplasyInput<'_>,
  mutations: &'a BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  branch: &'a Branch,
  class: MutationClass,
) -> impl Iterator<Item = &'a Sub> + 'a {
  events_of_class(input, mutations, branch, class).filter_map(|event| match event {
    MutationEvent::Substitution(sub) => Some(sub),
    MutationEvent::Insertion(_) | MutationEvent::Deletion(_) => None,
  })
}

fn substitution_order(a: &Sub, b: &Sub) -> Ordering {
  (a.pos(), a.reff(), a.qry()).cmp(&(b.pos(), b.reff(), b.qry()))
}

fn indel_order(a: &IndelKey, b: &IndelKey) -> Ordering {
  (a.range, a.kind, &a.sequence).cmp(&(b.range, b.kind, &b.sequence))
}

fn substitution_stats(
  input: &HomoplasyInput<'_>,
  branches: &[Branch],
  genome_length: usize,
) -> Result<SubstitutionStats, OperationError> {
  let class = MutationClass::Substitution;
  let mutations = input.bridged_mutations;
  let events = || {
    branches.iter().flat_map(move |branch| {
      substitutions(input, mutations, branch, class).map(move |sub| (sub.clone(), branch))
    })
  };
  let all = RecurrenceTable::new(events().map(|(sub, branch)| (sub, branch.node)), substitution_order);
  let terminal = RecurrenceTable::new(
    events()
      .filter(|(_, branch)| branch.terminal)
      .map(|(sub, branch)| (sub, branch.node)),
    substitution_order,
  );

  let mut hits_per_site: BTreeMap<usize, usize> = BTreeMap::new();
  for (sub, _) in events() {
    *hits_per_site.entry(sub.pos()).or_default() += 1;
  }
  let mut leaves: BTreeMap<GraphNodeKey, Vec<Sub>> = BTreeMap::new();
  for (sub, branch) in events() {
    if branch.terminal && hits_per_site[&sub.pos()] > 1 {
      leaves.entry(branch.node).or_default().push(sub);
    }
  }

  let total_branch_length = branches.iter().map(|branch| input.branch_lengths[&branch.edge]).sum();
  let terminal_branch_length = branches
    .iter()
    .filter(|branch| branch.terminal)
    .map(|branch| odd_substitution_probability(input.branch_lengths[&branch.edge]))
    .sum();

  Ok(SubstitutionStats {
    total_branch_length,
    terminal_branch_length,
    all,
    terminal,
    sites: site_histogram(genome_length, &hits_per_site)?,
    leaves,
  })
}

fn odd_substitution_probability(branch_length: f64) -> f64 {
  -0.5 * (-2.0 * branch_length).exp_m1()
}

fn ambiguous_stats(input: &HomoplasyInput<'_>, branches: &[Branch]) -> AmbiguousStats {
  let class = MutationClass::Ambiguous;
  let mutations = input.raw_mutations;
  let events = || {
    branches.iter().flat_map(move |branch| {
      substitutions(input, mutations, branch, class).map(move |sub| (sub.clone(), branch))
    })
  };
  let all = RecurrenceTable::new(events().map(|(sub, branch)| (sub, branch.node)), substitution_order);

  let mut branches_per_site: BTreeMap<usize, usize> = BTreeMap::new();
  let mut leaves: BTreeMap<GraphNodeKey, usize> = BTreeMap::new();
  for (sub, branch) in events() {
    *branches_per_site.entry(sub.pos()).or_default() += 1;
    if branch.terminal {
      *leaves.entry(branch.node).or_default() += 1;
    }
  }
  let mut sites: Vec<SiteBranches> = branches_per_site
    .into_iter()
    .map(|(position, branches)| SiteBranches { position, branches })
    .collect();
  sites.sort_by(|a, b| b.branches.cmp(&a.branches).then_with(|| a.position.cmp(&b.position)));

  AmbiguousStats { all, sites, leaves }
}

fn indel_stats(input: &HomoplasyInput<'_>, branches: &[Branch]) -> IndelStats {
  let mutations = input.raw_mutations;
  let events = || {
    branches.iter().flat_map(move |branch| {
      events_of_class(input, mutations, branch, MutationClass::Indel)
        .filter_map(IndelKey::from_event)
        .map(move |key| (key, branch))
    })
  };
  let all = RecurrenceTable::new(events().map(|(key, branch)| (key, branch.node)), indel_order);
  let terminal = RecurrenceTable::new(
    events()
      .filter(|(_, branch)| branch.terminal)
      .map(|(key, branch)| (key, branch.node)),
    indel_order,
  );

  let multiplicities = all.multiplicities();
  let mut leaves: BTreeMap<GraphNodeKey, usize> = BTreeMap::new();
  for (key, branch) in events() {
    if branch.terminal && multiplicities[&key] > 1 {
      *leaves.entry(branch.node).or_default() += 1;
    }
  }

  IndelStats { all, terminal, leaves }
}

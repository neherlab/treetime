use crate::alphabet::alphabet::Alphabet;
use crate::make_error;
use crate::partition::fitch::passes::create_fitch_partition;
use crate::seq::alignment::NodeSeqInput;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::clades::{TipIds, graph_clade_keys};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node_names::NodeNames;
use treetime_primitives::Seq;

pub fn fitch_edge_counts(
  graph: &Graph,
  names: &NodeNames,
  seqs: &BTreeMap<&str, &Seq>,
  alphabet: &Alphabet,
) -> Result<BTreeMap<GraphEdgeKey, usize>, Report> {
  Ok(
    fitch_edge_subs(graph, names, seqs, alphabet)?
      .into_iter()
      .map(|(edge, subs)| (edge, subs.len()))
      .collect(),
  )
}

pub fn fitch_edge_subs(
  graph: &Graph,
  names: &NodeNames,
  seqs: &BTreeMap<&str, &Seq>,
  alphabet: &Alphabet,
) -> Result<BTreeMap<GraphEdgeKey, Vec<Sub>>, Report> {
  let inputs = graph
    .get_leaves()
    .map(|leaf| {
      let name = &names[leaf.key()];
      let Some(&seq) = seqs.get(name) else {
        return make_error!("Tree tip '{name}' has no sequence");
      };
      Ok((
        leaf.key(),
        NodeSeqInput {
          name: name.to_owned(),
          seq: Some(seq.clone()),
        },
      ))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  let fitch = create_fitch_partition(graph, alphabet.clone(), inputs)?;
  Ok(
    fitch
      .edges
      .iter()
      .map(|(&edge, obs)| (edge, obs.fitch_subs().to_vec()))
      .collect(),
  )
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SplitComparison {
  pub shared: usize,
  pub only_first: usize,
  pub only_second: usize,
}

pub struct TreeWithCounts<'a> {
  pub graph: &'a Graph,
  pub names: &'a NodeNames,
  pub edge_counts: &'a BTreeMap<GraphEdgeKey, usize>,
}

pub fn compare_splits(first: &TreeWithCounts<'_>, second: &TreeWithCounts<'_>) -> Result<SplitComparison, Report> {
  let tips = |tree: &TreeWithCounts<'_>| {
    tree
      .graph
      .get_leaves()
      .map(|leaf| tree.names[leaf.key()].to_owned())
      .sorted()
      .collect_vec()
  };
  let first_tips = tips(first);
  if first_tips != tips(second) {
    return make_error!("The two trees have different tip sets");
  }
  let mut ids = TipIds::default();
  for name in &first_tips {
    ids.id(name);
  }
  let first_splits = informative_splits(first, &mut ids, first_tips.len())?;
  let second_splits = informative_splits(second, &mut ids, first_tips.len())?;
  Ok(SplitComparison {
    shared: first_splits.intersection(&second_splits).count(),
    only_first: first_splits.difference(&second_splits).count(),
    only_second: second_splits.difference(&first_splits).count(),
  })
}

fn informative_splits<'a>(
  tree: &TreeWithCounts<'a>,
  ids: &mut TipIds<'a>,
  n_tips: usize,
) -> Result<BTreeSet<Vec<usize>>, Report> {
  let clades = graph_clade_keys(tree.graph, tree.names, ids)?;
  let mut splits = BTreeSet::new();
  for (node, tips) in clades {
    let Some((_, edge)) = tree.graph.node_parent(node)? else {
      continue;
    };
    if tips.len() < 2 || tips.len() + 1 >= n_tips || tree.edge_counts.get(&edge).is_none_or(|&count| count == 0) {
      continue;
    }
    let split = if tips.first() == Some(&0) {
      (0..n_tips).filter(|id| tips.binary_search(id).is_err()).collect_vec()
    } else {
      tips
    };
    splits.insert(split);
  }
  Ok(splits)
}

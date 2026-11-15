use crate::partition::assembly::AlignmentGroup;
use crate::partition::marginal::reconstruction::ReconstructionRef;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::set::PartitionSet;
use eyre::Report;
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use smallvec::SmallVec;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::EdgeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;

pub(crate) struct TimelessNode {
  pub(crate) held: GraphEdgeKey,
  pub(crate) others: Vec<GraphEdgeKey>,
}

pub(crate) trait EdgeContributions: Sync {
  fn with_edge<R>(
    &self,
    edge_key: GraphEdgeKey,
    f: impl FnOnce(&[OptimizationContribution]) -> Result<R, Report>,
  ) -> Result<R, Report>;

  fn timeless_nodes(&self, _graph: &Graph) -> Result<Vec<TimelessNode>, Report> {
    Ok(vec![])
  }

  fn frozen_edges(&self) -> Option<&BTreeSet<GraphEdgeKey>> {
    None
  }
}

impl EdgeContributions for PartitionSet {
  fn with_edge<R>(
    &self,
    edge_key: GraphEdgeKey,
    f: impl FnOnce(&[OptimizationContribution]) -> Result<R, Report>,
  ) -> Result<R, Report> {
    let contributions: SmallVec<[OptimizationContribution; 1]> = self
      .iter()
      .filter(|partition| partition.carries(edge_key))
      .map(|partition| partition.create_edge_contribution(edge_key))
      .try_collect()?;
    f(&contributions)
  }

  fn timeless_nodes(&self, graph: &Graph) -> Result<Vec<TimelessNode>, Report> {
    let schedule = GraphSchedule::new(graph)?;
    let passes_through = |key| {
      self.iter().all(|partition| {
        let inheritance = partition.inheritance();
        !inheritance.holds(&schedule, key)
          || (inheritance.parent(&schedule, key).is_some() && inheritance.children(&schedule, key).count() == 1)
      })
    };
    let edges = |pairs: &[(GraphNodeKey, GraphEdgeKey)]| pairs.iter().map(|&(_, edge)| edge).collect_vec();
    Ok(
      schedule
        .preorder()
        .iter()
        .filter(|&&key| !schedule.children(key).is_empty() && passes_through(key))
        .filter_map(|&key| match (schedule.parents(key), schedule.children(key)) {
          (parents, &[(_, child)]) if parents.len() > 1 => Some(TimelessNode {
            held: child,
            others: edges(parents),
          }),
          (&[(_, parent)], children) if children.len() > 1 => Some(TimelessNode {
            held: parent,
            others: edges(children),
          }),
          _ => None,
        })
        .collect(),
    )
  }
}

pub(crate) struct EdgeGroup<'a> {
  pub(crate) partitions: &'a PartitionSet,
  pub(crate) frozen: BTreeSet<GraphEdgeKey>,
}

impl EdgeContributions for EdgeGroup<'_> {
  fn with_edge<R>(
    &self,
    edge_key: GraphEdgeKey,
    f: impl FnOnce(&[OptimizationContribution]) -> Result<R, Report>,
  ) -> Result<R, Report> {
    self.partitions.with_edge(edge_key, f)
  }

  fn timeless_nodes(&self, graph: &Graph) -> Result<Vec<TimelessNode>, Report> {
    self.partitions.timeless_nodes(graph)
  }

  fn frozen_edges(&self) -> Option<&BTreeSet<GraphEdgeKey>> {
    Some(&self.frozen)
  }
}

pub(crate) fn series_groups(graph: &Graph, partitions: &PartitionSet) -> Result<Vec<BTreeSet<GraphEdgeKey>>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  let mut conflicts: BTreeMap<GraphEdgeKey, BTreeSet<GraphEdgeKey>> = BTreeMap::new();
  for &key in schedule.preorder() {
    if schedule.parents(key).len() + schedule.children(key).len() == 2 {
      continue;
    }
    for partition in partitions.iter() {
      let inheritance = partition.inheritance();
      let children = inheritance
        .children(&schedule, key)
        .map(|(_, _, edge)| edge)
        .collect_vec();
      let pair = match (inheritance.parent(&schedule, key), children.as_slice()) {
        (Some((_, _, parent_edge)), &[child_edge]) => Some((parent_edge, child_edge)),
        (None, &[first, second]) => Some((first, second)),
        _ => None,
      };
      if let Some((a, b)) = pair {
        conflicts.entry(a).or_default().insert(b);
        conflicts.entry(b).or_default().insert(a);
      }
    }
  }
  let all: BTreeSet<GraphEdgeKey> = graph.get_edges().map(|edge| edge.key()).collect();
  if conflicts.is_empty() {
    return Ok(vec![all]);
  }
  let mut color: BTreeMap<GraphEdgeKey, usize> = BTreeMap::new();
  for &edge in &all {
    let used: BTreeSet<usize> = conflicts
      .get(&edge)
      .into_iter()
      .flatten()
      .filter_map(|other| color.get(other).copied())
      .collect();
    let free = (0..=used.len()).find(|c| !used.contains(c)).unwrap_or(0);
    color.insert(edge, free);
  }
  let count = color.values().max().map_or(1, |max| max + 1);
  Ok(
    (0..count)
      .map(|c| color.iter().filter(|&(_, &k)| k == c).map(|(&edge, _)| edge).collect())
      .collect(),
  )
}

#[cfg(test)]
impl EdgeContributions for EdgeMap<Vec<OptimizationContribution>> {
  fn with_edge<R>(
    &self,
    edge_key: GraphEdgeKey,
    f: impl FnOnce(&[OptimizationContribution]) -> Result<R, Report>,
  ) -> Result<R, Report> {
    f(&self[&edge_key])
  }
}

#[cfg(test)]
pub(crate) fn gather_edge_contributions(
  graph: &Graph,
  partitions: &PartitionSet,
) -> Result<EdgeMap<Vec<OptimizationContribution>>, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(
      |edge_ref| -> Result<(GraphEdgeKey, Vec<OptimizationContribution>), Report> {
        let edge_key = edge_ref.key();
        let contributions = partitions
          .iter()
          .filter(|partition| partition.carries(edge_key))
          .map(|partition| partition.create_edge_contribution(edge_key))
          .try_collect()?;
        Ok((edge_key, contributions))
      },
    )
    .collect()
}

pub(crate) struct InitialGuessEvidence {
  pub total_length: usize,
  pub sub_counts: EdgeMap<usize>,
  pub effective_lengths: EdgeMap<usize>,
}

impl InitialGuessEvidence {
  pub(crate) fn gather(graph: &Graph, partitions: &PartitionSet) -> Result<Self, Report> {
    Ok(Self {
      total_length: partitions.total_length(),
      sub_counts: gather_edge_sums(graph, partitions, |partition, edge_key| {
        Ok(partition.edge_subs(graph, edge_key)?.len())
      })?,
      effective_lengths: gather_edge_sums(graph, partitions, |partition, edge_key| {
        partition.edge_effective_length(graph, edge_key)
      })?,
    })
  }
}

pub(crate) fn gather_edge_indel_counts(graph: &Graph, partitions: &PartitionSet) -> Result<EdgeMap<usize>, Report> {
  Ok(sum_alignment_counts(
    graph,
    &gather_alignment_edge_indel_counts(graph, partitions)?,
  ))
}

pub(crate) fn gather_alignment_edge_indel_counts(
  graph: &Graph,
  partitions: &PartitionSet,
) -> Result<BTreeMap<usize, EdgeMap<usize>>, Report> {
  let reconstructions = partitions.iter().collect_vec();
  partitions
    .alignment_groups()?
    .iter()
    .map(|group| {
      let members = group.members(&reconstructions);
      let counts = graph
        .get_edges()
        .collect::<Vec<_>>()
        .into_par_iter()
        .map(|edge_ref| {
          let edge_key = edge_ref.key();
          Ok((edge_key, alignment_edge_indel_count(group, &members, edge_key)?))
        })
        .collect::<Result<_, Report>>()?;
      Ok((group.alignment(), counts))
    })
    .collect()
}

pub(crate) fn sum_alignment_counts(graph: &Graph, counts: &BTreeMap<usize, EdgeMap<usize>>) -> EdgeMap<usize> {
  graph
    .get_edges()
    .map(|edge_ref| {
      let edge_key = edge_ref.key();
      (edge_key, counts.values().map(|counts| counts[&edge_key]).sum())
    })
    .collect()
}

fn gather_edge_sums(
  graph: &Graph,
  partitions: &PartitionSet,
  value: impl Fn(ReconstructionRef<'_>, GraphEdgeKey) -> Result<usize, Report> + Sync,
) -> Result<EdgeMap<usize>, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| -> Result<(GraphEdgeKey, usize), Report> {
      let edge_key = edge_ref.key();
      let sum = partitions
        .iter()
        .filter(|partition| partition.carries(edge_key))
        .map(|partition| value(partition, edge_key))
        .sum::<Result<usize, Report>>()?;
      Ok((edge_key, sum))
    })
    .collect()
}

fn alignment_edge_indel_count(
  group: &AlignmentGroup,
  members: &[&ReconstructionRef<'_>],
  edge_key: GraphEdgeKey,
) -> Result<usize, Report> {
  if !members.iter().any(|member| member.carries(edge_key)) {
    return Ok(0);
  }
  if let [member] = members
    && group.layout.is_whole()
  {
    return Ok(member.edge_indel_count(edge_key));
  }
  let indels = members.iter().map(|member| member.edge_indels(edge_key)).collect_vec();
  Ok(group.layout.assemble_indels(&indels)?.len())
}

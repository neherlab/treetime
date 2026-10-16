use crate::branch_lengths::one_mutation;
use crate::optimize::branch_length::{is_valid_branch_length_value, validate_branch_length_value};
use crate::optimize::gather::EdgeContributions;
use crate::optimize::indel::{EdgeIndels, IndelTerm};
use crate::optimize::method_brent::brent_sqrt_inner;
use crate::optimize::method_newton::newton_sqrt_inner;
use crate::optimize::params::{BranchOptMethod, ExistingBranchLengths};
use crate::optimize::zero_boundary::{is_zero_branch_optimal, min_branch_length_for_indels, reconcile_zero_boundary};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::{make_error, make_internal_report, make_report};
use eyre::{Report, WrapErr};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::{BTreeMap, BTreeSet};
use std::mem::take;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_utils::cast::ToF64Exact;

pub(crate) fn run_optimize_mixed(
  graph: &Graph,
  total_length: usize,
  contributions: &impl EdgeContributions,
  indels: &EdgeIndels,
  method: BranchOptMethod,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot optimize branch lengths");
  }

  let one_mutation = one_mutation(total_length);

  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.key();
    let branch_length = branch_lengths[&edge_key]
      .ok_or_else(|| make_report!("Cannot optimize edge {edge_key} with a missing branch length"))?;
    validate_branch_length_value(branch_length).wrap_err_with(|| format!("Cannot optimize edge {edge_key}"))?;
  }

  let chains = SeriesChain::capture_all(graph, branch_lengths);
  let held: BTreeSet<GraphEdgeKey> = chains.iter().flat_map(SeriesChain::held).collect();
  let chain_indels: BTreeMap<GraphEdgeKey, IndelTerm> = chains
    .iter()
    .map(|chain| (chain.optimized, chain.indel_term(indels, branch_lengths)))
    .collect();

  let branch_lengths_in = &*branch_lengths;
  let optimized: Vec<(GraphEdgeKey, f64)> = graph
    .get_edges()
    .filter(|edge_ref| !held.contains(&edge_ref.key()))
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| -> Result<(GraphEdgeKey, f64), Report> {
      let edge_key = edge_ref.key();
      let branch_length = branch_lengths_in[&edge_key]
        .ok_or_else(|| make_internal_report!("Validated edge {edge_key} lost its branch length"))?;
      let edge_indels = chain_indels
        .get(&edge_key)
        .copied()
        .unwrap_or_else(|| indels.term(edge_key));
      let edge = EdgeOptimization {
        branch_length,
        indels: edge_indels,
        one_mutation,
        method,
      };
      let new_branch_length = contributions.with_edge(edge_key, |contributions| optimize_edge(contributions, &edge))?;
      Ok((edge_key, new_branch_length))
    })
    .collect::<Result<Vec<_>, Report>>()?;

  for (edge_key, new_branch_length) in optimized {
    branch_lengths.insert(edge_key, Some(new_branch_length));
  }

  for chain in chains {
    let indels = chain_indels[&chain.optimized];
    chain.restore(indels, branch_lengths);
  }

  Ok(())
}

struct SeriesChain {
  edges: Vec<GraphEdgeKey>,
  ratios: Vec<f64>,
  optimized: GraphEdgeKey,
}

impl SeriesChain {
  fn capture_all(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Vec<Self> {
    let mut chain_of: BTreeMap<GraphEdgeKey, usize> = BTreeMap::new();
    let mut groups: Vec<Vec<GraphEdgeKey>> = vec![];
    for node in graph.get_nodes() {
      let incident = [node.inbound(), node.outbound()].concat();
      if node.is_leaf() || incident.len() != 2 {
        continue;
      }
      match (chain_of.get(&incident[0]).copied(), chain_of.get(&incident[1]).copied()) {
        (None, None) => {
          incident.iter().for_each(|&e| _ = chain_of.insert(e, groups.len()));
          groups.push(incident);
        },
        (Some(g), None) | (None, Some(g)) => {
          for e in incident {
            if chain_of.insert(e, g).is_none() {
              groups[g].push(e);
            }
          }
        },
        (Some(g0), Some(g1)) => {
          let moved = take(&mut groups[g1]);
          moved.iter().for_each(|&e| _ = chain_of.insert(e, g0));
          groups[g0].extend(moved);
        },
      }
    }

    groups
      .into_iter()
      .filter(|edges| !edges.is_empty())
      .map(|edges| {
        let lengths: Vec<f64> = edges.iter().map(|e| branch_lengths[e].unwrap_or(0.0)).collect();
        let total: f64 = lengths.iter().sum();
        let ratios = if total > 0.0 {
          lengths.iter().map(|bl| bl / total).collect()
        } else {
          vec![1.0 / edges.len().to_f64_exact(); edges.len()]
        };
        let longest = (0..edges.len()).fold(0, |best, i| if lengths[i] > lengths[best] { i } else { best });
        Self {
          optimized: edges[longest],
          edges,
          ratios,
        }
      })
      .collect()
  }

  fn held(&self) -> impl Iterator<Item = GraphEdgeKey> + '_ {
    self.edges.iter().copied().filter(|&e| e != self.optimized)
  }

  fn indel_term(&self, indels: &EdgeIndels, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> IndelTerm {
    let held_length = self.held().map(|e| branch_lengths[&e].unwrap_or(0.0)).sum();
    indels.chain_term(&self.edges, held_length)
  }

  fn restore(self, indels: IndelTerm, branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>) {
    let total = if branch_lengths[&self.optimized] == Some(0.0) && indels.count == 0 {
      0.0
    } else {
      self.edges.iter().map(|e| branch_lengths[e].unwrap_or(0.0)).sum()
    };
    for (edge, ratio) in self.edges.into_iter().zip(self.ratios) {
      branch_lengths.insert(edge, Some(total * ratio));
    }
  }
}

pub(crate) fn initial_guess_mixed(
  graph: &Graph,
  total_length: usize,
  indels: &EdgeIndels,
  sub_counts: &BTreeMap<GraphEdgeKey, usize>,
  effective_lengths: &BTreeMap<GraphEdgeKey, usize>,
  existing: ExistingBranchLengths,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot compute initial guess");
  }

  let one_mutation = one_mutation(total_length);

  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.key();

    let edge_indels = indels.term(edge_key);

    if existing == ExistingBranchLengths::Keep {
      if let Some(bl) = branch_lengths.get(&edge_key).copied().flatten() {
        if is_valid_branch_length_value(bl) && (bl > 0.0 || !edge_indels.requires_positive_length()) {
          continue;
        }
      }
    }

    let sub_count: usize = sub_counts[&edge_key];

    let effective_length: usize = effective_lengths[&edge_key];

    let branch_length = if effective_length > 0 {
      let sub_estimate = sub_count.to_f64_exact() / effective_length.to_f64_exact();
      if sub_estimate == 0.0 && edge_indels.count > 0 {
        edge_indels.length_estimate().unwrap_or(one_mutation)
      } else {
        sub_estimate
      }
    } else if edge_indels.count > 0 {
      one_mutation
    } else {
      one_mutation * 0.1
    };
    branch_lengths.insert(edge_key, Some(branch_length));
  }

  Ok(())
}

fn optimize_edge(contributions: &[OptimizationContribution], edge: &EdgeOptimization) -> Result<f64, Report> {
  let EdgeOptimization {
    mut branch_length,
    indels,
    one_mutation,
    method,
  } = *edge;

  if branch_length == 0.0 && indels.requires_positive_length() {
    branch_length = indels
      .length_estimate()
      .map_or(one_mutation, |length| length.max(one_mutation));
  }

  if branch_length == 0.0 && !contributions.iter().all(|c| c.all_sites_valid_at_zero()) {
    branch_length = one_mutation;
  }

  if indels.count == 0 && is_zero_branch_optimal(contributions) {
    return Ok(0.0);
  }

  let min_branch_length = min_branch_length_for_indels(indels, one_mutation);

  let new_branch_length = match method {
    BranchOptMethod::Brent => brent_sqrt_inner(branch_length, contributions, indels, min_branch_length, one_mutation),
    BranchOptMethod::Newton => newton_sqrt_inner(branch_length, contributions, indels, min_branch_length, one_mutation),
  }?;

  let new_branch_length =
    reconcile_zero_boundary(new_branch_length, branch_length, contributions, indels, one_mutation)?;

  Ok(new_branch_length)
}

#[derive(Clone, Copy)]
struct EdgeOptimization {
  branch_length: f64,
  indels: IndelTerm,
  one_mutation: f64,
  method: BranchOptMethod,
}

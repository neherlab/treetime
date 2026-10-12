use crate::branch_lengths::one_mutation;
use crate::optimize::branch_length::{is_valid_branch_length_value, validate_branch_length_value};
use crate::optimize::gather::EdgeContributions;
use crate::optimize::indel::estimate_indel_rate;
use crate::optimize::method_brent::brent_sqrt_inner;
use crate::optimize::method_newton::newton_sqrt_inner;
use crate::optimize::params::{BranchOptMethod, ExistingBranchLengths};
use crate::optimize::zero_boundary::{is_zero_branch_optimal, min_branch_length_for_indels, reconcile_zero_boundary};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::{make_error, make_internal_report, make_report};
use eyre::{Report, WrapErr};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

pub(crate) fn run_optimize_mixed(
  graph: &Graph,
  total_length: usize,
  contributions: &impl EdgeContributions,
  indel_counts: &BTreeMap<GraphEdgeKey, usize>,
  method: BranchOptMethod,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot optimize branch lengths");
  }

  let indel_rate = estimate_indel_rate(graph, indel_counts, branch_lengths);
  run_optimize_mixed_inner(
    graph,
    total_length,
    contributions,
    indel_counts,
    method,
    indel_rate,
    false,
    branch_lengths,
  )
}

#[expect(
  clippy::too_many_arguments,
  reason = "each argument is an independent input of this step; a parameter struct would be built only for this call"
)]
pub(crate) fn run_optimize_mixed_inner(
  graph: &Graph,
  total_length: usize,
  contributions: &impl EdgeContributions,
  indel_counts: &BTreeMap<GraphEdgeKey, usize>,
  method: BranchOptMethod,
  indel_rate: f64,
  no_indels: bool,
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
  // The likelihood depends only on the sum of edges joined by a node without data of its own and with two edges (a
  // bifurcating root or a single-child node). Optimizing one edge of such a chain with the others held at their current
  // lengths therefore optimizes the sum; the held edges are skipped and the sum is redistributed by `restore`.
  let held: BTreeSet<GraphEdgeKey> = chains.iter().flat_map(SeriesChain::held).collect();
  let chain_indel_counts: BTreeMap<GraphEdgeKey, usize> = chains
    .iter()
    .map(|chain| (chain.optimized, chain.edges.iter().map(|e| indel_counts[e]).sum()))
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
      let indel_count = chain_indel_counts
        .get(&edge_key)
        .copied()
        .unwrap_or_else(|| indel_counts[&edge_key]);
      let edge = EdgeOptimization {
        branch_length,
        indel_count: if no_indels { 0 } else { indel_count },
        indel_rate,
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
    chain.restore(branch_lengths);
  }

  Ok(())
}

/// Edges in series: joined by nodes without data that have exactly two edges (bifurcating root, single-child nodes).
struct SeriesChain {
  edges: Vec<GraphEdgeKey>,
  ratios: Vec<f64>,
  /// Edge whose optimization yields the optimal sum (the longest one, so the zero clamp rarely binds).
  optimized: GraphEdgeKey,
}

impl SeriesChain {
  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
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
          let moved = std::mem::take(&mut groups[g1]);
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
          vec![1.0 / edges.len() as f64; edges.len()]
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

  fn restore(self, branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>) {
    // A zero optimum of the optimized edge only bounds the sum by the held lengths; collapse it and let the next
    // iteration, with all chain edges at zero, find the sum exactly.
    let total = if branch_lengths[&self.optimized] == Some(0.0) {
      0.0
    } else {
      self.edges.iter().map(|e| branch_lengths[e].unwrap_or(0.0)).sum()
    };
    for (edge, ratio) in self.edges.into_iter().zip(self.ratios) {
      branch_lengths.insert(edge, Some(total * ratio));
    }
  }
}

#[expect(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
#[expect(
  clippy::too_many_arguments,
  reason = "each argument is an independent input of this step; a parameter struct would be built only for this call"
)]
pub(crate) fn initial_guess_mixed(
  graph: &Graph,
  total_length: usize,
  indel_counts: &BTreeMap<GraphEdgeKey, usize>,
  sub_counts: &BTreeMap<GraphEdgeKey, usize>,
  effective_lengths: &BTreeMap<GraphEdgeKey, usize>,
  existing: ExistingBranchLengths,
  no_indels: bool,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot compute initial guess");
  }

  let one_mutation = one_mutation(total_length);
  let indel_rate = if no_indels {
    0.0
  } else {
    estimate_indel_rate(graph, indel_counts, branch_lengths)
  };

  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.key();

    let indel_count: usize = if no_indels { 0 } else { indel_counts[&edge_key] };

    if existing == ExistingBranchLengths::Keep {
      if let Some(bl) = branch_lengths.get(&edge_key).copied().flatten() {
        if is_valid_branch_length_value(bl) && (bl > 0.0 || indel_count == 0) {
          continue;
        }
      }
    }

    let sub_count: usize = sub_counts[&edge_key];

    let effective_length: usize = effective_lengths[&edge_key];

    let branch_length = if effective_length > 0 {
      let sub_estimate = sub_count as f64 / effective_length as f64;
      if sub_estimate == 0.0 && indel_count > 0 {
        if indel_rate > 0.0 {
          indel_count as f64 / indel_rate
        } else {
          one_mutation
        }
      } else {
        sub_estimate
      }
    } else if indel_count > 0 {
      one_mutation
    } else {
      one_mutation * 0.1
    };
    branch_lengths.insert(edge_key, Some(branch_length));
  }

  Ok(())
}

#[derive(Clone, Copy)]
struct EdgeOptimization {
  branch_length: f64,
  indel_count: usize,
  indel_rate: f64,
  one_mutation: f64,
  method: BranchOptMethod,
}

#[expect(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
fn optimize_edge(contributions: &[OptimizationContribution], edge: &EdgeOptimization) -> Result<f64, Report> {
  let EdgeOptimization {
    mut branch_length,
    indel_count,
    indel_rate,
    one_mutation,
    method,
  } = *edge;

  if branch_length == 0.0 && indel_count > 0 {
    branch_length = if indel_rate > 0.0 {
      (indel_count as f64 / indel_rate).max(one_mutation)
    } else {
      one_mutation
    };
  }

  if branch_length == 0.0 && !contributions.iter().all(|c| c.all_sites_valid_at_zero()) {
    branch_length = one_mutation;
  }

  if indel_count == 0 && is_zero_branch_optimal(contributions) {
    return Ok(0.0);
  }

  let min_branch_length = min_branch_length_for_indels(indel_count, one_mutation);

  let new_branch_length = match method {
    BranchOptMethod::BrentSqrt => brent_sqrt_inner(
      branch_length,
      contributions,
      indel_count,
      indel_rate,
      min_branch_length,
      one_mutation,
    ),
    BranchOptMethod::NewtonSqrt => newton_sqrt_inner(
      branch_length,
      contributions,
      indel_count,
      indel_rate,
      min_branch_length,
      one_mutation,
    ),
  }?;

  let new_branch_length = reconcile_zero_boundary(
    new_branch_length,
    branch_length,
    contributions,
    indel_count,
    indel_rate,
    one_mutation,
  )?;

  Ok(new_branch_length)
}

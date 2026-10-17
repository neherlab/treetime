use crate::branch_lengths::one_mutation;
use crate::optimize::branch_length::{is_valid_branch_length_value, validate_branch_length_value};
use crate::optimize::gather::{EdgeContributions, InitialGuessEvidence};
use crate::optimize::indel::{EdgeIndels, IndelTerm};
use crate::optimize::method_brent::solve_brent;
use crate::optimize::method_newton::solve_newton;
use crate::optimize::params::{BranchOptMethod, ExistingBranchLengths};
use crate::optimize::solver_outcome::{Convergence, EdgeSolution, UnconvergedEdges};
use crate::optimize::zero_boundary::{
  ZeroLength, is_zero_branch_optimal, min_branch_length_for_indels, reconcile_zero_boundary, start_branch_length,
};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::{make_error, make_report};
use eyre::{Report, WrapErr};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::{BTreeMap, BTreeSet};
use std::mem::take;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_utils::cast::ToF64Exact;

pub(crate) fn optimize_edge_lengths(
  graph: &Graph,
  total_length: usize,
  contributions: &impl EdgeContributions,
  indels: &EdgeIndels,
  method: BranchOptMethod,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<UnconvergedEdges, Report> {
  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot optimize branch lengths");
  }

  let one_mutation = one_mutation(total_length);
  let mut lengths = parse_optimizable_lengths(graph, branch_lengths)?;

  let chains = SeriesChain::capture_all(graph, &lengths);
  let held: BTreeSet<GraphEdgeKey> = chains.iter().flat_map(SeriesChain::held).collect();
  let chain_indels: BTreeMap<GraphEdgeKey, IndelTerm> = chains
    .iter()
    .map(|chain| (chain.optimized, chain.indel_term(indels, &lengths)))
    .collect();

  let lengths_in = &lengths;
  let optimized: Vec<(GraphEdgeKey, EdgeSolution)> = graph
    .get_edges()
    .filter(|edge_ref| !held.contains(&edge_ref.key()))
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| -> Result<(GraphEdgeKey, EdgeSolution), Report> {
      let edge_key = edge_ref.key();
      let edge_indels = chain_indels
        .get(&edge_key)
        .copied()
        .unwrap_or_else(|| indels.term(edge_key));
      let edge = EdgeOptimization {
        branch_length: lengths_in[&edge_key],
        indels: edge_indels,
        one_mutation,
        method,
      };
      let solution = contributions.with_edge(edge_key, |contributions| optimize_edge(contributions, &edge))?;
      Ok((edge_key, solution))
    })
    .collect::<Result<Vec<_>, Report>>()?;

  let unconverged = UnconvergedEdges {
    edges: optimized
      .iter()
      .filter(|(_, solution)| solution.convergence == Convergence::IterationLimit)
      .map(|&(edge_key, _)| edge_key)
      .collect(),
    optimized: optimized.len(),
  };
  lengths.extend(
    optimized
      .into_iter()
      .map(|(edge_key, solution)| (edge_key, solution.branch_length)),
  );

  for chain in chains {
    let indels = chain_indels[&chain.optimized];
    chain.restore(indels, &mut lengths);
  }

  branch_lengths.extend(lengths.into_iter().map(|(edge_key, length)| (edge_key, Some(length))));

  Ok(unconverged)
}

fn parse_optimizable_lengths(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<BTreeMap<GraphEdgeKey, f64>, Report> {
  graph
    .get_edges()
    .map(|edge_ref| {
      let edge_key = edge_ref.key();
      let branch_length = branch_lengths[&edge_key]
        .ok_or_else(|| make_report!("Cannot optimize edge {edge_key} with a missing branch length"))?;
      validate_branch_length_value(branch_length).wrap_err_with(|| format!("Cannot optimize edge {edge_key}"))?;
      Ok((edge_key, branch_length))
    })
    .collect()
}

struct SeriesChain {
  edges: Vec<GraphEdgeKey>,
  ratios: Vec<f64>,
  optimized: GraphEdgeKey,
}

impl SeriesChain {
  fn capture_all(graph: &Graph, lengths: &BTreeMap<GraphEdgeKey, f64>) -> Vec<Self> {
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
        let edge_lengths: Vec<f64> = edges.iter().map(|e| lengths[e]).collect();
        let total: f64 = edge_lengths.iter().sum();
        let ratios = if total > 0.0 {
          edge_lengths.iter().map(|bl| bl / total).collect()
        } else {
          vec![1.0 / edges.len().to_f64_exact(); edges.len()]
        };
        let longest = (0..edges.len()).fold(0, |best, i| if edge_lengths[i] > edge_lengths[best] { i } else { best });
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

  fn indel_term(&self, indels: &EdgeIndels, lengths: &BTreeMap<GraphEdgeKey, f64>) -> IndelTerm {
    let held_length = self.held().map(|e| lengths[&e]).sum();
    indels.chain_term(&self.edges, held_length)
  }

  fn restore(self, indels: IndelTerm, lengths: &mut BTreeMap<GraphEdgeKey, f64>) {
    let total = if lengths[&self.optimized] == 0.0 && indels.count == 0 {
      0.0
    } else {
      self.edges.iter().map(|e| lengths[e]).sum()
    };
    for (edge, ratio) in self.edges.into_iter().zip(self.ratios) {
      lengths.insert(edge, total * ratio);
    }
  }
}

pub(crate) fn guess_initial_branch_lengths(
  graph: &Graph,
  indels: &EdgeIndels,
  evidence: &InitialGuessEvidence,
  existing: ExistingBranchLengths,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  if evidence.total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot compute initial guess");
  }

  let one_mutation = one_mutation(evidence.total_length);

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

    let sub_count: usize = evidence.sub_counts[&edge_key];

    let effective_length: usize = evidence.effective_lengths[&edge_key];

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

fn optimize_edge(contributions: &[OptimizationContribution], edge: &EdgeOptimization) -> Result<EdgeSolution, Report> {
  let EdgeOptimization {
    branch_length,
    indels,
    one_mutation,
    method,
  } = *edge;

  let zero_length = ZeroLength::of_edge(contributions, indels);

  if indels.count == 0 && is_zero_branch_optimal(contributions, zero_length) {
    return Ok(EdgeSolution {
      branch_length: 0.0,
      convergence: Convergence::Converged,
    });
  }

  let start = start_branch_length(branch_length, indels, zero_length, one_mutation);
  let min_branch_length = min_branch_length_for_indels(indels, one_mutation);

  let solution = match method {
    BranchOptMethod::Brent => solve_brent(start, contributions, indels, min_branch_length, one_mutation),
    BranchOptMethod::Newton => solve_newton(start, contributions, indels, min_branch_length, one_mutation),
  }?;

  Ok(EdgeSolution {
    branch_length: reconcile_zero_boundary(
      solution.branch_length,
      start,
      contributions,
      indels,
      zero_length,
      one_mutation,
    )?,
    ..solution
  })
}

#[derive(Clone, Copy)]
struct EdgeOptimization {
  branch_length: f64,
  indels: IndelTerm,
  one_mutation: f64,
  method: BranchOptMethod,
}

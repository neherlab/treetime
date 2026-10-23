use crate::optimize::branch_length::validate_branch_length_value;
use crate::optimize::gather::{gather_alignment_edge_indel_counts, gather_edge_indel_counts, sum_alignment_counts};
use crate::optimize::likelihood::OptimizationMetrics;
use crate::partition::set::PartitionSet;
use eyre::Report;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use statrs::function::factorial::ln_factorial;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_primitives::LogLh;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::cast::usize_to_u64;
use treetime_utils::{make_error, make_report};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IndelPolicy {
  Include,
  Exclude,
}

impl IndelPolicy {
  pub(crate) fn edge_indels(
    self,
    graph: &Graph,
    partitions: &PartitionSet,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<EdgeIndels, Report> {
    match self {
      Self::Include => {
        let alignment_counts = gather_alignment_edge_indel_counts(graph, partitions)?;
        let alignment_rates = alignment_counts
          .iter()
          .map(|(&alignment, counts)| (alignment, estimate_indel_rate(graph, counts, branch_lengths)))
          .collect();
        Ok(EdgeIndels::Included {
          counts: sum_alignment_counts(graph, &alignment_counts),
          alignment_rates,
        })
      },
      Self::Exclude => Ok(EdgeIndels::Excluded),
    }
  }
}

#[derive(Clone, Debug)]
pub(crate) enum EdgeIndels {
  Included {
    counts: BTreeMap<GraphEdgeKey, usize>,
    alignment_rates: BTreeMap<usize, f64>,
  },
  Excluded,
}

impl EdgeIndels {
  pub(crate) fn recount(&self, graph: &Graph, partitions: &PartitionSet) -> Result<Self, Report> {
    Ok(match self {
      Self::Included { alignment_rates, .. } => Self::Included {
        counts: gather_edge_indel_counts(graph, partitions)?,
        alignment_rates: alignment_rates.clone(),
      },
      Self::Excluded => Self::Excluded,
    })
  }

  pub(crate) fn divide_rates(&mut self, divisor: f64) {
    if let Self::Included { alignment_rates, .. } = self {
      for rate in alignment_rates.values_mut() {
        *rate /= divisor;
      }
    }
  }

  pub(crate) fn alignment_rates(&self) -> BTreeMap<usize, f64> {
    match self {
      Self::Included { alignment_rates, .. } => alignment_rates.clone(),
      Self::Excluded => BTreeMap::new(),
    }
  }

  pub(crate) fn rate(&self) -> f64 {
    match self {
      Self::Included { alignment_rates, .. } => alignment_rates.values().sum(),
      Self::Excluded => 0.0,
    }
  }

  pub(crate) fn term(&self, edge_key: GraphEdgeKey) -> IndelTerm {
    match self {
      Self::Included { counts, .. } => IndelTerm::new(counts[&edge_key], self.rate()),
      Self::Excluded => IndelTerm::new(0, 0.0),
    }
  }

  pub(crate) fn chain_term(&self, edges: &[GraphEdgeKey], held_length: f64) -> IndelTerm {
    IndelTerm {
      count: edges.iter().map(|&edge_key| self.term(edge_key).count).sum(),
      rate: self.rate(),
      held_length,
    }
  }

  pub(crate) fn total_log_lh(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<LogLh, Report> {
    match self {
      Self::Included { counts, .. } => total_indel_log_lh(graph, counts, branch_lengths, self.rate()),
      Self::Excluded => Ok(LogLh::ZERO),
    }
  }

  pub(crate) fn has_indel_edge_with_zero_length(
    &self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> bool {
    graph.get_edges().any(|edge_ref| {
      let edge_key = edge_ref.key();
      branch_lengths[&edge_key].unwrap_or(0.0) == 0.0 && self.term(edge_key).requires_positive_length()
    })
  }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct IndelTerm {
  pub(crate) count: usize,
  pub(crate) rate: f64,
  pub(crate) held_length: f64,
}

impl IndelTerm {
  pub(crate) fn new(count: usize, rate: f64) -> Self {
    Self {
      count,
      rate,
      held_length: 0.0,
    }
  }

  pub(crate) fn metrics(&self, branch_length: f64) -> Result<OptimizationMetrics, Report> {
    poisson_indel_log_lh(self.count, self.rate, branch_length + self.held_length)
  }

  pub(crate) fn requires_positive_length(&self) -> bool {
    self.count > 0 && self.held_length == 0.0
  }

  pub(crate) fn length_estimate(&self) -> Option<f64> {
    (self.count > 0 && self.rate > 0.0).then(|| self.count.to_f64_exact() / self.rate)
  }
}

fn estimate_indel_rate(
  graph: &Graph,
  indel_counts: &BTreeMap<GraphEdgeKey, usize>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> f64 {
  let per_edge = graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| {
      let edge_key = edge_ref.key();
      let branch_length = branch_lengths[&edge_key].unwrap_or(0.0);
      let edge_indels = indel_counts[&edge_key];
      (edge_indels, branch_length)
    })
    .collect::<Vec<_>>();
  let total_indels = per_edge.iter().map(|(indels, _)| indels).sum::<usize>();
  let total_branch_length = per_edge.iter().map(|(_, branch_length)| branch_length).sum::<f64>();

  if total_branch_length > 0.0 && total_indels > 0 {
    total_indels.to_f64_exact() / total_branch_length
  } else {
    0.0
  }
}

fn total_indel_log_lh(
  graph: &Graph,
  indel_counts: &BTreeMap<GraphEdgeKey, usize>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  indel_rate: f64,
) -> Result<LogLh, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| -> Result<LogLh, Report> {
      let edge_key = edge_ref.key();
      let branch_length = branch_lengths[&edge_key].ok_or_else(|| {
        make_report!("Cannot evaluate indel likelihood for edge {edge_key} with a missing branch length")
      })?;
      validate_branch_length_value(branch_length)?;
      let indel_count: usize = indel_counts[&edge_key];
      if indel_count > 0 && branch_length <= 0.0 {
        Ok(LogLh::IMPOSSIBLE)
      } else {
        Ok(poisson_indel_log_lh(indel_count, indel_rate, branch_length)?.log_lh)
      }
    })
    .collect::<Vec<_>>()
    .into_iter()
    .sum()
}

pub(crate) fn poisson_indel_log_lh(k: usize, mu: f64, t: f64) -> Result<OptimizationMetrics, Report> {
  validate_branch_length_value(t)?;

  if mu == 0.0 {
    return Ok(if k == 0 {
      OptimizationMetrics::default()
    } else {
      OptimizationMetrics::new(LogLh::IMPOSSIBLE, 0.0, 0.0)
    });
  }

  if k == 0 {
    return Ok(OptimizationMetrics::new(LogLh::new(-mu * t), -mu, 0.0));
  }

  if t == 0.0 {
    return make_error!("Poisson indel likelihood requires a positive branch length when k > 0, got t={t}");
  }

  let k_f = k.to_f64_exact();
  let lambda = mu * t;
  let log_lh = k_f * lambda.ln() - lambda - ln_factorial(usize_to_u64(k));
  let derivative = k_f / t - mu;
  let second_derivative = -k_f / (t * t);

  Ok(OptimizationMetrics::new(
    LogLh::new(log_lh),
    derivative,
    second_derivative,
  ))
}

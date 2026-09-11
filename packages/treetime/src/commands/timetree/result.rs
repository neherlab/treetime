use crate::clock::clock_model::ClockModel;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::partition::timetree::partition::GraphTimetree;
use crate::seq::mutation::Mutation;
use crate::timetree::confidence::NodeConfidenceInterval;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_io::dates_csv::DatesMap;
use treetime_primitives::Seq;

#[derive(Serialize)]
pub struct TimetreeGraphData {
  pub clock_model: ClockModel,
  pub confidence_intervals: Option<Vec<NodeConfidenceInterval>>,
  pub dates: Option<DatesMap>,
  pub gtr: Option<GTR>,
  pub model_name: Option<GtrModelName>,
  pub mutation_counts: Option<BTreeMap<GraphEdgeKey, usize>>,
}

impl TimetreeGraphData {
  pub fn new(
    clock_model: ClockModel,
    confidence_intervals: Option<Vec<NodeConfidenceInterval>>,
    dates: Option<DatesMap>,
    gtr: Option<GTR>,
    model_name: Option<GtrModelName>,
    mutation_counts: Option<BTreeMap<GraphEdgeKey, usize>>,
  ) -> Self {
    Self {
      clock_model,
      confidence_intervals,
      dates,
      gtr,
      model_name,
      mutation_counts,
    }
  }
}

/// Nucleotide sequences and mutations gathered from the timetree partition for the tree writers.
///
/// Gathered once, serially, from the partition while it is in scope in the command, so the auspice,
/// phyloxml, MAT, and Newick-comment writers read plain value maps instead of reading the partition
/// during serialization.
#[derive(Debug, Default)]
pub struct TimetreeOutputMaps {
  /// Reconstructed nucleotide root sequence, or `None` when no partition exists.
  pub root_sequence: Option<Seq>,
  /// Reconstructed nucleotide sequence per node, for phyloxml clades.
  pub node_sequences: BTreeMap<GraphNodeKey, Seq>,
  /// Nucleotide mutations (substitutions followed by indels) per edge.
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}

/// Per-node timetree output as a value.
///
/// Holds the durable per-node results the timetree output writers consume: the estimated `time`
/// (numerical date), the cumulative divergence `div`, the exclusion flags, the input branch support
/// `confidence`, and the three-rate `rate_susceptibility_dates`. `name` and `desc` are carried for
/// writers that key by name. The values are gathered from the tree after the pipeline and topology
/// ordering complete, so the map is keyed by the final node set.
#[derive(Debug, Clone, Serialize)]
pub struct TimetreeNodeOut {
  pub name: Option<String>,
  pub desc: Option<String>,
  pub confidence: Option<f64>,
  pub time: Option<f64>,
  pub div: f64,
  pub is_outlier: bool,
  pub bad_branch: bool,
  pub rate_susceptibility_dates: Option<[f64; 3]>,
}

/// Per-edge timetree output as a value: the branch lengths and the relaxed-clock rate multiplier the
/// output writers read.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct TimetreeEdgeOut {
  pub branch_length: Option<f64>,
  pub time_length: Option<f64>,
  pub clock_branch_length: Option<f64>,
  pub gamma: f64,
}

impl TimetreeEdgeOut {
  /// The clock-constrained length once it has been committed, the ML or input length before that and
  /// after a topology change has invalidated it.
  pub fn profile_branch_length(&self) -> Option<f64> {
    self.clock_branch_length.or(self.branch_length)
  }
}

#[derive(Serialize)]
pub struct TimetreeResult {
  #[serde(skip)]
  pub graph: GraphTimetree<TimetreeGraphData>,
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, TimetreeNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, TimetreeEdgeOut>,
}

impl std::ops::Deref for TimetreeResult {
  type Target = TimetreeGraphData;

  fn deref(&self) -> &Self::Target {
    self.graph.data()
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_result_profile_branch_length_prefers_clock_length() {
    let edge = TimetreeEdgeOut {
      branch_length: Some(0.5),
      time_length: Some(2.0),
      clock_branch_length: Some(0.25),
      gamma: 1.0,
    };
    assert_eq!(Some(0.25), edge.profile_branch_length());
  }

  #[test]
  fn test_result_profile_branch_length_falls_back_to_branch_length() {
    let edge = TimetreeEdgeOut {
      branch_length: Some(0.5),
      time_length: Some(2.0),
      clock_branch_length: None,
      gamma: 1.0,
    };
    assert_eq!(Some(0.5), edge.profile_branch_length());
  }

  #[test]
  fn test_result_profile_branch_length_none_when_both_absent() {
    let edge = TimetreeEdgeOut {
      branch_length: None,
      time_length: None,
      clock_branch_length: None,
      gamma: 1.0,
    };
    assert_eq!(None, edge.profile_branch_length());
  }
}

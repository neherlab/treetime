use crate::partition::marginal::discrete::partition::PartitionMarginalDiscrete;
use crate::payload::ancestral::GraphAncestral;
use indexmap::IndexMap;
use itertools::Itertools;
use ndarray::Array1;
use serde::Serialize;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::{GraphNodeKey, Named};
use treetime_primitives::LogLh;
#[derive(Clone, Debug, Serialize)]
pub struct ConfidenceRow {
  /// Node name.
  pub node: String,
  /// Probability for each state (in state order).
  #[serde(serialize_with = "treetime_utils::array::serde::array1_as_vec")]
  pub profile: Array1<f64>,
}

/// Structured confidence output for all nodes.
#[derive(Clone, Debug, Serialize)]
pub struct MugrationConfidenceOutput {
  /// State names in order (column headers).
  pub states: Vec<String>,
  /// Confidence rows for each node.
  pub rows: Vec<ConfidenceRow>,
}

impl MugrationConfidenceOutput {
  pub fn new(graph: &GraphAncestral, partition: &PartitionMarginalDiscrete) -> Self {
    let states: Vec<String> = partition.states.iter().map(|s| s.to_owned()).collect();

    let rows: Vec<ConfidenceRow> = graph
      .get_nodes()
      .iter()
      .filter_map(|node| {
        let node_guard = node.read_arc();
        let node_key = node_guard.key();
        let payload = node_guard.payload().read_arc();
        let node_name = payload
          .name()
          .map_or_else(|| format!("node_{}", node_key.0), |n| n.as_ref().to_owned());

        partition.get_confidence(node_key).map(|profile| ConfidenceRow {
          node: node_name,
          profile,
        })
      })
      .collect();

    Self { states, rows }
  }

  /// Convert to map format for test comparison.
  pub fn to_map(&self) -> BTreeMap<String, Vec<String>> {
    self
      .rows
      .iter()
      .map(|row| {
        let formatted: Vec<String> = row.profile.iter().map(|p| format!("{p:.6}")).collect();
        (row.node.clone(), formatted)
      })
      .collect()
  }

  /// Render as CSV content.
  pub fn render_csv(&self) -> String {
    let mut out = String::new();
    // Header
    writeln!(out, "node,{}", self.states.join(",")).unwrap();
    // Data rows
    for row in &self.rows {
      let probs = row.profile.iter().map(|p| format!("{p:.6}")).join(",");
      writeln!(out, "{},{probs}", row.node).unwrap();
    }
    out
  }
}

/// Structured trait output for all nodes.
#[derive(Clone, Debug, Serialize)]
pub struct MugrationTraitsOutput {
  /// Attribute name (column header).
  pub attribute: String,
  /// Trait assignments keyed by node name (insertion order preserved).
  pub assignments: IndexMap<String, String>,
}

impl MugrationTraitsOutput {
  pub fn new(attribute: &str, assignments: IndexMap<String, String>) -> Self {
    Self {
      attribute: attribute.to_owned(),
      assignments,
    }
  }

  /// Render as CSV content.
  pub fn render_csv(&self) -> String {
    let mut out = String::new();
    // Header
    writeln!(out, "node,{}", self.attribute).unwrap();
    // Data rows (IndexMap preserves insertion order)
    for (node, trait_value) in &self.assignments {
      writeln!(out, "{node},{trait_value}").unwrap();
    }
    out
  }
}

#[derive(Debug, Serialize)]
pub struct MugrationGraphData {
  pub traits: MugrationTraitsOutput,
  pub confidence: MugrationConfidenceOutput,
  pub log_lh: LogLh,
  pub partition: Arc<PartitionMarginalDiscrete>,
}

/// Per-node mugration output as a value: the name and input branch support the output writers read.
/// Trait assignments and confidence stay sourced from the discrete partition.
#[derive(Debug, Clone, Serialize)]
pub struct MugrationNodeOut {
  pub name: Option<String>,
  pub confidence: Option<f64>,
}

/// Per-edge mugration output as a value: the branch length the output writers read.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct EdgeOut {
  pub branch_length: Option<f64>,
}

/// Mugration result as a value.
///
/// The discrete inference result (`discrete`), the reconstructed attribute name, and the marginal
/// log likelihood are reachable directly off the result. `discrete` shares its inference partition
/// with the copy the output writers still read from the graph, so no per-node profile data is
/// duplicated. `graph` carries the tree and the metadata the writers read; the writers move onto the
/// result value in a later step, after which `graph` and the partition-in-graph go away.
#[derive(Debug, serde::Serialize)]
pub struct MugrationResult {
  #[serde(skip)]
  pub graph: GraphAncestral<MugrationGraphData>,
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, MugrationNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
  #[serde(skip)]
  pub discrete: Arc<PartitionMarginalDiscrete>,
  #[serde(skip)]
  pub attribute: String,
  #[serde(skip)]
  pub log_lh: LogLh,
}

impl std::ops::Deref for MugrationResult {
  type Target = MugrationGraphData;

  fn deref(&self) -> &Self::Target {
    self.graph.data()
  }
}

impl MugrationResult {
  pub fn new(graph: GraphAncestral, partition: PartitionMarginalDiscrete, attribute: &str, log_lh: LogLh) -> Self {
    let assignments = extract_trait_assignments(&graph, &partition);
    let traits = MugrationTraitsOutput::new(attribute, assignments);
    let confidence = MugrationConfidenceOutput::new(&graph, &partition);

    // Gather the per-node name/confidence and per-edge branch length off the tree into keyed value
    // maps the output writers consume. Trait assignments and entropy stay sourced from the discrete
    // partition; only the name, input-branch-support, and branch-length reads move off the payload.
    // Topology ordering reorders children only, so these maps match what the post-order writers read.
    let nodes: BTreeMap<GraphNodeKey, MugrationNodeOut> = graph
      .get_nodes()
      .iter()
      .map(|node| {
        let node = node.read_arc();
        let payload = node.payload().read_arc();
        (
          node.key(),
          MugrationNodeOut {
            name: payload.name.clone(),
            confidence: payload.confidence,
          },
        )
      })
      .collect();
    let edges: BTreeMap<GraphEdgeKey, EdgeOut> = graph
      .get_edges()
      .iter()
      .map(|edge| {
        let edge = edge.read_arc();
        let branch_length = edge.payload().read_arc().branch_length;
        (edge.key(), EdgeOut { branch_length })
      })
      .collect();

    let partition = Arc::new(partition);
    let data = MugrationGraphData {
      traits,
      confidence,
      log_lh,
      partition: Arc::clone(&partition),
    };
    Self {
      graph: graph.map_data(data),
      nodes,
      edges,
      discrete: partition,
      attribute: attribute.to_owned(),
      log_lh,
    }
  }

  pub fn trait_assignments(&self) -> &IndexMap<String, String> {
    &self.graph.data().traits.assignments
  }
}

fn extract_trait_assignments(
  graph: &GraphAncestral,
  partition: &PartitionMarginalDiscrete,
) -> IndexMap<String, String> {
  graph
    .get_nodes()
    .iter()
    .filter_map(|node| {
      let node_guard = node.read_arc();
      let node_key = node_guard.key();
      let payload = node_guard.payload().read_arc();
      let node_name = payload
        .name()
        .map_or_else(|| format!("node_{}", node_key.0), |n| n.as_ref().to_owned());

      partition
        .get_reconstructed_trait(node_key)
        .map(|trait_value| (node_name, trait_value))
    })
    .collect()
}

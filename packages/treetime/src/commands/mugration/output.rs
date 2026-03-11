use crate::representation::partition::discrete::PartitionDiscrete;
use crate::representation::payload::ancestral::GraphAncestral;
use indexmap::IndexMap;
use itertools::Itertools;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt::Write;
use treetime_graph::node::Named;
use treetime_utils::array::serde::{array1_as_vec, array1_from_vec};

/// GTR model summary for mugration output.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MugrationGtrOutput {
  /// Attribute being reconstructed.
  pub attribute: String,
  /// Number of discrete states.
  pub n_states: usize,
  /// State names in order.
  pub states: Vec<String>,
  /// Equilibrium frequencies.
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub pi: Array1<f64>,
  /// Substitution rate.
  pub mu: f64,
}

impl MugrationGtrOutput {
  pub fn new(partition: &PartitionDiscrete, attribute: &str) -> Self {
    Self {
      attribute: attribute.to_owned(),
      n_states: partition.n_states(),
      states: partition.states.iter().map(|s| s.to_owned()).collect(),
      pi: partition.gtr.pi.clone(),
      mu: partition.gtr.mu,
    }
  }
}

/// Confidence profile for a single node.
#[derive(Clone, Debug)]
pub struct ConfidenceRow {
  /// Node name.
  pub node: String,
  /// Probability for each state (in state order).
  pub profile: Array1<f64>,
}

/// Structured confidence output for all nodes.
#[derive(Clone, Debug)]
pub struct MugrationConfidenceOutput {
  /// State names in order (column headers).
  pub states: Vec<String>,
  /// Confidence rows for each node.
  pub rows: Vec<ConfidenceRow>,
}

impl MugrationConfidenceOutput {
  pub fn new(graph: &GraphAncestral, partition: &PartitionDiscrete) -> Self {
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
          profile: profile.clone(),
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
#[derive(Clone, Debug)]
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

/// Result of mugration execution.
///
/// Contains all structured results needed by output builders and tests.
#[derive(Debug)]
pub struct MugrationResult {
  /// GTR model summary.
  pub gtr: MugrationGtrOutput,
  /// Reconstructed trait assignments.
  pub traits: MugrationTraitsOutput,
  /// Confidence profiles for all nodes.
  pub confidence: MugrationConfidenceOutput,
  /// Total log likelihood.
  pub log_lh: f64,
  /// The graph with ancestral state data (for tree rendering).
  pub graph: GraphAncestral,
  /// The partition with discrete reconstruction (for tree rendering).
  pub partition: PartitionDiscrete,
}

impl MugrationResult {
  pub fn new(graph: GraphAncestral, partition: PartitionDiscrete, attribute: &str, log_lh: f64) -> Self {
    let gtr = MugrationGtrOutput::new(&partition, attribute);
    let assignments = extract_trait_assignments(&graph, &partition);
    let traits = MugrationTraitsOutput::new(attribute, assignments);
    let confidence = MugrationConfidenceOutput::new(&graph, &partition);

    Self {
      gtr,
      traits,
      confidence,
      log_lh,
      graph,
      partition,
    }
  }

  /// Get trait assignments as a map (for backward compatibility with tests).
  pub fn trait_assignments(&self) -> &IndexMap<String, String> {
    &self.traits.assignments
  }
}

fn extract_trait_assignments(graph: &GraphAncestral, partition: &PartitionDiscrete) -> IndexMap<String, String> {
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

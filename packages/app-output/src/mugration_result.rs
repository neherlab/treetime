use indexmap::IndexMap;
use itertools::Itertools;
use ndarray::Array1;
use serde::Serialize;
use std::collections::BTreeMap;
use std::fmt::Write;
use treetime::mugration::pipeline::MugrationOutput;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;

// Serializable projection of the mugration core `MugrationOutput` for the application output layer.
//
// The core owns the one canonical `MugrationOutput` value; these types are the per-file, per-wire
// projection the writers serialize. Each is derived from the core output plus the parse-time value
// maps (names, input-tree branch support, branch lengths) the writers key by.

/// One node's confidence profile row.
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
  pub fn new(output: &MugrationOutput, names: &BTreeMap<GraphNodeKey, Option<String>>) -> Self {
    let states: Vec<String> = output.states.iter().map(|s| s.to_owned()).collect();

    let rows: Vec<ConfidenceRow> = output
      .graph
      .get_nodes()
      .filter_map(|node| {
        let node_key = node.key();
        let node_name = node_name_or_fallback(names, node_key);

        output.confidences[&node_key].clone().map(|profile| ConfidenceRow {
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

  #[allow(
    clippy::unwrap_used,
    reason = "unwrap on a value an upstream invariant guarantees is present"
  )]
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

  #[allow(
    clippy::unwrap_used,
    reason = "unwrap on a value an upstream invariant guarantees is present"
  )]
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

/// Per-node mugration output as a value: the name and input branch support the output writers read.
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

/// Mugration result projection as a value.
///
/// `nodes`/`edges` carry the per-node and per-edge metadata the writers read; `traits` and
/// `confidence` carry the reconstructed discrete-trait projections. The tree topology and the
/// reconstructed value maps stay in the core [`MugrationOutput`].
#[derive(Debug, Serialize)]
pub struct MugrationResult {
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, MugrationNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
  #[serde(skip)]
  pub traits: MugrationTraitsOutput,
  #[serde(skip)]
  pub confidence: MugrationConfidenceOutput,
}

impl MugrationResult {
  pub fn new(
    output: &MugrationOutput,
    input_confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    attribute: &str,
  ) -> Self {
    // Gather the per-node name/input-branch-support and per-edge branch length off the tree into keyed
    // value maps the output writers consume. The name comes from the passed name map, and
    // input-branch-support and branch length come from the passed value maps. Topology ordering
    // reorders children only, so these maps match what the post-order writers read.
    let nodes: BTreeMap<GraphNodeKey, MugrationNodeOut> = output
      .graph
      .get_nodes()
      .map(|node| {
        let key = node.key();
        (
          key,
          MugrationNodeOut {
            name: names.get(&key).cloned().flatten(),
            confidence: input_confidences.get(&key).copied().flatten(),
          },
        )
      })
      .collect();
    let assignments = extract_trait_assignments(output, names);
    let traits = MugrationTraitsOutput::new(attribute, assignments);
    let confidence = MugrationConfidenceOutput::new(output, names);

    let edges: BTreeMap<GraphEdgeKey, EdgeOut> = output
      .graph
      .get_edges()
      .map(|edge| {
        let key = edge.key();
        let branch_length = branch_lengths.get(&key).copied().flatten();
        (key, EdgeOut { branch_length })
      })
      .collect();

    Self {
      nodes,
      edges,
      traits,
      confidence,
    }
  }

  pub fn trait_assignments(&self) -> &IndexMap<String, String> {
    &self.traits.assignments
  }
}

/// Gather the reconstructed trait assignments keyed by node name, sourced from the core output's
/// per-node reconstructed-trait map. The reads are keyed by node key and independent of node ordering.
fn extract_trait_assignments(
  output: &MugrationOutput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> IndexMap<String, String> {
  output
    .graph
    .get_nodes()
    .filter_map(|node| {
      let node_key = node.key();
      let node_name = node_name_or_fallback(names, node_key);

      output.reconstructed_traits[&node_key]
        .clone()
        .map(|trait_value| (node_name, trait_value))
    })
    .collect()
}

/// Resolve a node's label from the threaded name map, falling back to `node_{key}` when unnamed.
fn node_name_or_fallback(names: &BTreeMap<GraphNodeKey, Option<String>>, node_key: GraphNodeKey) -> String {
  names[&node_key]
    .clone()
    .unwrap_or_else(|| format!("node_{}", node_key.0))
}

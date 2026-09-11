use crate::gtr::gtr::GTR;
use crate::partition::marginal::discrete::partition::PartitionMarginalDiscrete;
use crate::partition::storage::discrete::DiscreteStates;
use crate::payload::ancestral::GraphAncestral;
use indexmap::IndexMap;
use itertools::Itertools;
use ndarray::Array1;
use serde::Serialize;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
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
  pub fn new(
    graph: &GraphAncestral,
    partition: &PartitionMarginalDiscrete,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
  ) -> Self {
    let states: Vec<String> = partition.states.iter().map(|s| s.to_owned()).collect();

    let rows: Vec<ConfidenceRow> = graph
      .get_nodes()
      .iter()
      .filter_map(|node| {
        let node_key = node.read_arc().key();
        let node_name = node_name_or_fallback(names, node_key);

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

/// Discrete traits and confidence profiles gathered from the mugration partition for the output writers.
///
/// Gathered once, serially, from the discrete partition while it is in scope in the command, so the
/// auspice, phyloxml, Newick-comment, augur, and GTR writers read plain value maps instead of reading
/// the partition during serialization. The maps are keyed by node key over every node and carry the raw
/// per-key profile, so `build_confidence_map`/`compute_entropy` reproduce the current output exactly.
#[derive(Debug)]
pub struct MugrationOutputMaps {
  /// Reconstructed discrete trait per node (argmax state name), or `None` when the node has no profile.
  pub reconstructed_traits: BTreeMap<GraphNodeKey, Option<String>>,
  /// Confidence profile per node (raw, unfiltered), or `None` when the node has no profile.
  pub confidences: BTreeMap<GraphNodeKey, Option<Array1<f64>>>,
  /// Discrete state names in order.
  pub states: DiscreteStates,
  /// The inferred discrete GTR model.
  pub gtr: GTR,
  /// Number of real states (excludes the missing-data marker).
  pub n_states: usize,
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
  pub fn new(
    graph: GraphAncestral,
    confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    partition: PartitionMarginalDiscrete,
    attribute: &str,
    log_lh: LogLh,
  ) -> Self {
    // Gather the per-node name/confidence and per-edge branch length off the tree into keyed value
    // maps the output writers consume. Trait assignments and entropy stay sourced from the discrete
    // partition; the name comes from the passed name map, and input-branch-support and branch-length
    // reads move off the payload. Topology ordering reorders children only, so these maps match what
    // the post-order writers read.
    let nodes: BTreeMap<GraphNodeKey, MugrationNodeOut> = graph
      .get_nodes()
      .iter()
      .map(|node| {
        let key = node.read_arc().key();
        (
          key,
          MugrationNodeOut {
            name: names.get(&key).cloned().flatten(),
            confidence: confidences.get(&key).copied().flatten(),
          },
        )
      })
      .collect();
    let assignments = extract_trait_assignments(&graph, &partition, names);
    let traits = MugrationTraitsOutput::new(attribute, assignments);
    let confidence = MugrationConfidenceOutput::new(&graph, &partition, names);

    let edges: BTreeMap<GraphEdgeKey, EdgeOut> = graph
      .get_edges()
      .iter()
      .map(|edge| {
        let key = edge.read_arc().key();
        let branch_length = branch_lengths.get(&key).copied().flatten();
        (key, EdgeOut { branch_length })
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
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> IndexMap<String, String> {
  graph
    .get_nodes()
    .iter()
    .filter_map(|node| {
      let node_key = node.read_arc().key();
      let node_name = node_name_or_fallback(names, node_key);

      partition
        .get_reconstructed_trait(node_key)
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

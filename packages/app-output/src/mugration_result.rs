use indexmap::IndexMap;
use itertools::Itertools;
use ndarray::Array1;
use serde::Serialize;
use std::collections::BTreeMap;
use std::fmt::Write;
use treetime::mugration::pipeline::MugrationOutput;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;

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

#[derive(Clone, Debug, Serialize)]
pub struct MugrationConfidenceOutput {
  states: Vec<String>,
  rows: Vec<ConfidenceRow>,
}

impl MugrationConfidenceOutput {
  fn new(output: &MugrationOutput, names: &BTreeMap<GraphNodeKey, Option<String>>) -> Self {
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
  pub fn render_csv(&self) -> String {
    let mut out = String::new();
    writeln!(out, "node,{}", self.states.join(",")).unwrap();
    for row in &self.rows {
      let probs = row.profile.iter().map(|p| format!("{p:.6}")).join(",");
      writeln!(out, "{},{probs}", row.node).unwrap();
    }
    out
  }
}

#[derive(Clone, Debug, Serialize)]
pub struct ConfidenceRow {
  node: String,
  #[serde(serialize_with = "treetime_utils::array::serde::array1_as_vec")]
  profile: Array1<f64>,
}

#[derive(Clone, Debug, Serialize)]
pub struct MugrationTraitsOutput {
  pub(crate) attribute: String,
  pub(crate) assignments: IndexMap<String, String>,
}

impl MugrationTraitsOutput {
  fn new(attribute: &str, assignments: IndexMap<String, String>) -> Self {
    Self {
      attribute: attribute.to_owned(),
      assignments,
    }
  }

  #[allow(
    clippy::unwrap_used,
    reason = "unwrap on a value an upstream invariant guarantees is present"
  )]
  pub fn render_csv(&self) -> String {
    let mut out = String::new();
    writeln!(out, "node,{}", self.attribute).unwrap();
    for (node, trait_value) in &self.assignments {
      writeln!(out, "{node},{trait_value}").unwrap();
    }
    out
  }
}

#[derive(Debug, Clone, Serialize)]
pub struct MugrationNodeOut {
  pub(crate) name: Option<String>,
  pub(crate) confidence: Option<f64>,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct EdgeOut {
  branch_length: Option<f64>,
}

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

fn node_name_or_fallback(names: &BTreeMap<GraphNodeKey, Option<String>>, node_key: GraphNodeKey) -> String {
  names[&node_key]
    .clone()
    .unwrap_or_else(|| format!("node_{}", node_key.0))
}

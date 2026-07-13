//! Auspice v2 JSON adapter for the format-neutral TreeIR graph.
//!
//! Output semantics follow augur's `export_v2.py`:
//! - `node_attrs.div` is cumulative divergence from the root (root = 0), precomputed
//!   on each IR node, rounded with [`format_number`] at precision 6.
//! - `node_attrs.num_date` carries the numeric date with optional `[lower, upper]`
//!   confidence, rounded at precision 3.
//! - `meta.colorings` are derived from which data classes the dataset populates:
//!   `num_date` (continuous) when dates are present, `bad_branch` (categorical) when
//!   the outlier flag is present, one categorical coloring per discrete trait, and
//!   the genotype coloring `gt` when branch mutations are present.
//! - `branch_attrs.mutations` groups per-branch substitutions by gene track.

use crate::auspice::{AuspiceGraphContext, AuspiceRead, AuspiceTreeContext, AuspiceWrite};
use crate::auspice_types::{
  AuspiceColoring, AuspiceDisplayDefaults, AuspiceNumDate, AuspiceTree, AuspiceTreeBranchAttrs, AuspiceTreeData,
  AuspiceTreeMeta, AuspiceTreeNode, AuspiceTreeNodeAttr, AuspiceTreeNodeAttrs,
};
use crate::tree_ir::mutation::TreeIrSub;
use crate::tree_ir::types::{TreeIrData, TreeIrEdge, TreeIrNode, TreeIrTrait};
use eyre::Report;
use serde_json::{Value, json};
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::Named;
use treetime_utils::make_error;

const COLORING_GENOTYPE: &str = "gt";
const COLORING_NUM_DATE: &str = "num_date";
const COLORING_BAD_BRANCH: &str = "bad_branch";

/// Round a float the way augur's `format_number` does: keep `precision` significant
/// digits in the fractional part while preserving all integer digits.
///
/// Reimplements `augur.export_v2.format_number`: total significant figures equal the
/// number of integer digits plus `precision`. Formatting through scientific notation
/// and parsing back yields the same `f64` augur stores after `float(f"{n:.{k}g}")`.
pub fn format_number(n: f64, precision: i32) -> f64 {
  if n == 0.0 || !n.is_finite() {
    return n;
  }
  let integral = n.abs().trunc();
  let significand = if integral >= 1.0 {
    integral.log10().floor() as i32 + 1
  } else {
    0
  };
  let sig_figs = (significand + precision).max(1) as usize;
  format!("{n:.*e}", sig_figs - 1)
    .parse()
    .expect("a float formatted in scientific notation must parse back")
}

/// Converter writing a TreeIR graph to Auspice v2 JSON.
pub struct TreeIrAuspiceWriter {
  has_bad_branch: bool,
}

impl AuspiceWrite<TreeIrNode, TreeIrEdge, TreeIrData> for TreeIrAuspiceWriter {
  fn new(graph: &Graph<TreeIrNode, TreeIrEdge, TreeIrData>) -> Result<Self, Report> {
    let has_bad_branch = graph.data().read_arc().has_bad_branch;
    Ok(Self { has_bad_branch })
  }

  fn auspice_data_from_graph_data(
    &self,
    graph: &Graph<TreeIrNode, TreeIrEdge, TreeIrData>,
  ) -> Result<AuspiceTreeData, Report> {
    let data = graph.data().read_arc();

    let mut colorings = vec![];
    if data.has_dates {
      colorings.push(coloring(COLORING_NUM_DATE, "Date", "continuous"));
    }
    if data.has_bad_branch {
      colorings.push(coloring(COLORING_BAD_BRANCH, "Excluded", "categorical"));
    }
    for attr in &data.trait_attrs {
      colorings.push(coloring(attr, attr, "categorical"));
    }
    if data.has_mutations {
      colorings.push(coloring(COLORING_GENOTYPE, "Genotype", "categorical"));
    }

    // Default coloring: a discrete trait if present, else the outlier flag, else date.
    let color_by = data
      .trait_attrs
      .first()
      .cloned()
      .or_else(|| data.has_bad_branch.then(|| COLORING_BAD_BRANCH.to_owned()))
      .or_else(|| data.has_dates.then(|| COLORING_NUM_DATE.to_owned()));

    let mut filters = data.trait_attrs.clone();
    if data.has_bad_branch {
      filters.push(COLORING_BAD_BRANCH.to_owned());
    }

    let root_sequence = (!data.root_sequence.is_empty()).then(|| data.root_sequence.clone());

    Ok(AuspiceTreeData {
      version: Some("v2".to_owned()),
      meta: AuspiceTreeMeta {
        title: data.title.clone().or_else(|| Some("TreeTime analysis".to_owned())),
        description: data.description.clone(),
        panels: vec!["tree".to_owned()],
        colorings,
        filters,
        display_defaults: AuspiceDisplayDefaults {
          color_by,
          ..AuspiceDisplayDefaults::default()
        },
        ..AuspiceTreeMeta::default()
      },
      root_sequence,
      other: Value::default(),
    })
  }

  fn auspice_node_from_graph_components(
    &mut self,
    context: &AuspiceGraphContext<TreeIrNode, TreeIrEdge, TreeIrData>,
  ) -> Result<AuspiceTreeNode, Report> {
    let &AuspiceGraphContext {
      node_key, node, edge, ..
    } = context;

    let name = node
      .name()
      .map_or_else(|| format!("node_{}", node_key.as_usize()), |n| n.as_ref().to_owned());

    let div = match node.div {
      Some(div) => {
        if !div.is_finite() {
          return make_error!("Node '{name}' has non-finite div={div}");
        }
        Some(format_number(div, 6))
      },
      None => None,
    };

    let num_date = match node.date {
      Some(date) => {
        if !date.is_finite() {
          return make_error!("Node '{name}' has non-finite date={date}");
        }
        let confidence = match node.date_confidence {
          Some([lo, hi]) => {
            if !lo.is_finite() || !hi.is_finite() {
              return make_error!("Node '{name}' has non-finite date confidence [{lo}, {hi}]");
            }
            Some([format_number(lo, 3), format_number(hi, 3)])
          },
          None => None,
        };
        Some(AuspiceNumDate {
          value: format_number(date, 3),
          confidence,
        })
      },
      None => None,
    };

    let bad_branch = self
      .has_bad_branch
      .then(|| AuspiceTreeNodeAttr::new(if node.bad_branch { "Yes" } else { "No" }));

    let other = build_trait_attrs(&node.traits);

    let mutations = edge.map_or_else(BTreeMap::new, |edge| group_mutations(&edge.mutations));

    Ok(AuspiceTreeNode {
      name,
      branch_attrs: AuspiceTreeBranchAttrs {
        mutations,
        labels: None,
        other: Value::default(),
      },
      node_attrs: AuspiceTreeNodeAttrs {
        div,
        num_date,
        bad_branch,
        clade_membership: None,
        region: None,
        country: None,
        division: None,
        other,
      },
      children: vec![],
      other: Value::default(),
    })
  }
}

/// Converter reading a TreeIR graph from Auspice v2 JSON.
pub struct TreeIrAuspiceReader {
  trait_attrs: Vec<String>,
}

impl AuspiceRead<TreeIrNode, TreeIrEdge, TreeIrData> for TreeIrAuspiceReader {
  fn new(tree: &AuspiceTree) -> Result<Self, Report> {
    Ok(Self {
      trait_attrs: trait_attrs_from_meta(&tree.data.meta),
    })
  }

  fn auspice_data_to_graph_data(&mut self, tree: &AuspiceTree) -> Result<TreeIrData, Report> {
    let meta = &tree.data.meta;
    let keys: Vec<&str> = meta.colorings.iter().map(|c| c.key.as_str()).collect();
    Ok(TreeIrData {
      title: meta.title.clone(),
      description: meta.description.clone(),
      root_sequence: tree.data.root_sequence.clone().unwrap_or_default(),
      trait_attrs: self.trait_attrs.clone(),
      has_dates: keys.contains(&COLORING_NUM_DATE),
      has_bad_branch: keys.contains(&COLORING_BAD_BRANCH),
      has_mutations: keys.contains(&COLORING_GENOTYPE),
    })
  }

  fn auspice_node_to_graph_components(
    &mut self,
    context: &AuspiceTreeContext,
  ) -> Result<(TreeIrNode, TreeIrEdge), Report> {
    let node = context.node;

    let date = node.node_attrs.num_date.as_ref().map(|nd| nd.value);
    let date_confidence = node.node_attrs.num_date.as_ref().and_then(|nd| nd.confidence);
    let bad_branch = node
      .node_attrs
      .bad_branch
      .as_ref()
      .is_some_and(|attr| attr.value == "Yes");

    let mut traits = BTreeMap::new();
    for attr in &self.trait_attrs {
      if let Some(value) = node.node_attrs.other.get(attr) {
        if let Some(trait_) = parse_trait(value) {
          traits.insert(attr.clone(), trait_);
        }
      }
    }

    let ir_node = TreeIrNode {
      name: Some(node.name.clone()),
      div: node.node_attrs.div,
      date,
      date_confidence,
      bad_branch,
      traits,
      ..TreeIrNode::default()
    };

    let mut mutations = vec![];
    for (gene, muts) in &node.branch_attrs.mutations {
      for m in muts {
        mutations.push(TreeIrSub::from_auspice_string(gene, m)?);
      }
    }

    let ir_edge = TreeIrEdge {
      branch_length: context.branch_length(),
      mutations,
      ..TreeIrEdge::default()
    };

    Ok((ir_node, ir_edge))
  }
}

fn coloring(key: &str, title: &str, type_: &str) -> AuspiceColoring {
  AuspiceColoring {
    key: key.to_owned(),
    title: title.to_owned(),
    type_: type_.to_owned(),
    ..AuspiceColoring::default()
  }
}

fn build_trait_attrs(traits: &BTreeMap<String, TreeIrTrait>) -> Value {
  if traits.is_empty() {
    return Value::default();
  }
  let mut map = serde_json::Map::new();
  for (attr, trait_) in traits {
    let mut obj = serde_json::Map::new();
    obj.insert("value".to_owned(), json!(trait_.value));
    if !trait_.confidence.is_empty() {
      let conf: BTreeMap<String, f64> = trait_
        .confidence
        .iter()
        .map(|(k, v)| (k.clone(), format_number(*v, 3)))
        .collect();
      obj.insert("confidence".to_owned(), json!(conf));
    }
    if let Some(entropy) = trait_.entropy {
      obj.insert("entropy".to_owned(), json!(format_number(entropy, 3)));
    }
    map.insert(attr.clone(), Value::Object(obj));
  }
  Value::Object(map)
}

fn group_mutations(subs: &[TreeIrSub]) -> BTreeMap<String, Vec<String>> {
  let mut grouped: BTreeMap<String, Vec<String>> = BTreeMap::new();
  for sub in subs {
    grouped
      .entry(sub.gene.clone())
      .or_default()
      .push(sub.to_auspice_string());
  }
  grouped
}

fn trait_attrs_from_meta(meta: &AuspiceTreeMeta) -> Vec<String> {
  let excluded = [COLORING_GENOTYPE, COLORING_NUM_DATE, COLORING_BAD_BRANCH, "author"];
  meta
    .colorings
    .iter()
    .map(|c| c.key.clone())
    .filter(|key| !excluded.contains(&key.as_str()))
    .collect()
}

fn parse_trait(value: &Value) -> Option<TreeIrTrait> {
  let value_str = value.get("value")?.as_str()?.to_owned();
  let confidence = value
    .get("confidence")
    .and_then(Value::as_object)
    .map(|obj| {
      obj
        .iter()
        .filter_map(|(k, v)| v.as_f64().map(|f| (k.clone(), f)))
        .collect()
    })
    .unwrap_or_default();
  let entropy = value.get("entropy").and_then(Value::as_f64);
  Some(TreeIrTrait {
    value: value_str,
    confidence,
    entropy,
  })
}

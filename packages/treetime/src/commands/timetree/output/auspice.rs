use crate::commands::timetree::output::confidence::NodeConfidenceInterval;
use crate::representation::payload::timetree::{EdgeTimetree, NodeTimetree};
use eyre::{Report, WrapErr};
use log::info;
use maplit::btreemap;
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::Path;
use treetime_graph::edge::HasBranchLength;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNodeKey, Named};
use treetime_io::auspice::{AuspiceGraphContext, AuspiceWrite, auspice_from_graph_with};
use treetime_io::auspice_types::{
  AuspiceColoring, AuspiceDisplayDefaults, AuspiceNumDate, AuspiceTreeBranchAttrs, AuspiceTreeData, AuspiceTreeMeta,
  AuspiceTreeNode, AuspiceTreeNodeAttr, AuspiceTreeNodeAttrs,
};
use treetime_io::json::{JsonPretty, json_write_file};

/// Write auspice v2 JSON for timetree results.
///
/// Produces `auspice_tree.json` matching v0's `create_auspice_json()` output:
/// node dates (with optional confidence intervals), cumulative divergence,
/// bad branch status. Mutations are omitted due to the partition-layer access
/// gap (see `M-core-branch-mutations-no-unified-api`).
pub fn write_auspice_json(
  graph: &Graph<NodeTimetree, EdgeTimetree, ()>,
  confidence_intervals: Option<&[NodeConfidenceInterval]>,
  outdir: &Path,
) -> Result<(), Report> {
  let ci_map = confidence_intervals.map_or_else(BTreeMap::new, build_ci_map);
  let mut writer = TimetreeAuspiceWriter {
    divs: btreemap! {},
    ci_map,
  };
  let tree = auspice_from_graph_with(&mut writer, graph)?;
  let filepath = outdir.join("auspice_tree.json");
  json_write_file(&filepath, &tree, JsonPretty(true)).wrap_err("Failed to write auspice JSON")?;
  info!("Wrote auspice JSON to {filepath}", filepath = filepath.display());
  Ok(())
}

/// Auspice v2 JSON writer for timetree output.
///
/// Accumulates cumulative divergence during pre-order traversal
/// and looks up confidence intervals from a precomputed map.
struct TimetreeAuspiceWriter {
  divs: BTreeMap<GraphNodeKey, f64>,
  ci_map: BTreeMap<String, ConfidenceBounds>,
}

impl AuspiceWrite<NodeTimetree, EdgeTimetree, ()> for TimetreeAuspiceWriter {
  fn new(_graph: &Graph<NodeTimetree, EdgeTimetree, ()>) -> Result<Self, Report> {
    Ok(Self {
      divs: btreemap! {},
      ci_map: BTreeMap::new(),
    })
  }

  fn auspice_data_from_graph_data(
    &self,
    _graph: &Graph<NodeTimetree, EdgeTimetree, ()>,
  ) -> Result<AuspiceTreeData, Report> {
    let colorings = vec![
      AuspiceColoring {
        title: "Date".to_owned(),
        type_: "continuous".to_owned(),
        key: "num_date".to_owned(),
        ..AuspiceColoring::default()
      },
      AuspiceColoring {
        title: "Excluded".to_owned(),
        type_: "categorical".to_owned(),
        key: "bad_branch".to_owned(),
        ..AuspiceColoring::default()
      },
    ];

    Ok(AuspiceTreeData {
      version: Some("v2".to_owned()),
      meta: AuspiceTreeMeta {
        title: Some("TreeTime timetree analysis".to_owned()),
        updated: Some(chrono::Utc::now().format("%Y-%m-%d").to_string()),
        panels: vec!["tree".to_owned()],
        colorings,
        display_defaults: AuspiceDisplayDefaults {
          color_by: Some("num_date".to_owned()),
          distance_measure: Some("num_date".to_owned()),
          ..AuspiceDisplayDefaults::default()
        },
        filters: vec!["bad_branch".to_owned()],
        ..AuspiceTreeMeta::default()
      },
      root_sequence: None,
      other: Value::default(),
    })
  }

  fn auspice_node_from_graph_components(
    &mut self,
    context: &AuspiceGraphContext<NodeTimetree, EdgeTimetree, ()>,
  ) -> Result<AuspiceTreeNode, Report> {
    let &AuspiceGraphContext {
      node_key,
      node,
      parent_key,
      edge,
      ..
    } = context;

    let name = node
      .name()
      .map_or_else(|| format!("node_{}", node_key.as_usize()), |n| n.as_ref().to_owned());

    // Cumulative divergence: parent_div + branch_length
    let parent_div = parent_key.and_then(|pk| self.divs.get(&pk).copied()).unwrap_or(0.0);
    let branch_length = edge.map_or(0.0, |e| e.branch_length().unwrap_or(0.0));
    let div = parent_div + branch_length;
    self.divs.insert(node_key, div);

    // Numeric date with optional confidence interval
    let num_date = node.time.map(|t| {
      let confidence = self.ci_map.get(&name).map(|ci| [ci.lower, ci.upper]);
      AuspiceNumDate { value: t, confidence }
    });

    let bad_branch = Some(AuspiceTreeNodeAttr::new(if node.bad_branch { "Yes" } else { "No" }));

    Ok(AuspiceTreeNode {
      name,
      branch_attrs: AuspiceTreeBranchAttrs {
        mutations: BTreeMap::new(),
        labels: None,
        other: Value::default(),
      },
      node_attrs: AuspiceTreeNodeAttrs {
        div: Some(div),
        num_date,
        bad_branch,
        clade_membership: None,
        region: None,
        country: None,
        division: None,
        other: Value::default(),
      },
      children: vec![],
      other: Value::default(),
    })
  }
}

struct ConfidenceBounds {
  lower: f64,
  upper: f64,
}

fn build_ci_map(intervals: &[NodeConfidenceInterval]) -> BTreeMap<String, ConfidenceBounds> {
  intervals
    .iter()
    .map(|ci| {
      (
        ci.name.clone(),
        ConfidenceBounds {
          lower: ci.lower,
          upper: ci.upper,
        },
      )
    })
    .collect()
}

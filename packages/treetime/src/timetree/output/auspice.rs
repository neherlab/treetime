use crate::partition::traits::BranchTopology;
use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
use crate::timetree::confidence::NodeConfidenceInterval;
use eyre::Report;
use serde_json::Value;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNodeKey, Named};
use treetime_io::auspice::{AuspiceGraphContext, AuspiceWrite, auspice_from_graph_with};
use treetime_io::auspice_types::{
  AuspiceColoring, AuspiceDisplayDefaults, AuspiceNumDate, AuspiceTree, AuspiceTreeBranchAttrs, AuspiceTreeData,
  AuspiceTreeMeta, AuspiceTreeNode, AuspiceTreeNodeAttr, AuspiceTreeNodeAttrs,
};
use treetime_utils::make_error;

/// Build partial auspice v2 tree data from a timetree graph without file I/O.
///
/// Produces node dates (with optional confidence intervals), cumulative divergence,
/// and bad-branch status. Branch mutations, branch-support confidence, and genome
/// annotations are not yet included (see `kb/issues/N-timetree-auspice-json-incomplete.md`).
///
/// When `mutation_counts` is provided, cumulative divergence is computed by
/// accumulating per-edge mutation counts during traversal instead of reading
/// `node.div` (which stores cumulative subs/site from branch lengths).
pub fn build_timetree_auspice(
  graph: &Graph<NodeTimetree, EdgeTimetree, ()>,
  confidence_intervals: Option<&[NodeConfidenceInterval]>,
  mutation_counts: Option<&BTreeMap<GraphEdgeKey, usize>>,
) -> Result<AuspiceTree, Report> {
  let ci_map = confidence_intervals.map_or_else(BTreeMap::new, build_ci_map);
  let mut writer = TimetreeAuspiceWriter {
    ci_map,
    mutation_counts: mutation_counts.cloned(),
    accumulated_divs: BTreeMap::new(),
  };
  auspice_from_graph_with(&mut writer, graph)
}

struct TimetreeAuspiceWriter {
  ci_map: BTreeMap<GraphNodeKey, ConfidenceBounds>,
  mutation_counts: Option<BTreeMap<GraphEdgeKey, usize>>,
  accumulated_divs: BTreeMap<GraphNodeKey, f64>,
}

impl AuspiceWrite<NodeTimetree, EdgeTimetree, ()> for TimetreeAuspiceWriter {
  fn new(_graph: &Graph<NodeTimetree, EdgeTimetree, ()>) -> Result<Self, Report> {
    make_error!(
      "TimetreeAuspiceWriter must be constructed directly with a CI map. \
       Use build_timetree_auspice() or auspice_from_graph_with() instead of auspice_from_graph()."
    )
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
          color_by: Some("bad_branch".to_owned()),
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
      node_key, node, graph, ..
    } = context;

    let name = node
      .name()
      .map_or_else(|| format!("node_{}", node_key.as_usize()), |n| n.as_ref().to_owned());

    let div = match &self.mutation_counts {
      Some(counts) => {
        let parent_div = match graph.node_parent(node_key)? {
          Some((_parent_key, edge_key)) => {
            let parent_div = context
              .parent_key
              .and_then(|pk| self.accumulated_divs.get(&pk).copied())
              .unwrap_or(0.0);
            let edge_muts = counts.get(&edge_key).copied().unwrap_or_default();
            parent_div + edge_muts as f64
          },
          None => 0.0,
        };
        self.accumulated_divs.insert(node_key, parent_div);
        parent_div
      },
      None => node.div,
    };

    if !div.is_finite() {
      return make_error!("Node '{name}' has non-finite div={div}");
    }
    if let Some(t) = node.time {
      if !t.is_finite() {
        return make_error!("Node '{name}' has non-finite time={t}");
      }
    }

    let num_date = node.time.map(|t| {
      let confidence = self.ci_map.get(&node_key).map(|ci| {
        if !ci.lower.is_finite() || !ci.upper.is_finite() {
          log::warn!(
            "Node '{name}' has non-finite CI bounds: [{lower}, {upper}]",
            lower = ci.lower,
            upper = ci.upper
          );
        }
        debug_assert!(
          ci.lower <= ci.upper,
          "CI bounds inverted for node '{name}': [{lower}, {upper}]",
          lower = ci.lower,
          upper = ci.upper
        );
        [ci.lower, ci.upper]
      });
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

fn build_ci_map(intervals: &[NodeConfidenceInterval]) -> BTreeMap<GraphNodeKey, ConfidenceBounds> {
  intervals
    .iter()
    .map(|ci| {
      (
        ci.key,
        ConfidenceBounds {
          lower: ci.lower,
          upper: ci.upper,
        },
      )
    })
    .collect()
}

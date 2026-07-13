//! Projection of a timetree analysis result into the format-neutral TreeIR graph.
//!
//! TreeTime's timetree command writes Auspice v2 JSON through the shared TreeIR
//! output path. This module mirrors the timetree domain graph into a TreeIR graph,
//! carrying the node-level data Auspice represents: cumulative divergence, numeric
//! date with confidence interval, and the bad-branch (temporal outlier) flag.
//!
//! Divergence accumulation follows augur's `node_div`: when per-edge mutation counts
//! are supplied (divergence measured in mutations) the cumulative divergence is the
//! running sum of those counts from the root; otherwise the node's stored cumulative
//! substitutions-per-site value is used directly. The root has divergence zero.

use crate::partition::timetree::GraphTimetree;
use crate::timetree::confidence::NodeConfidenceInterval;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};
use treetime_graph::node::{GraphNodeKey, Named};
use treetime_io::graph::TreeIrGraph;
use treetime_io::tree_ir::types::{TreeIrData, TreeIrEdge, TreeIrNode};

/// Build a TreeIR graph from a timetree result for Auspice v2 output.
pub fn build_timetree_ir(
  graph: &GraphTimetree,
  confidence_intervals: Option<&[NodeConfidenceInterval]>,
  mutation_counts: Option<&BTreeMap<GraphEdgeKey, usize>>,
) -> Result<TreeIrGraph, Report> {
  let ci_map: BTreeMap<GraphNodeKey, [f64; 2]> = confidence_intervals
    .map(|cis| cis.iter().map(|ci| (ci.key, [ci.lower, ci.upper])).collect())
    .unwrap_or_default();

  let mut ir = TreeIrGraph::with_data(TreeIrData {
    title: Some("TreeTime timetree analysis".to_owned()),
    has_dates: true,
    has_bad_branch: true,
    ..TreeIrData::default()
  });

  let mut key_map: BTreeMap<GraphNodeKey, GraphNodeKey> = BTreeMap::new();
  let mut div_map: BTreeMap<GraphNodeKey, f64> = BTreeMap::new();

  graph.iter_depth_first_preorder_forward(|node| {
    let dkey = node.key;
    let parent = node.parent_keys.first().copied();

    let div = match mutation_counts {
      Some(counts) => match parent {
        Some((pkey, ekey)) => {
          div_map.get(&pkey).copied().unwrap_or(0.0) + *counts.get(&ekey).unwrap_or(&0) as f64
        },
        None => 0.0,
      },
      None => node.payload.div,
    };
    div_map.insert(dkey, div);

    let date = node.payload.time;
    let ir_node = TreeIrNode {
      name: node.payload.name().map(|n| n.as_ref().to_owned()),
      div: Some(div),
      date,
      date_confidence: date.and_then(|_| ci_map.get(&dkey).copied()),
      bad_branch: node.payload.bad_branch,
      ..TreeIrNode::default()
    };
    let ir_key = ir.add_node(ir_node);
    key_map.insert(dkey, ir_key);

    if let Some((pkey, _ekey)) = parent {
      let ir_parent = key_map[&pkey];
      let branch_length = node.parents.first().and_then(|(_, edge)| edge.read_arc().branch_length());
      ir.add_edge(ir_parent, ir_key, TreeIrEdge {
        branch_length,
        ..TreeIrEdge::default()
      })?;
    }

    Ok(())
  })?;

  ir.build()?;
  Ok(ir)
}

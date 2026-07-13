use crate::commands::shared::ir_projection::subs_to_ir;
use crate::partition::timetree::GraphTimetree;
use crate::partition::traits::PartitionBranchOps;
use crate::timetree::confidence::NodeConfidenceInterval;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};
use treetime_graph::node::{GraphNodeKey, Named};
use treetime_io::graph::TreeIrGraph;
use treetime_io::tree_ir::types::{TreeIrData, TreeIrEdge, TreeIrNode};

pub fn build_timetree_ir(
  graph: &GraphTimetree,
  confidence_intervals: Option<&[NodeConfidenceInterval]>,
  mutation_counts: Option<&BTreeMap<GraphEdgeKey, usize>>,
  partition: Option<&dyn PartitionBranchOps>,
) -> Result<TreeIrGraph, Report> {
  let ci_map: BTreeMap<GraphNodeKey, [f64; 2]> = confidence_intervals
    .map(|cis| cis.iter().map(|ci| (ci.key, [ci.lower, ci.upper])).collect())
    .unwrap_or_default();

  let mut ir = TreeIrGraph::with_data(TreeIrData {
    title: Some("TreeTime timetree analysis".to_owned()),
    has_dates: true,
    has_bad_branch: true,
    has_mutations: partition.is_some(),
    ..TreeIrData::default()
  });

  let mut key_map: BTreeMap<GraphNodeKey, GraphNodeKey> = BTreeMap::new();
  let mut div_map: BTreeMap<GraphNodeKey, f64> = BTreeMap::new();

  graph.iter_depth_first_preorder_forward(|node| {
    let dkey = node.key;
    let parent = node.parent_keys.first().copied();

    let div = match mutation_counts {
      Some(counts) => match parent {
        Some((pkey, ekey)) => div_map.get(&pkey).copied().unwrap_or(0.0) + *counts.get(&ekey).unwrap_or(&0) as f64,
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

    if let Some((pkey, ekey)) = parent {
      let ir_parent = key_map[&pkey];
      let branch_length = node
        .parents
        .first()
        .and_then(|(_, edge)| edge.read_arc().branch_length());
      let mutations = match partition {
        Some(p) => subs_to_ir(&p.edge_subs(graph, ekey)?),
        None => vec![],
      };
      ir.add_edge(
        ir_parent,
        ir_key,
        TreeIrEdge {
          branch_length,
          mutations,
          ..TreeIrEdge::default()
        },
      )?;
    }

    Ok(())
  })?;

  ir.build()?;
  Ok(ir)
}

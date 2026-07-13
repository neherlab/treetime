use crate::partition::discrete_states::DiscreteStates;
use crate::partition::marginal_discrete::PartitionMarginalDiscrete;
use crate::partition::traits::PartitionBranchOps;
use crate::seq::mutation::Sub;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_io::graph::TreeIrGraph;
use treetime_io::tree_ir::mutation::{NUC_GENE, TreeIrSub};
use treetime_io::tree_ir::types::{TreeIrData, TreeIrEdge, TreeIrNode, TreeIrTrait};

pub fn sub_to_ir(sub: &Sub) -> TreeIrSub {
  TreeIrSub {
    gene: NUC_GENE.to_owned(),
    position: sub.pos() + 1,
    parent: sub.reff(),
    child: sub.qry(),
  }
}

pub fn subs_to_ir(subs: &[Sub]) -> Vec<TreeIrSub> {
  subs.iter().map(sub_to_ir).collect()
}

/// Build a TreeIR graph from a domain graph, extracting per-edge nucleotide
/// mutations from a partition. This is the common projection for commands that
/// have ancestral reconstruction (ancestral, optimize, timetree).
pub fn build_ir_with_mutations<N, E, D>(
  graph: &Graph<N, E, D>,
  partition: &dyn PartitionBranchOps,
  data: TreeIrData,
) -> Result<TreeIrGraph, Report>
where
  N: GraphNode + Named,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  let mut ir = TreeIrGraph::with_data(data);
  let mut key_map: BTreeMap<GraphNodeKey, GraphNodeKey> = BTreeMap::new();

  graph.iter_depth_first_preorder_forward(|node| {
    let dkey = node.key;
    let parent = node.parent_keys.first().copied();

    let ir_node = TreeIrNode {
      name: node.payload.name().map(|n| n.as_ref().to_owned()),
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
      let subs = partition.edge_subs(graph, ekey)?;
      ir.add_edge(
        ir_parent,
        ir_key,
        TreeIrEdge {
          branch_length,
          mutations: subs_to_ir(&subs),
          ..TreeIrEdge::default()
        },
      )?;
    }

    Ok(())
  })?;

  ir.build()?;
  Ok(ir)
}

/// Build a TreeIR graph from a domain graph without mutations. For commands
/// that lack ancestral reconstruction (clock, prune without alignment).
pub fn build_ir_topology_only<N, E, D>(
  graph: &Graph<N, E, D>,
  data: TreeIrData,
  node_mapper: impl Fn(GraphNodeKey, &N) -> TreeIrNode,
) -> Result<TreeIrGraph, Report>
where
  N: GraphNode + Named,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  let mut ir = TreeIrGraph::with_data(data);
  let mut key_map: BTreeMap<GraphNodeKey, GraphNodeKey> = BTreeMap::new();

  graph.iter_depth_first_preorder_forward(|node| {
    let dkey = node.key;
    let parent = node.parent_keys.first().copied();

    let ir_node = node_mapper(dkey, &*node.payload);
    let ir_key = ir.add_node(ir_node);
    key_map.insert(dkey, ir_key);

    if let Some((pkey, _ekey)) = parent {
      let ir_parent = key_map[&pkey];
      let branch_length = node
        .parents
        .first()
        .and_then(|(_, edge)| edge.read_arc().branch_length());
      ir.add_edge(
        ir_parent,
        ir_key,
        TreeIrEdge {
          branch_length,
          ..TreeIrEdge::default()
        },
      )?;
    }

    Ok(())
  })?;

  ir.build()?;
  Ok(ir)
}

/// Build a TreeIR graph from a mugration result with discrete trait
/// assignments and confidence.
pub fn build_ir_mugration<N, E, D>(
  graph: &Graph<N, E, D>,
  partition: &PartitionMarginalDiscrete,
  attribute: &str,
) -> Result<TreeIrGraph, Report>
where
  N: GraphNode + Named,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  let data = TreeIrData {
    trait_attrs: vec![attribute.to_owned()],
    ..TreeIrData::default()
  };

  let mut ir = TreeIrGraph::with_data(data);
  let mut key_map: BTreeMap<GraphNodeKey, GraphNodeKey> = BTreeMap::new();

  graph.iter_depth_first_preorder_forward(|node| {
    let dkey = node.key;
    let parent = node.parent_keys.first().copied();

    let mut traits = BTreeMap::new();
    if let Some(trait_value) = partition.get_reconstructed_trait(dkey) {
      let confidence = partition
        .get_confidence(dkey)
        .map(|profile| build_confidence_map(&partition.states, &profile))
        .unwrap_or_default();
      let entropy = partition.get_confidence(dkey).map(|profile| compute_entropy(&profile));
      traits.insert(
        attribute.to_owned(),
        TreeIrTrait {
          value: trait_value,
          confidence,
          entropy,
        },
      );
    }

    let ir_node = TreeIrNode {
      name: node.payload.name().map(|n| n.as_ref().to_owned()),
      traits,
      ..TreeIrNode::default()
    };
    let ir_key = ir.add_node(ir_node);
    key_map.insert(dkey, ir_key);

    if let Some((pkey, _ekey)) = parent {
      let ir_parent = key_map[&pkey];
      let branch_length = node
        .parents
        .first()
        .and_then(|(_, edge)| edge.read_arc().branch_length());
      ir.add_edge(
        ir_parent,
        ir_key,
        TreeIrEdge {
          branch_length,
          ..TreeIrEdge::default()
        },
      )?;
    }

    Ok(())
  })?;

  ir.build()?;
  Ok(ir)
}

fn build_confidence_map(states: &DiscreteStates, profile: &ndarray::Array1<f64>) -> BTreeMap<String, f64> {
  let mut pairs: Vec<(&str, f64)> = states.iter().zip(profile.iter()).map(|(s, &p)| (s, p)).collect();
  pairs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
  pairs
    .into_iter()
    .filter(|(_, p)| *p > 0.001)
    .map(|(s, p)| (s.to_owned(), p))
    .collect()
}

fn compute_entropy(profile: &ndarray::Array1<f64>) -> f64 {
  const TINY: f64 = 1e-12;
  -profile.iter().map(|&p| p * (p + TINY).ln()).sum::<f64>()
}

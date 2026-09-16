use crate::commands::prune::result::{PruneNodeOut, PruneOutputMaps};
use crate::commands::shared::tree_output::{
  NUC_TRACK, auspice_data, auspice_from_graph, cumulative_branch_length_from, generation_date, mat_from_graph,
  node_name_value, sequence_auspice_node, write_tree_outputs,
};
use crate::seq::mutation::Mutation;
use eyre::Report;
use maplit::btreemap;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::auspice_types::AuspiceTree;
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::CommentProviders;
use treetime_io::usher_mat::UsherTree;

pub fn write_prune_tree_outputs(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, PruneNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &PruneOutputMaps,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let updated = generation_date();
  let names: BTreeMap<GraphNodeKey, Option<String>> =
    nodes.iter().map(|(key, node)| (*key, node.name.clone())).collect();
  write_tree_outputs(
    graph,
    &names,
    branch_lengths,
    branch_lengths,
    outputs,
    providers,
    "prune",
    || prune_to_auspice(graph, nodes, branch_lengths, maps, &updated),
    || prune_to_mat(graph, &names, branch_lengths, maps),
  )
}

pub(crate) fn prune_to_auspice(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, PruneNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &PruneOutputMaps,
  updated: &str,
) -> Result<AuspiceTree, Report> {
  let root_sequences = prune_root_sequences(maps);
  let data = auspice_data(
    "TreeTime prune analysis",
    updated,
    vec![],
    vec![],
    None,
    None,
    Some(root_sequences),
    maps.root_sequence.is_some(),
  );
  auspice_from_graph(graph, data, |context| {
    let out = &nodes[&context.node_key];
    let name = node_name_value(context.node_key, out.name.as_deref());
    let div = cumulative_branch_length_from(graph, branch_lengths, context.node_key)?;
    sequence_auspice_node(
      &name,
      div,
      out.confidence,
      prune_mutations(maps, context.edge_key),
      None,
      None,
    )
  })
}

pub(crate) fn prune_to_mat(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  nwk_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &PruneOutputMaps,
) -> Result<UsherTree, Report> {
  let reference = prune_root_sequences(maps).remove(NUC_TRACK);
  mat_from_graph(
    graph,
    names,
    nwk_weights,
    reference.as_deref(),
    |_node_key, edge_key| Ok(prune_mutations(maps, Some(edge_key))),
  )
}

fn prune_root_sequences(maps: &PruneOutputMaps) -> BTreeMap<String, String> {
  maps
    .root_sequence
    .as_ref()
    .map(|sequence| btreemap! { NUC_TRACK.to_owned() => sequence.to_string() })
    .unwrap_or_default()
}

fn prune_mutations(maps: &PruneOutputMaps, edge_key: Option<GraphEdgeKey>) -> Vec<Mutation> {
  edge_key
    .and_then(|edge_key| maps.edge_mutations.get(&edge_key).cloned())
    .unwrap_or_default()
}

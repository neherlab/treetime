use crate::commands::clock::run::ClockNodeOut;
use crate::commands::shared::tree_output::{
  COLORING_BAD_BRANCH, COLORING_NUM_DATE, auspice_data, auspice_from_graph, auspice_node, coloring, finite_number,
  generation_date, mutation_free_mat, node_name_value, write_tree_outputs,
};
use eyre::Report;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::auspice_types::AuspiceTree;
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::CommentProviders;
use treetime_io::usher_mat::UsherTree;

pub fn write_clock_tree_outputs(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, ClockNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
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
    "clock",
    || clock_to_auspice(graph, nodes, &updated),
    || clock_to_mat(graph, &names, branch_lengths),
  )
}

pub(crate) fn clock_to_auspice(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, ClockNodeOut>,
  updated: &str,
) -> Result<AuspiceTree, Report> {
  let data = auspice_data(
    "TreeTime clock analysis",
    updated,
    vec![
      coloring(COLORING_NUM_DATE, "Date", "continuous"),
      coloring(COLORING_BAD_BRANCH, "Excluded", "categorical"),
    ],
    vec![COLORING_BAD_BRANCH.to_owned()],
    Some(COLORING_BAD_BRANCH.to_owned()),
    None,
    None,
    false,
  );
  auspice_from_graph(graph, data, |context| {
    let out = &nodes[&context.node_key];
    let name = node_name_value(context.node_key, out.name.as_deref());
    Ok(auspice_node(
      name.clone(),
      finite_number(Some(out.div), 6, "clock", &name, "div")?,
      finite_number(out.time, 3, "clock", &name, "date")?,
      None,
      Some(out.bad_branch || out.is_outlier),
      BTreeMap::new(),
      BTreeMap::new(),
      None,
    ))
  })
}

pub(crate) fn clock_to_mat(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  nwk_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<UsherTree, Report> {
  mutation_free_mat(graph, names, nwk_weights)
}

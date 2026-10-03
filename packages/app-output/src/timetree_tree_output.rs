use crate::timetree_result::{TimetreeEdgeOut, TimetreeNodeOut, TimetreeOutputMaps};
use crate::tree_output::{
  COLORING_BAD_BRANCH, COLORING_NUM_DATE, MatOutput, NUC_TRACK, auspice_data, auspice_from_graph, auspice_node,
  coloring, ensure_finite, finite_number, format_number, generation_date, group_mutations, mat_from_graph,
  node_name_value, write_tree_outputs,
};
use eyre::Report;
use maplit::btreemap;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime::clock::divergence::root_to_node_divergences;
use treetime::progress::LogSink;
use treetime::seq::mutation::Mutation;
use treetime::timetree::confidence::NodeConfidenceInterval;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::auspice_types::AuspiceTree;
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::CommentProviders;

pub fn write_timetree_tree_outputs(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, TimetreeNodeOut>,
  edges: &BTreeMap<GraphEdgeKey, TimetreeEdgeOut>,
  maps: &TimetreeOutputMaps,
  confidence_intervals: Option<&[NodeConfidenceInterval]>,
  mutation_counts: Option<&BTreeMap<GraphEdgeKey, usize>>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let updated = generation_date();
  let names: BTreeMap<GraphNodeKey, Option<String>> =
    nodes.iter().map(|(key, node)| (*key, node.name.clone())).collect();
  let nwk_weights: BTreeMap<GraphEdgeKey, Option<f64>> = edges
    .iter()
    .map(|(key, edge)| (*key, edge.date_branch_length))
    .collect();
  let graphviz_weights: BTreeMap<GraphEdgeKey, Option<f64>> =
    edges.iter().map(|(key, edge)| (*key, edge.branch_length)).collect();
  write_tree_outputs(
    graph,
    &names,
    &nwk_weights,
    &graphviz_weights,
    outputs,
    providers,
    "timetree",
    || timetree_to_auspice(graph, nodes, maps, confidence_intervals, mutation_counts, &updated),
    || timetree_to_mat(graph, &names, &nwk_weights, maps),
    log,
  )
}

pub(crate) fn timetree_to_auspice(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, TimetreeNodeOut>,
  maps: &TimetreeOutputMaps,
  confidence_intervals: Option<&[NodeConfidenceInterval]>,
  mutation_counts: Option<&BTreeMap<GraphEdgeKey, usize>>,
  updated: &str,
) -> Result<AuspiceTree, Report> {
  let root_sequences = timetree_root_sequences(maps);
  let data = auspice_data(
    "TreeTime timetree analysis",
    updated,
    vec![
      coloring(COLORING_NUM_DATE, "Date", "continuous"),
      coloring(COLORING_BAD_BRANCH, "Excluded", "categorical"),
    ],
    vec![COLORING_BAD_BRANCH.to_owned()],
    Some(COLORING_BAD_BRANCH.to_owned()),
    BTreeMap::new(),
    Some(root_sequences),
    maps.root_sequence.is_some(),
  )?;
  let mutation_divergences = mutation_counts
    .map(|counts| mutation_divergences(graph, counts))
    .transpose()?;
  auspice_from_graph(graph, data, |context| {
    let out = &nodes[&context.node_key];
    let name = node_name_value(context.node_key, out.name.as_deref());
    let div = mutation_divergences
      .as_ref()
      .map_or(out.div, |divergences| divergences[&context.node_key]);
    let confidence = timetree_date_confidence(graph, context.node_key, &name, confidence_intervals)?;
    Ok(auspice_node(
      name.clone(),
      finite_number(Some(div), 6, "timetree", &name, "div")?,
      finite_number(out.time, 3, "timetree", &name, "date")?,
      confidence,
      Some(out.bad_branch || out.is_outlier),
      BTreeMap::new(),
      group_mutations(timetree_mutations(maps, context.edge_key))?,
      None,
    ))
  })
}

pub(crate) fn timetree_to_mat(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  nwk_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &TimetreeOutputMaps,
) -> Result<MatOutput, Report> {
  let reference = timetree_root_sequences(maps).remove(NUC_TRACK);
  mat_from_graph(
    graph,
    names,
    nwk_weights,
    reference.as_deref(),
    |_node_key, edge_key| Ok(timetree_mutations(maps, Some(edge_key))),
  )
}

fn timetree_root_sequences(maps: &TimetreeOutputMaps) -> BTreeMap<String, String> {
  maps
    .root_sequence
    .as_ref()
    .map(|sequence| btreemap! { NUC_TRACK.to_owned() => sequence.to_string() })
    .unwrap_or_default()
}

fn timetree_mutations(maps: &TimetreeOutputMaps, edge_key: Option<GraphEdgeKey>) -> Vec<Mutation> {
  edge_key
    .and_then(|edge_key| maps.edge_mutations.get(&edge_key).cloned())
    .unwrap_or_default()
}

#[expect(
  clippy::as_conversions,
  reason = "a mutation count is far below 2^53, so the conversion to f64 is exact"
)]
fn mutation_divergences(
  graph: &Graph,
  mutation_counts: &BTreeMap<GraphEdgeKey, usize>,
) -> Result<BTreeMap<GraphNodeKey, f64>, Report> {
  root_to_node_divergences(graph, |edge_key| {
    mutation_counts.get(&edge_key).copied().unwrap_or_default() as f64
  })
}

fn timetree_date_confidence(
  graph: &Graph,
  node_key: GraphNodeKey,
  node_name: &str,
  confidence_intervals: Option<&[NodeConfidenceInterval]>,
) -> Result<Option<[f64; 2]>, Report> {
  let confidence = confidence_intervals
    .and_then(|intervals| intervals.iter().find(|interval| interval.key == node_key))
    .map(|interval| [interval.lower, interval.upper]);
  if let Some([lower, upper]) = confidence {
    ensure_finite(lower, "timetree", node_name, "date confidence lower bound")?;
    ensure_finite(upper, "timetree", node_name, "date confidence upper bound")?;
    Ok(Some([format_number(lower, 3), format_number(upper, 3)]))
  } else {
    Ok(None)
  }
}

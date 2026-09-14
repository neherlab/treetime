use crate::commands::shared::tree_output::{
  APPLIES_BRANCH, APPLIES_NODE, COLORING_BAD_BRANCH, COLORING_NUM_DATE, DT_BOOLEAN, DT_DOUBLE, NUC_TRACK,
  REF_BAD_BRANCH, REF_DATE_INFERRED, REF_DIV, REF_GAMMA, auspice_data, auspice_from_graph, auspice_node, coloring,
  empty_phyloxml_clade, ensure_finite, ensure_optional_finite, finite_number, format_number, generation_date,
  group_mutations, input_branch_confidence, mat_from_graph, mutation_property, node_name_value, phyloxml_from_graph,
  phyloxml_sequences, property, write_tree_outputs,
};
use crate::commands::timetree::result::{TimetreeEdgeOut, TimetreeNodeOut, TimetreeOutputMaps};
use crate::seq::mutation::Mutation;
use crate::timetree::confidence::NodeConfidenceInterval;
use eyre::Report;
use maplit::btreemap;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::auspice_types::AuspiceTree;
use treetime_io::dates_csv::DatesMap;
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::CommentProviders;
use treetime_io::phyloxml::{Phyloxml, PhyloxmlDate};
use treetime_io::usher_mat::UsherTree;

pub fn write_timetree_tree_outputs(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, TimetreeNodeOut>,
  edges: &BTreeMap<GraphEdgeKey, TimetreeEdgeOut>,
  maps: &TimetreeOutputMaps,
  confidence_intervals: Option<&[NodeConfidenceInterval]>,
  mutation_counts: Option<&BTreeMap<GraphEdgeKey, usize>>,
  dates: Option<&DatesMap>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let updated = generation_date();
  let names: BTreeMap<GraphNodeKey, Option<String>> =
    nodes.iter().map(|(key, node)| (*key, node.name.clone())).collect();
  // R2: the Newick/Nexus weight and the embedded MAT Newick weight are the branch time length, while
  // the Graphviz weight stays the substitution branch length. These diverge for timetree, so the two
  // writer paths take distinct edge-weight maps.
  let nwk_weights: BTreeMap<GraphEdgeKey, Option<f64>> =
    edges.iter().map(|(key, edge)| (*key, edge.time_length)).collect();
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
    || timetree_to_phyloxml(graph, nodes, edges, maps, confidence_intervals, mutation_counts, dates),
    || timetree_to_mat(graph, &names, &nwk_weights, maps),
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
    None,
    Some(root_sequences),
    maps.root_sequence.is_some(),
  );
  auspice_from_graph(graph, data, |context| {
    let out = &nodes[&context.node_key];
    let name = node_name_value(context.node_key, out.name.as_deref());
    let div = timetree_divergence(graph, context.node_key, out.div, mutation_counts)?;
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

pub(crate) fn timetree_to_phyloxml(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, TimetreeNodeOut>,
  edges: &BTreeMap<GraphEdgeKey, TimetreeEdgeOut>,
  maps: &TimetreeOutputMaps,
  confidence_intervals: Option<&[NodeConfidenceInterval]>,
  mutation_counts: Option<&BTreeMap<GraphEdgeKey, usize>>,
  dates: Option<&DatesMap>,
) -> Result<Phyloxml, Report> {
  phyloxml_from_graph(graph, "TreeTime timetree analysis", |context| {
    let out = &nodes[&context.node_key];
    let gamma = context.edge_key.map_or(1.0, |edge_key| edges[&edge_key].gamma);
    let name = node_name_value(context.node_key, out.name.as_deref());
    let divergence = timetree_divergence(graph, context.node_key, out.div, mutation_counts)?;
    ensure_finite(divergence, "timetree", &name, "divergence")?;
    ensure_optional_finite(out.time, "timetree", &name, "date")?;
    ensure_finite(gamma, "timetree", &name, "gamma")?;
    let mut properties = vec![
      property(REF_DIV, DT_DOUBLE, APPLIES_NODE, &divergence.to_string()),
      property(
        REF_BAD_BRANCH,
        DT_BOOLEAN,
        APPLIES_NODE,
        if out.bad_branch || out.is_outlier {
          "true"
        } else {
          "false"
        },
      ),
    ];
    if timetree_date_is_inferred(graph, context.node_key, out.name.as_deref(), out.time, dates) == Some(true) {
      properties.push(property(REF_DATE_INFERRED, DT_BOOLEAN, APPLIES_NODE, "true"));
    }
    if context.edge_key.is_some() {
      properties.push(property(REF_GAMMA, DT_DOUBLE, APPLIES_BRANCH, &gamma.to_string()));
    }
    for mutation in timetree_mutations(maps, context.edge_key) {
      properties.push(mutation_property(&mutation)?);
    }
    let date = out.time.map(|value| {
      let confidence =
        confidence_intervals.and_then(|intervals| intervals.iter().find(|interval| interval.key == context.node_key));
      PhyloxmlDate {
        desc: None,
        value: Some(value),
        minimum: confidence.map(|interval| interval.lower),
        maximum: confidence.map(|interval| interval.upper),
        unit: Some("year".to_owned()),
      }
    });
    let mut clade = empty_phyloxml_clade(
      out.name.clone(),
      context.edge_key.and_then(|edge_key| edges[&edge_key].branch_length),
    );
    clade.confidence = input_branch_confidence(out.confidence, "timetree", &name)?;
    clade.date = date;
    clade.property = properties;
    clade.sequence = phyloxml_sequences(&timetree_node_sequences(maps, context.node_key));
    Ok(clade)
  })
}

pub(crate) fn timetree_to_mat(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  nwk_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &TimetreeOutputMaps,
) -> Result<UsherTree, Report> {
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

fn timetree_node_sequences(maps: &TimetreeOutputMaps, node_key: GraphNodeKey) -> BTreeMap<String, String> {
  maps
    .node_sequences
    .get(&node_key)
    .map(|sequence| btreemap! { NUC_TRACK.to_owned() => sequence.to_string() })
    .unwrap_or_default()
}

fn timetree_mutations(maps: &TimetreeOutputMaps, edge_key: Option<GraphEdgeKey>) -> Vec<Mutation> {
  edge_key
    .and_then(|edge_key| maps.edge_mutations.get(&edge_key).cloned())
    .unwrap_or_default()
}

fn timetree_divergence(
  graph: &Graph,
  node_key: GraphNodeKey,
  div: f64,
  mutation_counts: Option<&BTreeMap<GraphEdgeKey, usize>>,
) -> Result<f64, Report> {
  mutation_counts.map_or(Ok(div), |counts| {
    let mut key = node_key;
    let mut count = 0;
    while let Some((parent, edge)) = graph.node_parent(key)? {
      count += counts.get(&edge).copied().unwrap_or_default();
      key = parent;
    }
    Ok(count as f64)
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

fn timetree_date_is_inferred(
  graph: &Graph,
  node_key: GraphNodeKey,
  name: Option<&str>,
  time: Option<f64>,
  dates: Option<&DatesMap>,
) -> Option<bool> {
  let dates = dates?;
  Some(
    name.and_then(|name| dates.get(name)).and_then(Option::as_ref).is_none()
      && time.is_some()
      && graph.get_node(node_key).is_some(),
  )
}

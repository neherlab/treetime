use crate::commands::clock::run::ClockNodeOut;
use crate::commands::shared::tree_output::{
  APPLIES_NODE, COLORING_BAD_BRANCH, COLORING_NUM_DATE, DT_BOOLEAN, DT_DOUBLE, REF_BAD_BRANCH, REF_DIV, auspice_data,
  auspice_from_graph, auspice_node, coloring, ensure_finite, ensure_optional_finite, finite_number, generation_date,
  mutation_free_mat, node_name_value, phyloxml_date, phyloxml_from_graph, property, write_tree_outputs,
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
use treetime_io::phyloxml::{Phyloxml, PhyloxmlClade};
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
    || clock_to_phyloxml(graph, nodes, branch_lengths),
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

pub(crate) fn clock_to_phyloxml(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, ClockNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<Phyloxml, Report> {
  phyloxml_from_graph(graph, "TreeTime clock analysis", |context| {
    let out = &nodes[&context.node_key];
    let name = node_name_value(context.node_key, out.name.as_deref());
    ensure_optional_finite(out.time, "clock", &name, "date")?;
    ensure_finite(out.div, "clock", &name, "divergence")?;
    let property = vec![
      property(REF_DIV, DT_DOUBLE, APPLIES_NODE, &out.div.to_string()),
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
    Ok(PhyloxmlClade {
      name: out.name.clone(),
      branch_length_elem: context.edge_key.and_then(|edge_key| branch_lengths[&edge_key]),
      branch_length_attr: None,
      confidence: vec![],
      width: None,
      color: None,
      node_id: None,
      taxonomy: vec![],
      sequence: vec![],
      events: None,
      binary_characters: None,
      distribution: vec![],
      date: out.time.map(phyloxml_date),
      reference: vec![],
      property,
      clade: vec![],
      other: BTreeMap::new(),
    })
  })
}

pub(crate) fn clock_to_mat(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  nwk_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<UsherTree, Report> {
  mutation_free_mat(graph, names, nwk_weights)
}

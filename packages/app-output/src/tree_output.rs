use crate::annotated_graph::{AnnotatedGraph, AnnotatedTreeView};
use crate::auspice::auspice_tree;
use crate::nwk_comments::nwk_node_comments;
use crate::output_plan::{CommandKind, ResolvedOutputs, TreeWriteKind};
use crate::usher_mat::mat_tree;
use chrono::Utc;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime::progress::LogSink;
use treetime::progress_warn;
use treetime_io::graphviz::graphviz_write_file;
use treetime_io::nex::nex_write_file;
use treetime_io::nwk::{NwkNodeComments, NwkStyle, NwkWriteOptions, nwk_write_file};
use treetime_io::usher_mat::usher_mat_pb_write_file;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub fn write_graph_outputs(
  graph: &AnnotatedGraph<'_>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
) -> Result<(), Report> {
  for (kind, path) in outputs {
    match kind {
      TreeWriteKind::GraphJson => json_write_file(path, graph.graph, JsonPretty(true))?,
      TreeWriteKind::Dot => graphviz_write_file(path, graph.graph, graph.names, graph.divergence_branch_lengths)?,
      TreeWriteKind::Nwk(_)
      | TreeWriteKind::Nexus(_)
      | TreeWriteKind::Auspice
      | TreeWriteKind::MatPb
      | TreeWriteKind::MatJson => {},
    }
  }
  Ok(())
}

pub fn tree_view_for_outputs<'a>(
  graph: &'a AnnotatedGraph<'a>,
  outputs: &ResolvedOutputs,
) -> Result<Option<AnnotatedTreeView<'a>>, Report> {
  let paths = outputs.tree_based_paths();
  if paths.is_empty() {
    return Ok(None);
  }
  AnnotatedTreeView::new(graph).map(Some).wrap_err_with(|| {
    format!(
      "These outputs need a tree and were not written: {}",
      paths.iter().map(|path| format!("'{}'", path.display())).join(", ")
    )
  })
}

pub fn write_tree_outputs(
  tree: &AnnotatedTreeView<'_>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  command: CommandKind,
  log: &dyn LogSink,
) -> Result<(), Report> {
  write_tree_formats(tree, outputs, command, log)
    .wrap_err_with(|| format!("When writing the tree outputs of {}", command.stem()))
}

fn write_tree_formats(
  tree: &AnnotatedTreeView<'_>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  command: CommandKind,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let graph = tree.graph();
  let comments = if outputs.keys().copied().any(is_annotated_newick) {
    nwk_node_comments(tree)?
  } else {
    NwkNodeComments::new()
  };
  for (kind, path) in outputs {
    match kind {
      TreeWriteKind::Nwk(style) => nwk_write_file(
        path,
        tree.tree(),
        graph.names,
        graph.tree_branch_lengths(),
        &nwk_options(*style),
        &comments,
      )?,
      TreeWriteKind::Nexus(style) => nex_write_file(
        path,
        tree.tree(),
        graph.names,
        graph.tree_branch_lengths(),
        &nwk_options(*style),
        &comments,
      )?,
      TreeWriteKind::Auspice => json_write_file(
        path,
        &auspice_tree(tree, command, &generation_date())?,
        JsonPretty(true),
      )?,
      TreeWriteKind::MatPb | TreeWriteKind::MatJson | TreeWriteKind::GraphJson | TreeWriteKind::Dot => {},
    }
  }
  write_mat_outputs(
    tree,
    outputs.get(&TreeWriteKind::MatPb),
    outputs.get(&TreeWriteKind::MatJson),
    log,
  )
}

fn is_annotated_newick(kind: TreeWriteKind) -> bool {
  match kind {
    TreeWriteKind::Nwk(style) | TreeWriteKind::Nexus(style) => style != NwkStyle::Plain,
    TreeWriteKind::Auspice
    | TreeWriteKind::MatPb
    | TreeWriteKind::MatJson
    | TreeWriteKind::GraphJson
    | TreeWriteKind::Dot => false,
  }
}

fn nwk_options(style: NwkStyle) -> NwkWriteOptions {
  NwkWriteOptions {
    style,
    ..NwkWriteOptions::default()
  }
}

fn write_mat_outputs(
  tree: &AnnotatedTreeView<'_>,
  pb_path: Option<&PathBuf>,
  json_path: Option<&PathBuf>,
  log: &dyn LogSink,
) -> Result<(), Report> {
  if pb_path.is_none() && json_path.is_none() {
    return Ok(());
  }
  let mat = mat_tree(tree)?;
  if let Some(warning) = mat.gaps.warning() {
    progress_warn!(log, "{warning}");
  }
  if let Some(path) = pb_path {
    usher_mat_pb_write_file(path, &mat.tree)?;
  }
  if let Some(path) = json_path {
    json_write_file(path, &mat.tree, JsonPretty(true))?;
  }
  Ok(())
}

fn generation_date() -> String {
  Utc::now().format("%Y-%m-%d").to_string()
}

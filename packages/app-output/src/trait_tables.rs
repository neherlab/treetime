use crate::annotated_graph::{AnnotatedTreeView, TreeTraits};
use crate::output_plan::OutputSelection;
use crate::table_output::table_create;
use eyre::Report;
use std::io::Write;
use std::iter::once;
use std::path::Path;
use treetime_graph::assign_node_names::node_name_or_key;
use treetime_io::csv::CsvWriter;
use treetime_utils::make_internal_report;

pub fn write_traits_csv(tree: &AnnotatedTreeView<'_>, path: &Path) -> Result<(), Report> {
  let mut csv = table_create(OutputSelection::TraitsCsv, path)?;
  write_traits_rows(tree, &mut csv)?;
  csv.finish()
}

pub fn write_trait_confidence_csv(tree: &AnnotatedTreeView<'_>, path: &Path) -> Result<(), Report> {
  let mut csv = table_create(OutputSelection::ConfidenceCsv, path)?;
  write_trait_confidence_rows(tree, &mut csv)?;
  csv.finish()
}

pub(crate) fn write_traits_rows<W: Write>(tree: &AnnotatedTreeView<'_>, csv: &mut CsvWriter<W>) -> Result<(), Report> {
  let traits = tree_traits(tree)?;
  let graph = tree.graph();
  csv.write_record(["node", traits.attribute])?;
  graph.graph.get_nodes().try_for_each(|node| {
    let key = node.key();
    let Some(value) = &traits.values[&key] else {
      return Ok(());
    };
    let name = node_name_or_key(key, graph.names[&key].as_deref());
    csv.write_record([name.as_str(), value.as_str()])
  })
}

pub(crate) fn write_trait_confidence_rows<W: Write>(
  tree: &AnnotatedTreeView<'_>,
  csv: &mut CsvWriter<W>,
) -> Result<(), Report> {
  let traits = tree_traits(tree)?;
  let graph = tree.graph();
  csv.write_record(once("node").chain(traits.states.iter()))?;
  graph.graph.get_nodes().try_for_each(|node| {
    let key = node.key();
    let Some(profile) = &traits.profiles[&key] else {
      return Ok(());
    };
    let name = node_name_or_key(key, graph.names[&key].as_deref());
    csv.write_record(once(name).chain(profile.iter().map(|probability| format!("{probability:.6}"))))
  })
}

fn tree_traits<'a>(tree: &'a AnnotatedTreeView<'_>) -> Result<&'a TreeTraits<'a>, Report> {
  tree
    .graph()
    .traits
    .as_ref()
    .ok_or_else(|| make_internal_report!("Trait tables require reconstructed traits"))
}

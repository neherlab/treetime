use crate::annotated_graph::{AnnotatedTreeView, TreeTraits};
use crate::output_plan::OutputSelection;
use crate::table_output::table_create;
use deser::Serialize;
use eyre::Report;
use itertools::Itertools;
use std::borrow::Cow;
use std::io::Write;
use std::iter::once;
use std::path::Path;
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
  let ambiguity_key = traits.ambiguity_key(OutputSelection::TraitsCsv);
  csv.write_record(["node", traits.attribute].into_iter().chain(ambiguity_key.as_deref()))?;
  graph.graph.get_nodes().try_for_each(|node| {
    let key = node.key();
    let ambiguous = ambiguity_key
      .is_some()
      .then(|| traits.ambiguous_states(key, OutputSelection::TraitsCsv))
      .transpose()?
      .map(Option::unwrap_or_default);
    csv.write_record(
      [Cow::Borrowed(&graph.names[key]), Cow::Borrowed(traits.state_name(key)?)]
        .into_iter()
        .chain(ambiguous.map(Cow::Owned)),
    )
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
    let profile = traits.posteriors[&key].profile();
    let row = once(ConfidenceField::Name(&graph.names[key]))
      .chain(
        profile
          .iter()
          .map(|&probability| ConfidenceField::Probability(probability)),
      )
      .collect_vec();
    csv.write_row(&row)
  })
}

#[derive(Serialize)]
#[deser(untagged)]
enum ConfidenceField<'a> {
  Name(&'a str),
  Probability(f64),
}

fn tree_traits<'a>(tree: &'a AnnotatedTreeView<'_>) -> Result<&'a TreeTraits<'a>, Report> {
  tree
    .graph()
    .traits
    .as_ref()
    .ok_or_else(|| make_internal_report!("Trait tables require reconstructed traits"))
}

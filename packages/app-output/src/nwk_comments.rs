use crate::annotated_graph::AnnotatedTreeView;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime::seq::mutation::{Mutation, MutationEvent, mutation_event_strings};
use treetime_io::nwk::NwkNodeComments;

pub(crate) fn nwk_node_comments(tree: &AnnotatedTreeView<'_>) -> Result<NwkNodeComments, Report> {
  let graph = tree.graph();
  let mut comments = NwkNodeComments::new();
  for &key in tree.tree().preorder() {
    let mut node_comments = BTreeMap::new();
    if let Some(sequences) = &graph.sequences
      && let Some((_, edge_key)) = tree.tree().parent(key)
      && let Some(mutations) = mutation_comment(&sequences.edge_mutations[&edge_key])?
    {
      node_comments.insert("mutations".to_owned(), mutations);
    }
    if let Some(dates) = &graph.dates
      && let Some(date) = dates.num_date[&key]
    {
      node_comments.insert("date".to_owned(), format!("{date:.2}"));
    }
    if let Some(traits) = &graph.traits
      && let Some(value) = &traits.values[&key]
    {
      node_comments.insert(traits.attribute.to_owned(), value.clone());
    }
    comments.insert(key, node_comments);
  }
  Ok(comments)
}

fn mutation_comment(mutations: &[Mutation]) -> Result<Option<String>, Report> {
  if mutations.is_empty() {
    return Ok(None);
  }
  let strings: Vec<Vec<String>> = mutations
    .iter()
    .sorted_by_key(|mutation| match &mutation.event {
      MutationEvent::Substitution(substitution) => substitution.pos(),
      MutationEvent::Insertion(segment) | MutationEvent::Deletion(segment) => segment.range.0,
    })
    .map(|mutation| mutation_event_strings(&mutation.event))
    .try_collect()?;
  Ok(Some(strings.into_iter().flatten().join(",")))
}

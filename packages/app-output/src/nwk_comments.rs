use crate::annotated_graph::AnnotatedTreeView;
use eyre::Report;
use itertools::Itertools;
use treetime::seq::mutation::{Mutation, MutationEvent, mutation_event_strings};
use treetime_io::nwk::{NewickValue, NwkNodeComments};

pub(crate) fn nwk_node_comments(tree: &AnnotatedTreeView<'_>) -> Result<NwkNodeComments, Report> {
  let graph = tree.graph();
  let mut comments = NwkNodeComments::new();
  for &key in tree.tree().preorder() {
    let mut node_comments = vec![];
    if let Some(sequences) = &graph.sequences
      && let Some((_, edge_key)) = tree.tree().parent(key)
      && let Some(mutations) = mutation_comment(&sequences.edge_mutations[&edge_key])?
    {
      node_comments.push(("mutations".to_owned(), NewickValue::String(mutations)));
    }
    if let Some(dates) = &graph.dates
      && let Some(date) = dates.num_date[&key]
    {
      node_comments.push(("date".to_owned(), NewickValue::NumberText(format!("{date:.2}"))));
    }
    if let Some(traits) = &graph.traits
      && let Some(value) = traits.values[&key].as_ref().filter(|value| !value.is_empty())
    {
      node_comments.push((traits.attribute.to_owned(), NewickValue::String(value.clone())));
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

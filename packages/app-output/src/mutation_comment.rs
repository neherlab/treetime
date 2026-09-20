use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use std::collections::BTreeMap;
use treetime::seq::mutation::{Mutation, MutationEvent, mutation_event_strings};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::nwk::NodeCommentProvider;

pub struct EdgeMutationCommentProvider<'a> {
  edge_mutations: &'a BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  graph: &'a Graph,
}

impl<'a> EdgeMutationCommentProvider<'a> {
  pub fn new(edge_mutations: &'a BTreeMap<GraphEdgeKey, Vec<Mutation>>, graph: &'a Graph) -> Self {
    Self { edge_mutations, graph }
  }
}

impl NodeCommentProvider for EdgeMutationCommentProvider<'_> {
  fn node_comments(&self, key: GraphNodeKey) -> Result<BTreeMap<String, String>, Report> {
    let Some((_parent_key, edge_key)) = self.graph.node_parent(key)? else {
      return Ok(BTreeMap::new());
    };
    let mut mutations = self.edge_mutations[&edge_key].clone();
    if mutations.is_empty() {
      return Ok(BTreeMap::new());
    }
    mutations.sort_by_key(|mutation| match &mutation.event {
      MutationEvent::Substitution(substitution) => substitution.pos(),
      MutationEvent::Insertion(segment) | MutationEvent::Deletion(segment) => segment.range.0,
    });
    let mutations = mutations
      .iter()
      .map(|mutation| mutation_event_strings(&mutation.event))
      .collect::<Result<Vec<_>, _>>()?
      .into_iter()
      .flatten()
      .join(",");
    Ok(btreemap! {
      "mutations".to_owned() => mutations,
    })
  }
}

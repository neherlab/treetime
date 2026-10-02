#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{c, edge_mutation_map};
  use app_output::EdgeMutationCommentProvider;
  use eyre::Report;
  use helpers::leaf_key;

  use pretty_assertions::assert_eq;

  use treetime::seq::indel::InDel;
  use treetime::seq::mutation::Sub;

  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::{NodeCommentProvider, nwk_read_str};
  use treetime_primitives::Seq;

  #[test]
  fn test_mutation_comment_provider_formats_1_based_substitutions_and_indels() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1)root;")?;
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let edge_subs = &[(
      0,
      vec![
        Sub::new(c(b'A'), 0_usize, c(b'T'))?,
        Sub::new(c(b'G'), 5_usize, c(b'C'))?,
      ],
    )];
    let edge_indels = &[(0, vec![InDel::del((1, 3), Seq::try_from_str("CG")?)?])];
    let edge_mutations = edge_mutation_map(&graph, edge_subs, edge_indels)?;
    let provider = EdgeMutationCommentProvider::new(&edge_mutations, &graph);
    let comments = provider.node_comments(leaf_key(&graph))?;
    assert_eq!(comments.get("mutations").map(String::as_str), Some("A1T,C2-,G3-,G6C"));
    Ok(())
  }

  #[test]
  fn test_mutation_comment_provider_root_has_no_comments() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1)root;")?;
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let edge_subs = &[(0, vec![Sub::new(c(b'A'), 0_usize, c(b'T'))?])];
    let edge_mutations = edge_mutation_map(&graph, edge_subs, &[])?;
    let provider = EdgeMutationCommentProvider::new(&edge_mutations, &graph);
    let root_key = graph.get_roots().collect::<Vec<_>>()[0].key();
    let comments = provider.node_comments(root_key)?;
    assert!(comments.is_empty());
    Ok(())
  }

  #[test]
  fn test_mutation_comment_provider_no_mutations_returns_empty() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1)root;")?;
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let edge_subs = &[(0, vec![])];
    let edge_mutations = edge_mutation_map(&graph, edge_subs, &[])?;
    let provider = EdgeMutationCommentProvider::new(&edge_mutations, &graph);
    let comments = provider.node_comments(leaf_key(&graph))?;
    assert!(comments.is_empty());
    Ok(())
  }

  #[test]
  fn test_mutation_comment_provider_sorts_by_position() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1)root;")?;
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let edge_subs = &[(
      0,
      vec![
        Sub::new(c(b'C'), 50_usize, c(b'G'))?,
        Sub::new(c(b'A'), 10_usize, c(b'T'))?,
        Sub::new(c(b'G'), 30_usize, c(b'C'))?,
      ],
    )];
    let edge_mutations = edge_mutation_map(&graph, edge_subs, &[])?;
    let provider = EdgeMutationCommentProvider::new(&edge_mutations, &graph);
    let comments = provider.node_comments(leaf_key(&graph))?;
    assert_eq!(comments.get("mutations").map(String::as_str), Some("A11T,G31C,C51G"));
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn leaf_key(graph: &Graph) -> GraphNodeKey {
      graph.get_leaves().collect::<Vec<_>>()[0].key()
    }
  }
}

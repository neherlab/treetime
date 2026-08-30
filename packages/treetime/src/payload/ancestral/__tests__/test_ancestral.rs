#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::partition::storage::sparse::{SparseEdgePartition, SparseNodePartition};
  use crate::partition::traits::MutationCommentProvider;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::indel::InDel;
  use crate::seq::mutation::Sub;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::{NodeCommentProvider, nwk_read_str};
  use treetime_primitives::{AsciiChar, Seq};

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn make_test_partition(
    graph: &GraphAncestral,
    length: usize,
    edge_subs: &[(usize, Vec<Sub>)],
  ) -> Result<PartitionMarginalSparse, Report> {
    let alphabet = Alphabet::default();
    let mut ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
    for (_, subs) in edge_subs {
      for s in subs {
        if s.pos() < length {
          ref_seq[s.pos()] = s.reff();
        }
      }
    }

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: alphabet.clone(),
      length,
      root_sequence: ref_seq.clone(),
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    for node in graph.get_nodes() {
      let key = node.read_arc().key();
      let mut node_part = SparseNodePartition::empty(&alphabet);
      node_part.seq.sequence = ref_seq.clone();
      partition.nodes.insert(key, node_part);
    }

    let edges = graph.get_edges();
    for (idx, subs) in edge_subs {
      if let Some(edge) = edges.get(*idx) {
        let edge_key = edge.read_arc().key();
        let mut edge_part = SparseEdgePartition::with_fitch_subs(subs.clone());
        edge_part.set_ml_subs(subs.clone());
        partition.edges.insert(edge_key, edge_part);
      }
    }

    Ok(partition)
  }

  fn leaf_key(graph: &GraphAncestral) -> GraphNodeKey {
    graph.get_leaves()[0].read_arc().key()
  }

  #[test]
  fn test_mutation_comment_provider_formats_1_based_substitutions_and_indels() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1)root;")?;
    let mut partition = make_test_partition(
      &graph,
      100,
      &[(
        0,
        vec![
          Sub::new(c(b'A'), 0_usize, c(b'T'))?,
          Sub::new(c(b'G'), 5_usize, c(b'C'))?,
        ],
      )],
    )?;
    let edge_key = graph.get_edges()[0].read_arc().key();
    partition
      .edges
      .get_mut(&edge_key)
      .expect("fixture edge partition must exist")
      .indels = vec![InDel::del((1, 3), Seq::try_from_str("CG")?)?];
    let provider = MutationCommentProvider::new(&partition, &graph);
    let comments = provider.node_comments(leaf_key(&graph))?;
    assert_eq!(comments.get("mutations").map(String::as_str), Some("A1T,C2-,G3-,G6C"));
    Ok(())
  }

  #[test]
  fn test_mutation_comment_provider_root_has_no_comments() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1)root;")?;
    let partition = make_test_partition(&graph, 100, &[(0, vec![Sub::new(c(b'A'), 0_usize, c(b'T'))?])])?;
    let provider = MutationCommentProvider::new(&partition, &graph);
    let root_key = graph.get_roots()[0].read_arc().key();
    let comments = provider.node_comments(root_key)?;
    assert!(comments.is_empty());
    Ok(())
  }

  #[test]
  fn test_mutation_comment_provider_no_mutations_returns_empty() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1)root;")?;
    let partition = make_test_partition(&graph, 100, &[(0, vec![])])?;
    let provider = MutationCommentProvider::new(&partition, &graph);
    let comments = provider.node_comments(leaf_key(&graph))?;
    assert!(comments.is_empty());
    Ok(())
  }

  #[test]
  fn test_mutation_comment_provider_sorts_by_position() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1)root;")?;
    let partition = make_test_partition(
      &graph,
      100,
      &[(
        0,
        vec![
          Sub::new(c(b'C'), 50_usize, c(b'G'))?,
          Sub::new(c(b'A'), 10_usize, c(b'T'))?,
          Sub::new(c(b'G'), 30_usize, c(b'C'))?,
        ],
      )],
    )?;
    let provider = MutationCommentProvider::new(&partition, &graph);
    let comments = provider.node_comments(leaf_key(&graph))?;
    assert_eq!(comments.get("mutations").map(String::as_str), Some("A11T,G31C,C51G"));
    Ok(())
  }
}

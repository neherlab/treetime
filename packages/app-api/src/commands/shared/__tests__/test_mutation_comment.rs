#[cfg(test)]
mod tests {
  use treetime::alphabet::alphabet::Alphabet;
  use treetime::ancestral::pipeline::SparseReconstruction;
  use crate::commands::shared::mutation_comment::EdgeMutationCommentProvider;
  use treetime::gtr::get_gtr::{JC69Params, jc69};
  use treetime::partition::marginal::shared::update::MarginalEdges;
  use treetime::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use treetime::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
  use treetime::seq::indel::InDel;
  use treetime::seq::mutation::{Mutation, MutationTrack, Sub};
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::{NodeCommentProvider, nwk_read_str};
  use treetime_primitives::{AsciiChar, Seq};

  #[test]
  fn test_mutation_comment_provider_formats_1_based_substitutions_and_indels() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1)root;")?;
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
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
    let edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    partition
      .partition
      .obs_edges
      .get_mut(&edge_key)
      .expect("fixture edge partition must exist")
      .indels = vec![InDel::del((1, 3), Seq::try_from_str("CG")?)?];
    let edge_mutations = edge_mutation_map(&graph, &partition)?;
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
    let partition = make_test_partition(&graph, 100, &[(0, vec![Sub::new(c(b'A'), 0_usize, c(b'T'))?])])?;
    let edge_mutations = edge_mutation_map(&graph, &partition)?;
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
    let partition = make_test_partition(&graph, 100, &[(0, vec![])])?;
    let edge_mutations = edge_mutation_map(&graph, &partition)?;
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
    let edge_mutations = edge_mutation_map(&graph, &partition)?;
    let provider = EdgeMutationCommentProvider::new(&edge_mutations, &graph);
    let comments = provider.node_comments(leaf_key(&graph))?;
    assert_eq!(comments.get("mutations").map(String::as_str), Some("A11T,G31C,C51G"));
    Ok(())
  }

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  /// Gather the per-edge nucleotide mutation map the comment provider consumes off a completed
  /// sparse reconstruction, mirroring how the tree writers gather it in production.
  fn edge_mutation_map(
    graph: &Graph,
    partition: &SparseReconstruction,
  ) -> Result<BTreeMap<GraphEdgeKey, Vec<Mutation>>, Report> {
    graph
      .get_edges()
      .map(|edge| {
        let key = edge.key();
        Ok((key, partition.edge_mutations(key, &MutationTrack::Nucleotide)?))
      })
      .collect()
  }

  fn make_test_partition(
    graph: &Graph,
    length: usize,
    edge_subs: &[(usize, Vec<Sub>)],
  ) -> Result<SparseReconstruction, Report> {
    let alphabet = Alphabet::default();
    let mut ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
    for (_, subs) in edge_subs {
      for s in subs {
        if s.pos() < length {
          ref_seq[s.pos()] = s.reff();
        }
      }
    }

    let mut obs_nodes = btreemap! {};
    let mut node_states = btreemap! {};
    for node in graph.get_nodes() {
      let key = node.key();
      obs_nodes.insert(key, SparseNodeObs::new(&ref_seq, &alphabet));
      node_states.insert(key, SparseNodeState::leaf(&ref_seq));
    }

    let mut obs_edges = btreemap! {};
    // The MAP substitutions the comment provider reports come from the estimates map; the fixture seeds
    // it directly (comment tests do not run a marginal pass).
    let mut estimates = btreemap! {};
    let edges = graph.get_edges().collect::<Vec<_>>();
    for (idx, subs) in edge_subs {
      if let Some(edge) = edges.get(*idx) {
        let edge_key = edge.key();
        obs_edges.insert(edge_key, SparseEdgeObs::with_fitch_subs(subs.clone()));
        estimates.insert(edge_key, subs.clone());
      }
    }

    let partition = PartitionMarginalSparse {
      index: 0,
      alphabet,
      length,
      root_sequence: ref_seq,
      obs_nodes,
      obs_edges,
    };

    Ok(SparseReconstruction {
      partition,
      gtr: jc69(JC69Params::default())?,
      node_states,
      edges: MarginalEdges {
        estimates,
        ..MarginalEdges::default()
      },
    })
  }

  fn leaf_key(graph: &Graph) -> GraphNodeKey {
    graph.get_leaves().collect::<Vec<_>>()[0].key()
  }
}

#[cfg(test)]
mod tests {
  use app_output::EdgeMutationCommentProvider;
  use app_output::DateCommentProvider;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime::alphabet::alphabet::Alphabet;
  use treetime::ancestral::pipeline::SparseReconstruction;
  use treetime::gtr::get_gtr::{JC69Params, jc69};
  use treetime::partition::marginal::shared::update::MarginalEdges;
  use treetime::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use treetime::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
  use treetime::seq::mutation::{Mutation, MutationTrack, Sub};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nex::{NexWriteOptions, nex_write_str_with};
  use treetime_io::nwk::{CommentProviders, NodeCommentProvider, NwkStyle, nwk_read_str};
  use treetime_primitives::AsciiChar;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  /// Gather the per-edge nucleotide mutation map the comment provider consumes off a completed
  /// sparse reconstruction, mirroring how the timetree tree writers gather it in production.
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
    let mut ref_seq: treetime_primitives::Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
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
    // MAP substitutions the comment provider reports come from the estimates map; the fixture seeds it
    // directly since these output tests run no marginal pass.
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

  #[test]
  fn test_timetree_mutation_provider_produces_comments() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let partition = make_test_partition(
      &graph,
      100,
      &[(
        0,
        vec![
          Sub::new(c(b'A'), 54_usize, c(b'G'))?,
          Sub::new(c(b'T'), 92_usize, c(b'C'))?,
        ],
      )],
    )?;
    let edge_mutations = edge_mutation_map(&graph, &partition)?;
    let provider = EdgeMutationCommentProvider::new(&edge_mutations, &graph);
    let leaf_key = graph.get_leaves().collect::<Vec<_>>()[0].key();
    let comments = provider.node_comments(leaf_key)?;
    assert_eq!(comments.get("mutations").map(String::as_str), Some("A55G,T93C"));
    Ok(())
  }

  #[test]
  fn test_timetree_nexus_output_includes_mutations_and_date() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let partition = make_test_partition(
      &graph,
      100,
      &[(
        0,
        vec![
          Sub::new(c(b'A'), 54_usize, c(b'G'))?,
          Sub::new(c(b'T'), 92_usize, c(b'C'))?,
        ],
      )],
    )?;

    let date_times: BTreeMap<GraphNodeKey, f64> = graph.get_leaves().map(|leaf| (leaf.key(), 2003.84)).collect();

    let edge_mutations = edge_mutation_map(&graph, &partition)?;
    let provider = EdgeMutationCommentProvider::new(&edge_mutations, &graph);
    let date_provider = DateCommentProvider::new(&date_times);
    let providers = CommentProviders::new().with(&provider).with(&date_provider);
    let options = NexWriteOptions {
      style: NwkStyle::Beast,
      ..NexWriteOptions::default()
    };
    let time_lengths: BTreeMap<GraphEdgeKey, Option<f64>> = graph.get_edges().map(|edge| (edge.key(), None)).collect();
    let nexus = nex_write_str_with(&graph, &names, &time_lengths, &options, &providers)?;
    let expected = concat!(
      indoc! {r#"
        #NEXUS
        Begin Taxa;
          Dimensions NTax=1;
          TaxLabels A;
        End;
        Begin Trees;
          Tree tree1=(A[&date=2003.84,mutations="A55G,T93C"])root;
        End;
      "#},
      "\n"
    );
    assert_eq!(nexus, expected);
    Ok(())
  }
}

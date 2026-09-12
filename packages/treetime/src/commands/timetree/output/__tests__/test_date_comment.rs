#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::pipeline::SparseReconstruction;
  use crate::commands::shared::mutation_comment::EdgeMutationCommentProvider;
  use crate::commands::timetree::output::date_comment::DateCommentProvider;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
  use crate::partition::traits::PartitionBranchOps;
  use crate::seq::mutation::{Mutation, MutationTrack, Sub};
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nex::{NexWriteOptions, nex_write_str_with};
  use treetime_io::nwk::{CommentProviders, NodeCommentProvider, NwkParse, NwkStyle, nwk_read_str};
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
    let readout = partition.readout();
    graph
      .get_edges()
      .iter()
      .map(|edge| {
        let key = edge.read_arc().key();
        Ok((key, readout.edge_mutations(graph, key, MutationTrack::Nucleotide)?))
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
      let key = node.read_arc().key();
      obs_nodes.insert(key, SparseNodeObs::new(&ref_seq, &alphabet));
      node_states.insert(key, SparseNodeState::leaf(&ref_seq));
    }

    let mut obs_edges = btreemap! {};
    // MAP substitutions the comment provider reports come from the estimates map; the fixture seeds it
    // directly since these output tests run no marginal pass.
    let mut estimates = btreemap! {};
    let edges = graph.get_edges();
    for (idx, subs) in edge_subs {
      if let Some(edge) = edges.get(*idx) {
        let edge_key = edge.read_arc().key();
        obs_edges.insert(edge_key, SparseEdgeObs::with_fitch_subs(subs.clone()));
        estimates.insert(edge_key, subs.clone());
      }
    }

    let partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length,
      root_sequence: ref_seq,
      obs_nodes,
      obs_edges,
    };

    Ok(SparseReconstruction {
      partition,
      node_states,
      backward: btreemap! {},
      forward: btreemap! {},
      estimates,
    })
  }

  #[test]
  fn test_timetree_mutation_provider_produces_comments() -> Result<(), Report> {
    let NwkParse { graph, names, .. } = nwk_read_str("(A:0.1)root;")?;
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
    let leaf_key = graph.get_leaves()[0].read_arc().key();
    let comments = provider.node_comments(leaf_key)?;
    assert_eq!(comments.get("mutations").map(String::as_str), Some("A55G,T93C"));
    Ok(())
  }

  #[test]
  fn test_timetree_nexus_output_includes_mutations_and_date() -> Result<(), Report> {
    let NwkParse { graph, names, .. } = nwk_read_str("(A:0.1)root;")?;
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

    let date_times: BTreeMap<GraphNodeKey, f64> = graph
      .get_leaves()
      .iter()
      .map(|leaf| (leaf.read_arc().key(), 2003.84))
      .collect();

    let edge_mutations = edge_mutation_map(&graph, &partition)?;
    let provider = EdgeMutationCommentProvider::new(&edge_mutations, &graph);
    let date_provider = DateCommentProvider::new(&date_times);
    let providers = CommentProviders::new().with(&provider).with(&date_provider);
    let options = NexWriteOptions {
      style: NwkStyle::Beast,
      ..NexWriteOptions::default()
    };
    let time_lengths: BTreeMap<GraphEdgeKey, Option<f64>> = graph
      .get_edges()
      .iter()
      .map(|edge| (edge.read_arc().key(), None))
      .collect();
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

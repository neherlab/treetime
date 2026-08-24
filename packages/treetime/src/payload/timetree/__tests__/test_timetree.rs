#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::partition::sparse::{SparseEdgePartition, SparseNodePartition};
  use crate::partition::timetree::GraphTimetree;
  use crate::partition::traits::MutationCommentProvider;
  use crate::seq::mutation::Sub;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_io::nex::{NexWriteOptions, nex_write_str_with};
  use treetime_io::nwk::{CommentProviders, NodeCommentProvider, NwkStyle, nwk_read_str};
  use treetime_primitives::AsciiChar;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn make_test_partition(
    graph: &GraphTimetree,
    length: usize,
    edge_subs: &[(usize, Vec<Sub>)],
  ) -> Result<PartitionMarginalSparse, Report> {
    let alphabet = Alphabet::default();
    let mut ref_seq: treetime_primitives::Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
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

  #[test]
  fn test_timetree_mutation_provider_produces_comments() -> Result<(), Report> {
    let graph: GraphTimetree = nwk_read_str("(A:0.1)root;")?;
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
    let provider = MutationCommentProvider::new(&partition, &graph);
    let leaf_key = graph.get_leaves()[0].read_arc().key();
    let comments = provider.node_comments(leaf_key)?;
    assert_eq!(comments.get("mutations").map(String::as_str), Some("A55G,T93C"));
    Ok(())
  }

  #[test]
  fn test_timetree_nexus_output_includes_mutations_and_date() -> Result<(), Report> {
    let graph: GraphTimetree = nwk_read_str("(A:0.1)root;")?;
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

    for leaf in graph.get_leaves() {
      leaf.read_arc().payload().write_arc().time = Some(2003.84);
    }

    let provider = MutationCommentProvider::new(&partition, &graph);
    let providers = CommentProviders::new().with(&provider);
    let options = NexWriteOptions {
      style: NwkStyle::Beast,
      ..NexWriteOptions::default()
    };
    let nexus = nex_write_str_with(&graph, &options, &providers)?;
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

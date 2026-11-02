#[cfg(test)]
mod tests {
  use crate::optimize::__tests__::test_topology_cleanup::tests::helpers::sub;
  use crate::optimize::topology::collapse::collapse_short_internal_edges;
  use crate::seq::indel::InDel;
  use crate::seq::mutation::Sub;
  use eyre::Report;
  use helpers::{ShortEdgeCase, gap_change};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  const ONE_MUTATION: f64 = 1.0 / 100.0;

  #[rustfmt::skip]
  #[rstest]
  #[case::short_without_change(   (0.05 * ONE_MUTATION, vec![],                   vec![]            ), "(C,D,A,B)root;"   )]
  #[case::zero_without_change(    (0.0,                 vec![],                   vec![]            ), "(C,D,A,B)root;"   )]
  #[case::short_with_substitution((0.05 * ONE_MUTATION, vec![sub(b'A', 0, b'T')], vec![]            ), "((A,B)P,C,D)root;")]
  #[case::zero_with_substitution( (0.0,                 vec![sub(b'A', 0, b'T')], vec![]            ), "((A,B)P,C,D)root;")]
  #[case::zero_with_gap_change(   (0.0,                 vec![],                   vec![gap_change()]), "((A,B)P,C,D)root;")]
  #[case::at_threshold(           (0.1 * ONE_MUTATION,  vec![],                   vec![]            ), "((A,B)P,C,D)root;")]
  #[case::long_without_change(    (0.2 * ONE_MUTATION,  vec![],                   vec![]            ), "((A,B)P,C,D)root;")]
  #[trace]
  fn test_topology_cleanup_collapses_internal_edges_shorter_than_a_tenth_of_a_mutation_without_state_change(
    #[case] (length, estimate, indels): (f64, Vec<Sub>, Vec<InDel>),
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let mut case = ShortEdgeCase::new(length, estimate, indels)?;

    case.collapse()?;

    assert_eq!(expected, case.topology()?);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::optimize::__tests__::test_topology_cleanup::tests::helpers::{empty_sparse_recon, populate_test_nodes};
    use crate::partition::set::PartitionSet;
    use crate::test_utils::{deletion, find_edge_key, fitch_edge_obs, sparse_edge_obs};
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::EdgeMap;
    use treetime_graph::node_names::NodeNames;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read, nwk_write_str};
    use treetime_primitives::AsciiChar;

    pub(super) fn gap_change() -> InDel {
      deletion((3, 5), [AsciiChar::from_byte_unchecked(b'A'); 2].as_slice())
    }

    pub(super) struct ShortEdgeCase {
      graph: Graph,
      names: NodeNames,
      partitions: Option<PartitionSet>,
      branch_lengths: EdgeMap<Option<f64>>,
    }

    impl ShortEdgeCase {
      pub(super) fn new(length: f64, estimate: Vec<Sub>, indels: Vec<InDel>) -> Result<Self, Report> {
        let nwk_parsed = nwk_read(b"((A:0.1,B:0.1)P:0.1,C:0.1,D:0.1)root;".as_slice())?;
        let graph = nwk_parsed.graph;
        let names = nwk_parsed.names;
        let mut branch_lengths = nwk_parsed.branch_lengths;
        let mut reconstruction = empty_sparse_recon()?;
        populate_test_nodes(&mut reconstruction, &graph);
        for edge in graph.get_edges() {
          reconstruction
            .partition
            .obs_edges
            .insert(edge.key(), fitch_edge_obs(vec![]));
          reconstruction.edges.estimates.insert(edge.key(), vec![]);
        }
        let short_edge = find_edge_key(&graph, &names, "root", "P").unwrap();
        reconstruction
          .partition
          .obs_edges
          .insert(short_edge, sparse_edge_obs(estimate.clone(), indels));
        reconstruction.edges.estimates.insert(short_edge, estimate);
        branch_lengths.insert(short_edge, Some(length));
        Ok(Self {
          graph,
          names,
          partitions: Some(PartitionSet::Sparse(vec![reconstruction])),
          branch_lengths,
        })
      }

      pub(super) fn collapse(&mut self) -> Result<(), Report> {
        let collapsed = collapse_short_internal_edges(
          &mut self.graph,
          self.partitions.take().unwrap(),
          &mut self.branch_lengths,
        )?;
        self.partitions = Some(collapsed.partitions);
        Ok(())
      }

      pub(super) fn topology(&self) -> Result<String, Report> {
        let no_lengths: EdgeMap<Option<f64>> = self.graph.get_edges().map(|edge| (edge.key(), None)).collect();
        nwk_write_str(
          &GraphSchedule::new(&self.graph)?.into_tree()?,
          &self.names,
          &no_lengths,
          &NwkWriteOptions::default(),
          &NwkNodeComments::new(),
        )
      }
    }
  }
}

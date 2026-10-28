#[cfg(test)]
mod tests {
  use crate::optimize::__tests__::test_topology_cleanup::tests::helpers::sub;
  use crate::optimize::indel::EdgeIndels;
  use crate::optimize::topology::collapse::collapse_short_internal_edges;
  use crate::seq::mutation::Sub;
  use eyre::Report;
  use helpers::ShortEdgeCase;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  const ONE_MUTATION: f64 = 1.0 / 100.0;

  #[rustfmt::skip]
  #[rstest]
  #[case::short_without_substitution(0.05 * ONE_MUTATION, vec![],                   "(C,D,A,B)root;"    )]
  #[case::short_with_substitution(   0.05 * ONE_MUTATION, vec![sub(b'A', 0, b'T')], "((A,B)P,C,D)root;" )]
  #[case::long_without_substitution( 0.2 * ONE_MUTATION,  vec![],                   "((A,B)P,C,D)root;" )]
  #[case::zero_with_substitution(    0.0,                 vec![sub(b'A', 0, b'T')], "(C,D,A,B)root;"    )]
  #[trace]
  fn test_topology_cleanup_collapses_internal_edges_shorter_than_a_tenth_of_a_mutation_without_substitutions(
    #[case] length: f64,
    #[case] estimate: Vec<Sub>,
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let mut case = ShortEdgeCase::new(length, estimate)?;

    case.collapse(&EdgeIndels::Excluded)?;

    assert_eq!(expected, case.topology()?);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::optimize::__tests__::test_topology_cleanup::tests::helpers::{empty_sparse_recon, populate_test_nodes};
    use crate::partition::set::PartitionSet;
    use crate::test_utils::{find_edge_key, fitch_edge_obs};
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::EdgeMap;
    use treetime_graph::node_names::NodeNames;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read, nwk_write_str};

    pub(super) struct ShortEdgeCase {
      graph: Graph,
      names: NodeNames,
      partitions: Option<PartitionSet>,
      branch_lengths: EdgeMap<Option<f64>>,
    }

    impl ShortEdgeCase {
      pub(super) fn new(length: f64, estimate: Vec<Sub>) -> Result<Self, Report> {
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
          .insert(short_edge, fitch_edge_obs(estimate.clone()));
        reconstruction.edges.estimates.insert(short_edge, estimate);
        branch_lengths.insert(short_edge, Some(length));
        Ok(Self {
          graph,
          names,
          partitions: Some(PartitionSet::Sparse(vec![reconstruction])),
          branch_lengths,
        })
      }

      pub(super) fn collapse(&mut self, indels: &EdgeIndels) -> Result<(), Report> {
        let partitions = collapse_short_internal_edges(
          &mut self.graph,
          self.partitions.take().unwrap(),
          indels,
          &mut self.branch_lengths,
        )?;
        self.partitions = Some(partitions);
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

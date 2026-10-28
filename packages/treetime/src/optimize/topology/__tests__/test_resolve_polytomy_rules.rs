#[cfg(test)]
mod tests {
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::params::TopologyOps;
  use crate::optimize::topology::merge_shared_mutations::MergeLedger;
  use crate::optimize::topology::resolve_polytomy::resolve_polytomies;
  use eyre::Report;
  use helpers::{EdgeMutations, ResolveCase, sub};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::indels_excluded(IndelPolicy::Exclude, "(A,B,C)root;"          )]
  #[case::indels_included(IndelPolicy::Include, "(C,(A,B)NODE_0000000)root;" )]
  #[trace]
  fn test_resolve_polytomy_merges_siblings_sharing_only_an_indel_only_when_indels_are_included(
    #[case] indel_policy: IndelPolicy,
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let mut case = ResolveCase::new(
      "(A:0.1,B:0.1,C:0.1)root;",
      &[
        EdgeMutations::new("root", "A", vec![]).with_deletion(20),
        EdgeMutations::new("root", "B", vec![]).with_deletion(20),
        EdgeMutations::new("root", "C", vec![sub(b'G', 5, b'C')]),
      ],
    )?;

    case.resolve(indel_policy)?;

    assert_eq!(expected, case.topology()?);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::indels_excluded(IndelPolicy::Exclude, "((S,(C2,C3)V,C1)U,T)root;"         )]
  #[case::indels_included(IndelPolicy::Include, "((S,((C2,C3)V,C1)NODE_0000000)U,T)root;" )]
  #[trace]
  fn test_resolve_polytomy_retires_a_helper_edge_with_only_indels_only_when_indels_are_excluded(
    #[case] indel_policy: IndelPolicy,
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let mut case = ResolveCase::new(
      "(((C1:0.1,C2:0.1,C3:0.1)V:0.2,S:0.1)U:0.1,T:0.1)root;",
      &[
        EdgeMutations::new("U", "V", vec![sub(b'A', 0, b'T')]).with_deletion(20),
        EdgeMutations::new("V", "C1", vec![sub(b'T', 0, b'A')]),
        EdgeMutations::new("V", "C2", vec![sub(b'G', 5, b'C')]),
        EdgeMutations::new("V", "C3", vec![sub(b'C', 10, b'G')]),
      ],
    )?;

    case.resolve(indel_policy)?;

    assert_eq!(expected, case.topology()?);
    Ok(())
  }

  #[test]
  fn test_resolve_polytomy_keeps_a_polytomy_whose_children_all_share_a_substitution() -> Result<(), Report> {
    let mut case = ResolveCase::new(
      "(A:0.1,B:0.1,C:0.1)root;",
      &[
        EdgeMutations::new("root", "A", vec![sub(b'A', 0, b'T')]),
        EdgeMutations::new("root", "B", vec![sub(b'A', 0, b'T')]),
        EdgeMutations::new("root", "C", vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')]),
      ],
    )?;

    let changed = case.resolve(IndelPolicy::Include)?;

    assert_eq!((0, "(A,B,C)root;"), (changed, case.topology()?.as_str()));
    Ok(())
  }

  #[test]
  fn test_resolve_polytomy_hoist_from_a_bifurcation_leaves_no_single_child_node() -> Result<(), Report> {
    let mut case = ResolveCase::new(
      "(((C1:0.1,C2:0.1,C3:0.1)V:0.2,S:0.1)U:0.1,T:0.1)root;",
      &[
        EdgeMutations::new("U", "V", vec![sub(b'A', 0, b'T'), sub(b'C', 5, b'G')]),
        EdgeMutations::new("V", "C1", vec![sub(b'T', 0, b'A')]),
        EdgeMutations::new("V", "C2", vec![sub(b'T', 0, b'A')]),
      ],
    )?;

    case.resolve(IndelPolicy::Include)?;

    assert_eq!("((S,(C3,C1,C2)NODE_0000001)U,T)root;", case.topology()?);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::partition::marginal::reconstruction::SparseReconstruction;
    use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
    use crate::partition::storage::sparse::{SparseNodeObs, SparseNodeState};
    use crate::seq::mutation::Sub;
    use crate::test_utils::{deletion, find_edge_key, sparse_edge_obs, sparse_topology_partition};
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::EdgeMap;
    use treetime_graph::keymap;
    use treetime_graph::node_names::NodeNames;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read, nwk_write_str};
    use treetime_primitives::{AsciiChar, Seq};

    const LENGTH: usize = 100;

    pub(super) fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
      Sub::new(
        AsciiChar::from_byte_unchecked(reff),
        pos,
        AsciiChar::from_byte_unchecked(qry),
      )
      .unwrap()
    }

    pub(super) struct EdgeMutations {
      source: &'static str,
      target: &'static str,
      subs: Vec<Sub>,
      deletion_start: Option<usize>,
    }

    impl EdgeMutations {
      pub(super) fn new(source: &'static str, target: &'static str, subs: Vec<Sub>) -> Self {
        Self {
          source,
          target,
          subs,
          deletion_start: None,
        }
      }

      pub(super) fn with_deletion(self, start: usize) -> Self {
        Self {
          deletion_start: Some(start),
          ..self
        }
      }
    }

    pub(super) struct ResolveCase {
      graph: Graph,
      names: NodeNames,
      sparse: Vec<SparseReconstruction>,
      branch_lengths: EdgeMap<Option<f64>>,
    }

    impl ResolveCase {
      pub(super) fn new(newick: &str, edges: &[EdgeMutations]) -> Result<Self, Report> {
        let nwk_parsed = nwk_read(newick.as_bytes())?;
        let graph = nwk_parsed.graph;
        let names = nwk_parsed.names;
        let alphabet = Alphabet::new(AlphabetName::Nuc)?;
        let root_sequence: Seq = std::iter::repeat_with(|| AsciiChar::from_byte_unchecked(b'A'))
          .take(LENGTH)
          .collect();

        let mut obs_nodes = keymap! {};
        let mut node_states = keymap! {};
        for node in graph.get_nodes() {
          obs_nodes.insert(node.key(), SparseNodeObs::empty(&alphabet));
          node_states.insert(node.key(), SparseNodeState::leaf(&root_sequence));
        }
        let mut obs_edges = keymap! {};
        for edge in edges {
          let edge_key = find_edge_key(&graph, &names, edge.source, edge.target).unwrap();
          let indels = edge.deletion_start.map_or(vec![], |start| {
            vec![deletion((start, start + 3), Seq::try_from_str("GTA").unwrap())]
          });
          obs_edges.insert(edge_key, sparse_edge_obs(edge.subs.clone(), indels));
        }
        let partition = PartitionMarginalSparse {
          alphabet,
          length: LENGTH,
          total_length: LENGTH,
          root_sequence,
          obs_nodes,
          obs_edges,
          outside_at_no_data: false,
        };

        Ok(Self {
          graph,
          names,
          sparse: vec![sparse_topology_partition(partition, node_states)],
          branch_lengths: nwk_parsed.branch_lengths,
        })
      }

      pub(super) fn resolve(&mut self, indel_policy: IndelPolicy) -> Result<usize, Report> {
        let changed = resolve_polytomies(
          &mut self.graph,
          &mut self.names,
          &mut self.sparse,
          TopologyOps::default(),
          &mut self.branch_lengths,
          &mut MergeLedger::new(indel_policy),
        )?;
        self.graph.build()?;
        Ok(changed)
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

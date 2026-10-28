#[cfg(test)]
mod tests {
  use crate::optimize::__tests__::test_topology_cleanup::tests::helpers::sub;
  use eyre::Report;
  use helpers::LedgerCase;
  use pretty_assertions::assert_eq;
  use treetime_utils::o;

  #[test]
  fn test_optimize_merge_ledger_skips_a_group_whose_helper_edge_collapsed() -> Result<(), Report> {
    let mut case = LedgerCase::new()?;

    let after_merge = case.cleanup(&[])?;
    let after_collapse = case.cleanup(&[("P", "NODE_0000000")])?;

    let expected = (o!("((C,(A,B)NODE_0000000)P,E,F)G;"), o!("((C,A,B)P,E,F)G;"));
    assert_eq!(expected, (after_merge, after_collapse));
    Ok(())
  }

  #[test]
  fn test_optimize_merge_ledger_skips_the_group_under_the_grandparent_after_the_parent_collapsed() -> Result<(), Report>
  {
    let mut case = LedgerCase::new()?;
    case.cleanup(&[])?;
    case.cleanup(&[("P", "NODE_0000000")])?;

    let after_parent_collapse = case.cleanup(&[("G", "P")])?;

    assert_eq!("(E,F,C,A,B)G;", after_parent_collapse);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::optimize::__tests__::test_topology_cleanup::tests::helpers::{empty_sparse_recon, populate_test_nodes};
    use crate::optimize::indel::IndelPolicy;
    use crate::optimize::params::TopologyOps;
    use crate::optimize::run_loop::prune_and_merge_in_loop;
    use crate::optimize::topology::merge_shared_mutations::MergeLedger;
    use crate::partition::set::PartitionSet;
    use crate::partition::storage::sparse::SparseEdgeObs;
    use crate::test_utils::{find_edge_key, fitch_edge_obs};
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::EdgeMap;
    use treetime_graph::node_names::NodeNames;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read, nwk_write_str};

    pub(super) struct LedgerCase {
      graph: Graph,
      names: NodeNames,
      partitions: Option<PartitionSet>,
      branch_lengths: EdgeMap<Option<f64>>,
      ledger: MergeLedger,
    }

    impl LedgerCase {
      pub(super) fn new() -> Result<Self, Report> {
        let nwk_parsed = nwk_read(b"((A:0.1,B:0.1,C:0.1)P:0.1,E:0.1,F:0.1)G;".as_slice())?;
        let graph = nwk_parsed.graph;
        let names = nwk_parsed.names;
        let mut reconstruction = empty_sparse_recon()?;
        populate_test_nodes(&mut reconstruction, &graph);
        let edges = [
          ("G", "P", vec![]),
          ("P", "A", vec![sub(b'A', 0, b'T')]),
          ("P", "B", vec![sub(b'A', 0, b'T')]),
          ("P", "C", vec![sub(b'G', 5, b'C')]),
          ("G", "E", vec![sub(b'C', 10, b'G')]),
          ("G", "F", vec![sub(b'T', 15, b'A')]),
        ];
        for (source, target, subs) in edges {
          let edge_key = find_edge_key(&graph, &names, source, target).unwrap();
          let edge = if subs.is_empty() {
            SparseEdgeObs::default()
          } else {
            fitch_edge_obs(subs)
          };
          reconstruction.partition.obs_edges.insert(edge_key, edge);
        }
        Ok(Self {
          graph,
          names,
          partitions: Some(PartitionSet::Sparse(vec![reconstruction])),
          branch_lengths: nwk_parsed.branch_lengths,
          ledger: MergeLedger::new(IndelPolicy::Include),
        })
      }

      pub(super) fn cleanup(&mut self, zero_edges: &[(&str, &str)]) -> Result<String, Report> {
        let zero_edges = zero_edges
          .iter()
          .map(|(source, target)| find_edge_key(&self.graph, &self.names, source, target).unwrap())
          .collect::<Vec<_>>();
        let cleanup = prune_and_merge_in_loop(
          &mut self.graph,
          self.partitions.take().unwrap(),
          &zero_edges,
          TopologyOps::default(),
          &mut self.branch_lengths,
          &mut self.names,
          &mut self.ledger,
        )?;
        self.partitions = Some(cleanup.partitions);
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

#[cfg(test)]
mod tests {
  use crate::optimize::params::TopologyOps;
  use crate::optimize::topology::resolve_polytomy::resolve_polytomies;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::mutation::Sub;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;

  use helpers::{no_dense, reversion_present, sub, total_subs};
  use treetime_io::nwk::nwk_read_str;

  // root -> U -> V -> {C1, C2, C3}. V is the polytomy under test.
  const NWK: &str = "(((C1:0.1,C2:0.1,C3:0.1)V:0.2)U:0.1)root:0.0;";

  #[test]
  fn test_resolve_polytomy_merge_hoist_retire_worked_example() -> Result<(), Report> {
    // M_v = {A0T (p), C5G (q)}; C1 and C2 both revert p, C3 keeps it. The routine merges
    // C1+C2, hoists the reverting group, and retires the helper, reaching the parsimony
    // optimum of two mutations (q above, p only on the C3 lineage).
    let mut graph: GraphAncestral = nwk_read_str(NWK)?;
    let partition = helpers::make_partition(
      &graph,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T'), sub(b'C', 5, b'G')]),
        ("V", "C1", vec![sub(b'T', 0, b'A')]),
        ("V", "C2", vec![sub(b'T', 0, b'A')]),
        ("V", "C3", vec![]),
      ],
    );
    let sparse = vec![partition];

    let changed = resolve_polytomies(&mut graph, &sparse, &no_dense(), TopologyOps::default())?;
    assert!(changed > 0);

    let p = sparse[0].read_arc();
    assert_eq!(total_subs(&graph, &p), 2);
    assert!(!reversion_present(&graph, &p, &sub(b'T', 0, b'A')));

    for leaf in ["C1", "C2", "C3"] {
      assert!(
        find_node_key_by_name(&graph, leaf).is_some(),
        "leaf {leaf} must survive"
      );
    }
    Ok(())
  }

  #[test]
  fn test_resolve_polytomy_incompatible_splits_five_to_four() -> Result<(), Report> {
    // C1 reverts p1, C2 reverts p2 (different positions): the two required splits are
    // incompatible. One hoist takes the total from 5 to 4; the residual reversion is
    // irreducible homoplasy, and the routine stops there.
    let mut graph: GraphAncestral = nwk_read_str(NWK)?;
    let partition = helpers::make_partition(
      &graph,
      0,
      100,
      &[
        (
          "U",
          "V",
          vec![sub(b'A', 0, b'T'), sub(b'C', 5, b'G'), sub(b'G', 10, b'A')],
        ),
        ("V", "C1", vec![sub(b'T', 0, b'A')]),
        ("V", "C2", vec![sub(b'G', 5, b'C')]),
        ("V", "C3", vec![]),
      ],
    );
    let sparse = vec![partition];

    let before = total_subs(&graph, &sparse[0].read_arc());
    resolve_polytomies(&mut graph, &sparse, &no_dense(), TopologyOps::default())?;
    let after = total_subs(&graph, &sparse[0].read_arc());

    assert_eq!(before, 5);
    assert_eq!(after, 4);
    Ok(())
  }

  #[test]
  fn test_resolve_polytomy_retirement_preserves_preexisting_internal_node() -> Result<(), Report> {
    // W is a pre-existing internal node reached by a mutation-free edge from V. Helper
    // retirement must dissolve only nodes it created, never W, even though V->W is empty.
    let mut graph: GraphAncestral = nwk_read_str("((((X1:0.1,X2:0.1)W:0.0,C1:0.1,C2:0.1)V:0.2)U:0.1)root:0.0;")?;
    let partition = helpers::make_partition(
      &graph,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T')]),
        ("V", "C1", vec![sub(b'T', 0, b'A')]),
        ("V", "C2", vec![sub(b'T', 0, b'A')]),
        ("V", "W", vec![]),
      ],
    );
    let sparse = vec![partition];

    resolve_polytomies(&mut graph, &sparse, &no_dense(), TopologyOps::default())?;

    assert!(
      find_node_key_by_name(&graph, "W").is_some(),
      "pre-existing internal node W must survive helper retirement"
    );
    for leaf in ["X1", "X2", "C1", "C2"] {
      assert!(
        find_node_key_by_name(&graph, leaf).is_some(),
        "leaf {leaf} must survive"
      );
    }
    let p = sparse[0].read_arc();
    assert_eq!(total_subs(&graph, &p), 1);
    assert!(!reversion_present(&graph, &p, &sub(b'T', 0, b'A')));
    Ok(())
  }

  #[test]
  fn test_resolve_polytomy_root_polytomy_skipped() -> Result<(), Report> {
    // A polytomy at the root has no parent edge to revert, so no hoist fires. With no
    // shared substitutions there is nothing to do; the routine leaves the tree untouched.
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.1,C:0.1)root;")?;
    let partition = helpers::make_partition(
      &graph,
      0,
      100,
      &[
        ("root", "A", vec![sub(b'A', 0, b'T')]),
        ("root", "B", vec![sub(b'C', 5, b'G')]),
        ("root", "C", vec![sub(b'G', 10, b'A')]),
      ],
    );
    let sparse = vec![partition];
    let nodes_before = graph.get_nodes().len();

    let changed = resolve_polytomies(&mut graph, &sparse, &no_dense(), TopologyOps::default())?;

    assert_eq!(changed, 0);
    assert_eq!(graph.get_nodes().len(), nodes_before);
    assert_eq!(total_subs(&graph, &sparse[0].read_arc()), 3);
    Ok(())
  }

  #[test]
  fn test_resolve_polytomy_no_change_without_reversions() -> Result<(), Report> {
    // Distinct, non-shared, non-reverting child substitutions: nothing to merge or hoist.
    let mut graph: GraphAncestral = nwk_read_str(NWK)?;
    let partition = helpers::make_partition(
      &graph,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T')]),
        ("V", "C1", vec![sub(b'C', 5, b'G')]),
        ("V", "C2", vec![sub(b'G', 10, b'A')]),
        ("V", "C3", vec![sub(b'T', 15, b'A')]),
      ],
    );
    let sparse = vec![partition];
    let nodes_before = graph.get_nodes().len();

    let changed = resolve_polytomies(&mut graph, &sparse, &no_dense(), TopologyOps::default())?;

    assert_eq!(changed, 0);
    assert_eq!(graph.get_nodes().len(), nodes_before);
    assert_eq!(total_subs(&graph, &sparse[0].read_arc()), 4);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::partition::storage::sparse::{SparseEdgePartition, SparseNodePartition};
    use crate::test_utils::find_edge_key;
    use maplit::btreemap;
    use treetime_primitives::{AsciiChar, Seq, seq};

    pub fn c(b: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(b)
    }

    pub fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
      Sub::new(c(reff), pos, c(qry)).unwrap()
    }

    pub fn no_dense() -> Vec<Arc<RwLock<PartitionMarginalDense>>> {
      vec![]
    }

    pub fn total_subs(graph: &GraphAncestral, partition: &PartitionMarginalSparse) -> usize {
      graph
        .get_edges()
        .iter()
        .filter_map(|e| partition.edges.get(&e.read_arc().key()))
        .map(|e| e.fitch_subs().len())
        .sum()
    }

    pub fn reversion_present(graph: &GraphAncestral, partition: &PartitionMarginalSparse, needle: &Sub) -> bool {
      graph
        .get_edges()
        .iter()
        .filter_map(|e| partition.edges.get(&e.read_arc().key()))
        .any(|e| e.fitch_subs().contains(needle))
    }

    pub fn make_partition(
      graph: &GraphAncestral,
      index: usize,
      length: usize,
      edge_mutations: &[(&str, &str, Vec<Sub>)],
    ) -> Arc<RwLock<PartitionMarginalSparse>> {
      let mut partition = PartitionMarginalSparse {
        index,
        gtr: jc69(JC69Params::default()).unwrap(),
        alphabet: Alphabet::new(AlphabetName::Nuc).unwrap(),
        length,
        nodes: btreemap! {},
        edges: btreemap! {},
        root_sequence: seq![],
      };

      let mut ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
      for (_, _, subs) in edge_mutations {
        for s in subs {
          if s.pos() < length {
            ref_seq[s.pos()] = s.reff();
          }
        }
      }
      partition.root_sequence = ref_seq.clone();

      for node in graph.get_nodes() {
        let key = node.read_arc().key();
        let mut node_part = SparseNodePartition::empty(&partition.alphabet);
        node_part.seq.sequence = ref_seq.clone();
        partition.nodes.insert(key, node_part);
      }

      for (source, target, subs) in edge_mutations {
        let edge_key =
          find_edge_key(graph, source, target).unwrap_or_else(|| panic!("edge {source}->{target} missing"));
        partition
          .edges
          .insert(edge_key, SparseEdgePartition::with_fitch_subs(subs.clone()));
      }

      Arc::new(RwLock::new(partition))
    }
  }
}

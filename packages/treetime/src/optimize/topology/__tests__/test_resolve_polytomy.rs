#[cfg(test)]
mod tests {
  use crate::optimize::params::TopologyOps;
  use crate::optimize::topology::resolve_polytomy::resolve_polytomies;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::seq::mutation::Sub;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;

  use helpers::{reversion_present, sub, total_subs};
  use treetime_io::nwk::nwk_read_str;

  // root -> U -> V -> {C1, C2, C3}. V is the polytomy under test.
  const NWK: &str = "(((C1:0.1,C2:0.1,C3:0.1)V:0.2)U:0.1)root:0.0;";

  #[test]
  fn test_resolve_polytomy_merge_hoist_retire_worked_example() -> Result<(), Report> {
    // M_v = {A0T (p), C5G (q)}; C1 and C2 both revert p, C3 keeps it. The routine merges
    // C1+C2, hoists the reverting group, and retires the helper, reaching the parsimony
    // optimum of two mutations (q above, p only on the C3 lineage).
    let nwk_parsed = nwk_read_str(NWK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (partition, node_states) = helpers::make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T'), sub(b'C', 5, b'G')]),
        ("V", "C1", vec![sub(b'T', 0, b'A')]),
        ("V", "C2", vec![sub(b'T', 0, b'A')]),
        ("V", "C3", vec![]),
      ],
    );
    let mut sparse = vec![partition];
    let mut node_states = vec![node_states];

    let mut branch_lengths = branch_lengths;
    let changed = resolve_polytomies(
      &mut graph,
      &mut sparse,
      &mut node_states,
      TopologyOps::default(),
      &mut branch_lengths,
    )?;
    assert!(changed > 0);

    let p = &sparse[0];
    assert_eq!(total_subs(&graph, p), 2);
    assert!(!reversion_present(&graph, p, &sub(b'T', 0, b'A')));

    for leaf in ["C1", "C2", "C3"] {
      assert!(
        find_node_key_by_name(&graph, &names, leaf).is_some(),
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
    let nwk_parsed = nwk_read_str(NWK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (partition, node_states) = helpers::make_partition(
      &graph,
      &names,
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
    let mut sparse = vec![partition];
    let mut node_states = vec![node_states];

    let before = total_subs(&graph, &sparse[0]);
    let mut branch_lengths = branch_lengths;
    resolve_polytomies(
      &mut graph,
      &mut sparse,
      &mut node_states,
      TopologyOps::default(),
      &mut branch_lengths,
    )?;
    let after = total_subs(&graph, &sparse[0]);

    assert_eq!(before, 5);
    assert_eq!(after, 4);
    Ok(())
  }

  #[test]
  fn test_resolve_polytomy_bifurcating_root_cross_root_reversion() -> Result<(), Report> {
    // Bifurcating root: root -> {V, S}. At position 0 the majority under V is A, so V's parent
    // edge is empty there; the distinguishing G sits on the sibling edge root->S and on V's two
    // G-state children G1, G2, which the arbitrary root placement scatters as homoplasy. The
    // per-node scan sees nothing on V's own parent edge, but looking across the bifurcating root
    // exposes the reversion: the routine hoists the G-group above V and reaches the parsimony
    // optimum of one mutation (a single A->G change separating the A-clade {A1, A2}).
    //
    // Oracle: position 0 partitions the unrooted tree into one clade, so one change is the
    // analytical parsimony minimum, the same value the routine reaches on any rooting of this
    // tree (compare test_resolve_polytomy_merge_hoist_retire_worked_example, where the same
    // reversion sits on a genuine internal edge).
    let nwk_parsed = nwk_read_str("((G1:0.1,G2:0.1,A1:0.1,A2:0.1)V:0.1,S:0.1)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (partition, node_states) = helpers::make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        ("root", "S", vec![sub(b'A', 0, b'G')]),
        ("V", "G1", vec![sub(b'A', 0, b'G')]),
        ("V", "G2", vec![sub(b'A', 0, b'G')]),
        ("V", "A1", vec![]),
        ("V", "A2", vec![]),
      ],
    );
    let mut sparse = vec![partition];
    let mut node_states = vec![node_states];

    let before = total_subs(&graph, &sparse[0]);
    let mut branch_lengths = branch_lengths;
    let changed = resolve_polytomies(
      &mut graph,
      &mut sparse,
      &mut node_states,
      TopologyOps::default(),
      &mut branch_lengths,
    )?;
    let after = total_subs(&graph, &sparse[0]);

    assert_eq!(before, 3);
    assert!(changed > 0);
    assert_eq!(after, 1);
    assert!(!reversion_present(&graph, &sparse[0], &sub(b'A', 0, b'G')));

    for leaf in ["G1", "G2", "A1", "A2", "S"] {
      assert!(
        find_node_key_by_name(&graph, &names, leaf).is_some(),
        "leaf {leaf} must survive"
      );
    }
    Ok(())
  }

  #[test]
  fn test_resolve_polytomy_retirement_preserves_preexisting_internal_node() -> Result<(), Report> {
    // W is a pre-existing internal node reached by a mutation-free edge from V. Helper
    // retirement must dissolve only nodes it created, never W, even though V->W is empty.
    let nwk_parsed = nwk_read_str("((((X1:0.1,X2:0.1)W:0.0,C1:0.1,C2:0.1)V:0.2)U:0.1)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (partition, node_states) = helpers::make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T')]),
        ("V", "C1", vec![sub(b'T', 0, b'A')]),
        ("V", "C2", vec![sub(b'T', 0, b'A')]),
        ("V", "W", vec![]),
      ],
    );
    let mut sparse = vec![partition];
    let mut node_states = vec![node_states];

    let mut branch_lengths = branch_lengths;
    resolve_polytomies(
      &mut graph,
      &mut sparse,
      &mut node_states,
      TopologyOps::default(),
      &mut branch_lengths,
    )?;

    assert!(
      find_node_key_by_name(&graph, &names, "W").is_some(),
      "pre-existing internal node W must survive helper retirement"
    );
    for leaf in ["X1", "X2", "C1", "C2"] {
      assert!(
        find_node_key_by_name(&graph, &names, leaf).is_some(),
        "leaf {leaf} must survive"
      );
    }
    let p = &sparse[0];
    assert_eq!(total_subs(&graph, p), 1);
    assert!(!reversion_present(&graph, p, &sub(b'T', 0, b'A')));
    Ok(())
  }

  #[test]
  fn test_resolve_polytomy_root_polytomy_skipped() -> Result<(), Report> {
    // A polytomy at the root has no parent edge to revert, so no hoist fires. With no
    // shared substitutions there is nothing to do; the routine leaves the tree untouched.
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.1,C:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (partition, node_states) = helpers::make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        ("root", "A", vec![sub(b'A', 0, b'T')]),
        ("root", "B", vec![sub(b'C', 5, b'G')]),
        ("root", "C", vec![sub(b'G', 10, b'A')]),
      ],
    );
    let mut sparse = vec![partition];
    let mut node_states = vec![node_states];
    let nodes_before = graph.get_nodes().count();

    let mut branch_lengths = branch_lengths;
    let changed = resolve_polytomies(
      &mut graph,
      &mut sparse,
      &mut node_states,
      TopologyOps::default(),
      &mut branch_lengths,
    )?;

    assert_eq!(changed, 0);
    assert_eq!(graph.get_nodes().count(), nodes_before);
    assert_eq!(total_subs(&graph, &sparse[0]), 3);
    Ok(())
  }

  #[test]
  fn test_resolve_polytomy_no_change_without_reversions() -> Result<(), Report> {
    // Distinct, non-shared, non-reverting child substitutions: nothing to merge or hoist.
    let nwk_parsed = nwk_read_str(NWK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (partition, node_states) = helpers::make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T')]),
        ("V", "C1", vec![sub(b'C', 5, b'G')]),
        ("V", "C2", vec![sub(b'G', 10, b'A')]),
        ("V", "C3", vec![sub(b'T', 15, b'A')]),
      ],
    );
    let mut sparse = vec![partition];
    let mut node_states = vec![node_states];
    let nodes_before = graph.get_nodes().count();

    let mut branch_lengths = branch_lengths;
    let changed = resolve_polytomies(
      &mut graph,
      &mut sparse,
      &mut node_states,
      TopologyOps::default(),
      &mut branch_lengths,
    )?;

    assert_eq!(changed, 0);
    assert_eq!(graph.get_nodes().count(), nodes_before);
    assert_eq!(total_subs(&graph, &sparse[0]), 4);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};

    use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
    use crate::test_utils::find_edge_key;
    use maplit::btreemap;
    use treetime_primitives::{AsciiChar, Seq};

    pub fn c(b: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(b)
    }

    pub fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
      Sub::new(c(reff), pos, c(qry)).unwrap()
    }

    pub fn total_subs(graph: &Graph, recon: &PartitionMarginalSparse) -> usize {
      graph
        .get_edges()
        .filter_map(|e| recon.obs_edges.get(&e.key()))
        .map(|e| e.fitch_subs().len())
        .sum()
    }

    pub fn reversion_present(graph: &Graph, recon: &PartitionMarginalSparse, needle: &Sub) -> bool {
      graph
        .get_edges()
        .filter_map(|e| recon.obs_edges.get(&e.key()))
        .any(|e| e.fitch_subs().contains(needle))
    }

    pub fn make_partition(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      index: usize,
      length: usize,
      edge_mutations: &[(&str, &str, Vec<Sub>)],
    ) -> (PartitionMarginalSparse, BTreeMap<GraphNodeKey, SparseNodeState>) {
      let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();

      let mut ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
      for (_, _, subs) in edge_mutations {
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
        obs_nodes.insert(key, SparseNodeObs::empty(&alphabet));
        node_states.insert(key, SparseNodeState::leaf(&ref_seq));
      }

      let mut obs_edges = btreemap! {};
      for (source, target, subs) in edge_mutations {
        let edge_key =
          find_edge_key(graph, names, source, target).unwrap_or_else(|| panic!("edge {source}->{target} missing"));
        obs_edges.insert(edge_key, SparseEdgeObs::with_fitch_subs(subs.clone()));
      }

      let partition = PartitionMarginalSparse {
        index,
        alphabet,
        length,
        root_sequence: ref_seq,
        obs_nodes,
        obs_edges,
      };

      (partition, node_states)
    }
  }
}

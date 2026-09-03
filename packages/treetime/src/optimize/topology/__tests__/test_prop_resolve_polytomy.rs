#[cfg(test)]
mod tests {
  use crate::optimize::params::TopologyOps;
  use crate::optimize::topology::resolve_polytomy::resolve_polytomies;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::mutation::Sub;
  use parking_lot::RwLock;
  use proptest::prelude::*;
  use std::collections::BTreeSet;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;

  proptest! {
    /// The routine never increases the total mutation count, and when it changes anything
    /// the count strictly decreases: the first component of the (mutation count, node count)
    /// potential falls on every applied merge or hoist. This also exercises termination -
    /// a non-decreasing potential would hang the test.
    #[test]
    fn test_prop_resolve_polytomy_potential_decreases(
      n_children in 3_usize..7,
      k in 1_usize..5,
      revert_masks in prop::collection::vec(0_u32..32, 3..7),
      own_counts in prop::collection::vec(0_usize..3, 3..7),
    ) {
      let (mut graph, partition, before) = helpers::build_case(n_children, k, &revert_masks, &own_counts);
      let sparse = vec![partition];

      let changed = resolve_polytomies(&mut graph, &sparse, &helpers::no_dense(), TopologyOps::default()).unwrap();
      let after = helpers::total_subs(&graph, &sparse[0].read_arc());

      prop_assert!(after <= before, "mutation count increased: before={before} after={after}");
      if changed > 0 {
        prop_assert!(after < before, "changed but mutation count did not fall: before={before} after={after}");
      }
    }

    /// Structural invariants hold after resolution: leaves are preserved, branch lengths stay
    /// non-negative, and the result is still a single-rooted tree (each non-root node keeps
    /// exactly one parent edge).
    #[test]
    fn test_prop_resolve_polytomy_preserves_tree(
      n_children in 3_usize..7,
      k in 1_usize..5,
      revert_masks in prop::collection::vec(0_u32..32, 3..7),
      own_counts in prop::collection::vec(0_usize..3, 3..7),
    ) {
      let (mut graph, partition, _before) = helpers::build_case(n_children, k, &revert_masks, &own_counts);
      let leaves_before = helpers::leaf_names(&graph);
      let sparse = vec![partition];

      resolve_polytomies(&mut graph, &sparse, &helpers::no_dense(), TopologyOps::default()).unwrap();

      prop_assert_eq!(helpers::leaf_names(&graph), leaves_before);

      for edge in graph.get_edges() {
        if let Some(bl) = edge.read_arc().payload().read_arc().branch_length() {
          prop_assert!(bl >= 0.0, "negative branch length {bl}");
        }
      }

      let mut roots = 0;
      for node in graph.get_nodes() {
        let inbound = node.read_arc().inbound().len();
        if inbound == 0 {
          roots += 1;
        } else {
          prop_assert_eq!(inbound, 1, "non-root node has {} parents", inbound);
        }
      }
      prop_assert_eq!(roots, 1, "expected exactly one root");
    }
  }

  mod helpers {
    use super::*;
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::partition::marginal::dense::partition::PartitionMarginalDense;
    use crate::partition::storage::sparse::{SparseEdgePartition, SparseNodePartition};
    use crate::test_utils::find_edge_key;
    use itertools::Itertools;
    use maplit::btreemap;
    use treetime_io::nwk::nwk_read_str;
    use treetime_primitives::{AsciiChar, Seq, seq};

    fn c(b: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(b)
    }

    pub fn no_dense() -> Vec<Arc<RwLock<PartitionMarginalDense>>> {
      vec![]
    }

    pub fn leaf_names(graph: &GraphAncestral) -> BTreeSet<String> {
      graph
        .get_nodes()
        .iter()
        .filter(|n| n.read_arc().is_leaf())
        .filter_map(|n| n.read_arc().payload().read_arc().name.clone())
        .collect()
    }

    pub fn total_subs(graph: &GraphAncestral, partition: &PartitionMarginalSparse) -> usize {
      graph
        .get_edges()
        .iter()
        .filter_map(|e| partition.edges.get(&e.read_arc().key()))
        .map(|e| e.fitch_subs().len())
        .sum()
    }

    /// Build a `root -> U -> V -> {C0..}` case. The parent edge U->V carries `k`
    /// substitutions `A->C` at positions `0..k`. Each child reverts the M_v positions set in
    /// its mask and adds its own distinct substitutions at high positions. Returns the graph,
    /// the sparse partition, and the initial total mutation count.
    pub fn build_case(
      n_children: usize,
      k: usize,
      revert_masks: &[u32],
      own_counts: &[usize],
    ) -> (GraphAncestral, Arc<RwLock<PartitionMarginalSparse>>, usize) {
      let n_children = n_children.min(revert_masks.len()).min(own_counts.len()).max(3);
      let length = 200_usize;

      let child_names: Vec<String> = (0..n_children).map(|i| format!("C{i}")).collect();
      let inner = child_names.iter().map(|name| format!("{name}:0.1")).join(",");
      let newick = format!("((({inner})V:0.2)U:0.1)root:0.0;");
      let graph: GraphAncestral = nwk_read_str(&newick).unwrap();

      // M_v: k substitutions A->C at positions 0..k.
      let parent_subs: Vec<Sub> = (0..k).map(|pos| Sub::new(c(b'A'), pos, c(b'C')).unwrap()).collect();

      let mut edge_mutations: Vec<(String, String, Vec<Sub>)> = vec![("U".to_owned(), "V".to_owned(), parent_subs)];

      let mut own_pos = k + 1;
      for (i, name) in child_names.iter().enumerate() {
        let mask = revert_masks[i];
        let mut subs: Vec<Sub> = Vec::new();
        for pos in 0..k {
          if mask & (1 << pos) != 0 {
            // Revert the parent substitution: C->A.
            subs.push(Sub::new(c(b'C'), pos, c(b'A')).unwrap());
          }
        }
        for _ in 0..own_counts[i] {
          subs.push(Sub::new(c(b'G'), own_pos, c(b'T')).unwrap());
          own_pos += 2;
        }
        subs.sort_by_key(Sub::pos);
        edge_mutations.push(("V".to_owned(), name.clone(), subs));
      }

      let total: usize = edge_mutations.iter().map(|(_, _, subs)| subs.len()).sum();
      let partition = make_partition(&graph, length, &edge_mutations);
      (graph, partition, total)
    }

    fn make_partition(
      graph: &GraphAncestral,
      length: usize,
      edge_mutations: &[(String, String, Vec<Sub>)],
    ) -> Arc<RwLock<PartitionMarginalSparse>> {
      let mut partition = PartitionMarginalSparse {
        index: 0,
        gtr: jc69(JC69Params::default()).unwrap(),
        alphabet: Alphabet::new(AlphabetName::Nuc).unwrap(),
        length,
        nodes: btreemap! {},
        edges: btreemap! {},
        root_sequence: seq![],
      };

      let ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
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

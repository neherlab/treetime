#[cfg(test)]
mod tests {
  use crate::optimize::params::TopologyOps;
  use crate::optimize::topology::resolve_polytomy::resolve_polytomies;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::seq::mutation::Sub;
  use proptest::prelude::*;
  use std::collections::BTreeMap;
  use std::collections::BTreeSet;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;

  proptest! {
    #[test]
    fn test_prop_resolve_polytomy_potential_decreases(
      n_children in 3_usize..7,
      k in 1_usize..5,
      revert_masks in prop::collection::vec(0_u32..32, 3..7),
      own_counts in prop::collection::vec(0_usize..3, 3..7),
    ) {
      let (mut graph, names, partition, node_states, before, mut branch_lengths) = helpers::build_case(n_children, k, &revert_masks, &own_counts);
      let mut sparse = vec![partition];
    let mut node_states = vec![node_states];

      let changed = resolve_polytomies(&mut graph, &mut sparse, &mut node_states, TopologyOps::default(), &mut branch_lengths).unwrap();
      let after = helpers::total_subs(&graph, &sparse[0]);

      prop_assert!(after <= before, "mutation count increased: before={before} after={after}");
      if changed > 0 {
        prop_assert!(after < before, "changed but mutation count did not fall: before={before} after={after}");
      }
    }

    #[test]
    fn test_prop_resolve_polytomy_preserves_tree(
      n_children in 3_usize..7,
      k in 1_usize..5,
      revert_masks in prop::collection::vec(0_u32..32, 3..7),
      own_counts in prop::collection::vec(0_usize..3, 3..7),
    ) {
      let (mut graph, names, partition, node_states, _before, mut branch_lengths) = helpers::build_case(n_children, k, &revert_masks, &own_counts);
      let leaves_before = helpers::leaf_names(&names, &graph);
      let mut sparse = vec![partition];
    let mut node_states = vec![node_states];

      resolve_polytomies(&mut graph, &mut sparse, &mut node_states, TopologyOps::default(), &mut branch_lengths).unwrap();

      prop_assert_eq!(helpers::leaf_names(&names, &graph), leaves_before);

      for bl in branch_lengths.values().flatten() {
        prop_assert!(*bl >= 0.0, "negative branch length {bl}");
      }

      let mut roots = 0;
      for node in graph.get_nodes() {
        let inbound = node.inbound().len();
        if inbound == 0 {
          roots += 1;
        } else {
          prop_assert_eq!(inbound, 1, "non-root node has {} parents", inbound);
        }
      }
      prop_assert_eq!(roots, 1, "expected exactly one root");
    }

    #[test]
    fn test_prop_resolve_polytomy_bifurcating_root_reduces_to_bipartition_cost(
      g in 1_usize..5,
      a in 2_usize..5,
      own_counts in prop::collection::vec(0_usize..3, 9),
    ) {
      let (mut graph, partition, node_states, before, expected_after, mut branch_lengths) =
        helpers::build_bifurcating_case(g, a, &own_counts);
      let mut sparse = vec![partition];
    let mut node_states = vec![node_states];

      let changed = resolve_polytomies(&mut graph, &mut sparse, &mut node_states, TopologyOps::default(), &mut branch_lengths).unwrap();
      let after = helpers::total_subs(&graph, &sparse[0]);

      prop_assert_eq!(after, expected_after, "did not reach the bipartition cost");
      prop_assert!(after <= before);
      if changed > 0 {
        prop_assert!(after < before, "changed but count did not fall: before={} after={}", before, after);
      }
    }
  }

  mod helpers {
    use super::*;
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};

    use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
    use crate::test_utils::find_edge_key;
    use itertools::Itertools;
    use maplit::btreemap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_io::nwk::nwk_read_str;
    use treetime_primitives::{AsciiChar, Seq};

    fn c(b: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(b)
    }

    pub fn leaf_names(names: &BTreeMap<GraphNodeKey, Option<String>>, graph: &Graph) -> BTreeSet<String> {
      graph
        .get_nodes()
        .filter(|n| n.is_leaf())
        .filter_map(|n| names.get(&n.key()).cloned().flatten())
        .collect()
    }

    pub fn total_subs(graph: &Graph, recon: &PartitionMarginalSparse) -> usize {
      graph
        .get_edges()
        .filter_map(|e| recon.obs_edges.get(&e.key()))
        .map(|e| e.fitch_subs().len())
        .sum()
    }

    pub fn build_case(
      n_children: usize,
      k: usize,
      revert_masks: &[u32],
      own_counts: &[usize],
    ) -> (
      Graph,
      BTreeMap<GraphNodeKey, Option<String>>,
      PartitionMarginalSparse,
      BTreeMap<GraphNodeKey, SparseNodeState>,
      usize,
      BTreeMap<GraphEdgeKey, Option<f64>>,
    ) {
      let n_children = n_children.min(revert_masks.len()).min(own_counts.len()).max(3);
      let length = 200_usize;

      let child_names: Vec<String> = (0..n_children).map(|i| format!("C{i}")).collect();
      let inner = child_names.iter().map(|name| format!("{name}:0.1")).join(",");
      let newick = format!("((({inner})V:0.2)U:0.1)root:0.0;");
      let nwk_parsed = nwk_read_str(&newick).unwrap();
      let node_names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;

      let parent_subs: Vec<Sub> = (0..k).map(|pos| Sub::new(c(b'A'), pos, c(b'C')).unwrap()).collect();

      let mut edge_mutations: Vec<(String, String, Vec<Sub>)> = vec![("U".to_owned(), "V".to_owned(), parent_subs)];

      let mut own_pos = k + 1;
      for (i, name) in child_names.iter().enumerate() {
        let mask = revert_masks[i];
        let mut subs: Vec<Sub> = Vec::new();
        for pos in 0..k {
          if mask & (1 << pos) != 0 {
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
      let (partition, node_states) = make_partition(&graph, &node_names, length, &edge_mutations);
      (graph, node_names, partition, node_states, total, branch_lengths)
    }

    pub fn build_bifurcating_case(
      g: usize,
      a: usize,
      own_counts: &[usize],
    ) -> (
      Graph,
      PartitionMarginalSparse,
      BTreeMap<GraphNodeKey, SparseNodeState>,
      usize,
      usize,
      BTreeMap<GraphEdgeKey, Option<f64>>,
    ) {
      let length = 200_usize;
      let g_names: Vec<String> = (0..g).map(|i| format!("G{i}")).collect();
      let a_names: Vec<String> = (0..a).map(|i| format!("A{i}")).collect();
      let children = g_names
        .iter()
        .chain(&a_names)
        .map(|name| format!("{name}:0.1"))
        .join(",");
      let newick = format!("(({children})V:0.1,S:0.1)root:0.0;");
      let nwk_parsed = nwk_read_str(&newick).unwrap();
      let node_names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;

      let mut edge_mutations: Vec<(String, String, Vec<Sub>)> = vec![(
        "root".to_owned(),
        "S".to_owned(),
        vec![Sub::new(c(b'A'), 0_usize, c(b'G')).unwrap()],
      )];

      let mut own_pos = 10_usize;
      let mut own_total = 0_usize;
      for (i, name) in g_names.iter().chain(&a_names).enumerate() {
        let is_g = i < g;
        let mut subs: Vec<Sub> = Vec::new();
        if is_g {
          subs.push(Sub::new(c(b'A'), 0_usize, c(b'G')).unwrap());
        }
        for _ in 0..own_counts[i % own_counts.len()] {
          subs.push(Sub::new(c(b'G'), own_pos, c(b'T')).unwrap());
          own_pos += 2;
          own_total += 1;
        }
        subs.sort_by_key(Sub::pos);
        edge_mutations.push(("V".to_owned(), name.clone(), subs));
      }

      let before: usize = edge_mutations.iter().map(|(_, _, subs)| subs.len()).sum();
      let (partition, node_states) = make_partition(&graph, &node_names, length, &edge_mutations);
      (graph, partition, node_states, before, 1 + own_total, branch_lengths)
    }

    fn make_partition(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      length: usize,
      edge_mutations: &[(String, String, Vec<Sub>)],
    ) -> (PartitionMarginalSparse, BTreeMap<GraphNodeKey, SparseNodeState>) {
      let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();

      let ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();

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
        index: 0,
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

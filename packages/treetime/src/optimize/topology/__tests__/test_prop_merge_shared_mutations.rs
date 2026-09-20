#![allow(
  clippy::collection_is_never_read,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;

  use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
  use crate::seq::indel::InDel;
  use crate::seq::mutation::Sub;
  use crate::test_utils::find_edge_key;
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use std::collections::BTreeMap;
  use std::collections::BTreeSet;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::{AsciiChar, Seq};

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn sub_at(pos: usize) -> Sub {
    let nucs = [b'A', b'C', b'G', b'T'];
    Sub::new(c(nucs[pos % 4]), pos, c(nucs[(pos + 1) % 4])).unwrap()
  }

  proptest! {
    #[test]
    fn test_prop_merge_preserves_mutation_count(
      n_children in 3_usize..8,
      n_shared in 1_usize..6,
      unique_counts in prop::collection::vec(0_usize..6, 3..8),
      length in 100_usize..500,
    ) {
      let unique_counts = &unique_counts[..n_children.min(unique_counts.len())];
      let n_children = unique_counts.len();
      if n_children < 3 { return Ok(()); }

      let (mut graph, names, edge_mutations, length, mut branch_lengths) =
        helpers::build_polytomy(n_children, n_shared, unique_counts, length);

      let original_total: usize = edge_mutations.iter().map(|(_, _, subs)| subs.len()).sum();
      let partition = helpers::make_partition(&graph, &names, length, &edge_mutations);
      let mut partitions = vec![partition];

      merge_shared_mutation_branches(&mut graph, &mut partitions, &mut branch_lengths).unwrap();
      graph.build().unwrap();

      let p = &partitions[0];
      let total_after: usize = graph
        .get_edges()
        .filter_map(|e| p.obs_edges.get(&e.key()))
        .map(|e| e.fitch_subs().len())
        .sum();

      prop_assert!(
        total_after <= original_total,
        "mutations created: before={original_total} after={total_after}"
      );
    }

    #[test]
    fn test_prop_merge_branch_lengths_non_negative(
      n_children in 3_usize..8,
      n_shared in 1_usize..6,
      unique_counts in prop::collection::vec(0_usize..6, 3..8),
      length in 100_usize..500,
    ) {
      let unique_counts = &unique_counts[..n_children.min(unique_counts.len())];
      let n_children = unique_counts.len();
      if n_children < 3 { return Ok(()); }

      let (mut graph, names, edge_mutations, length, mut branch_lengths) =
        helpers::build_polytomy(n_children, n_shared, unique_counts, length);
      let partition = helpers::make_partition(&graph, &names, length, &edge_mutations);
      let mut partitions = vec![partition];

      merge_shared_mutation_branches(&mut graph, &mut partitions, &mut branch_lengths).unwrap();
      graph.build().unwrap();

      for (edge_key, bl) in &branch_lengths {
        if let Some(bl) = bl {
          prop_assert!(*bl >= 0.0, "negative branch length {bl} on edge {edge_key}");
        }
      }
    }

    #[test]
    fn test_prop_merge_idempotent(
      n_children in 3_usize..6,
      n_shared in 1_usize..4,
      unique_counts in prop::collection::vec(0_usize..4, 3..6),
      length in 100_usize..300,
    ) {
      let unique_counts = &unique_counts[..n_children.min(unique_counts.len())];
      let n_children = unique_counts.len();
      if n_children < 3 { return Ok(()); }

      let (mut graph, names, edge_mutations, length, mut branch_lengths) =
        helpers::build_polytomy(n_children, n_shared, unique_counts, length);
      let partition = helpers::make_partition(&graph, &names, length, &edge_mutations);
      let mut partitions = vec![partition];

      let merged_first = merge_shared_mutation_branches(&mut graph, &mut partitions, &mut branch_lengths).unwrap();
      graph.build().unwrap();

      let merged_second = merge_shared_mutation_branches(&mut graph, &mut partitions, &mut branch_lengths).unwrap();

      prop_assert_eq!(merged_second, 0);
    }

    #[test]
    fn test_prop_merge_preserves_leaves(n_children in 3_usize..8,
      n_shared in 1_usize..6,
      unique_counts in prop::collection::vec(0_usize..6, 3..8),
      length in 100_usize..500,
    ) {
      let unique_counts = &unique_counts[..n_children.min(unique_counts.len())];
      let n_children = unique_counts.len();
      if n_children < 3 { return Ok(()); }

      let (mut graph, names, edge_mutations, length, mut branch_lengths) =
        helpers::build_polytomy(n_children, n_shared, unique_counts, length);

      let leaves_before: BTreeSet<String> = graph
        .get_leaves()
        .filter_map(|n| names.get(&n.key()).cloned().flatten())
        .collect();

      let partition = helpers::make_partition(&graph, &names, length, &edge_mutations);
      let mut partitions = vec![partition];

      merge_shared_mutation_branches(&mut graph, &mut partitions, &mut branch_lengths).unwrap();
      graph.build().unwrap();

      let leaves_after: BTreeSet<String> = graph
        .get_leaves()
        .filter_map(|n| names.get(&n.key()).cloned().flatten())
        .collect();

      prop_assert_eq!(leaves_before, leaves_after);
    }

    #[test]
    fn test_prop_merge_no_remaining_shared_mutations(
      n_children in 3_usize..6,
      n_shared in 1_usize..4,
      unique_counts in prop::collection::vec(0_usize..4, 3..6),
      length in 100_usize..300,
    ) {
      let unique_counts = &unique_counts[..n_children.min(unique_counts.len())];
      let n_children = unique_counts.len();
      if n_children < 3 { return Ok(()); }

      let (mut graph, names, edge_mutations, length, mut branch_lengths) =
        helpers::build_polytomy(n_children, n_shared, unique_counts, length);
      let partition = helpers::make_partition(&graph, &names, length, &edge_mutations);
      let mut partitions = vec![partition];

      merge_shared_mutation_branches(&mut graph, &mut partitions, &mut branch_lengths).unwrap();
      graph.build().unwrap();

      let p = &partitions[0];
      for node_ref in graph.get_nodes() {
        let node = node_ref;
        let outbound = node.outbound();
        if outbound.len() <= 1 { continue; }

        let child_sub_sets: Vec<BTreeSet<Sub>> = outbound
          .iter()
          .filter_map(|ek| p.obs_edges.get(ek))
          .map(|e| e.fitch_subs().iter().cloned().collect::<BTreeSet<_>>())
          .collect();

        for i in 0..child_sub_sets.len() {
          for j in (i + 1)..child_sub_sets.len() {
            let shared: Vec<_> = child_sub_sets[i].intersection(&child_sub_sets[j]).collect();
            prop_assert!(
              shared.is_empty(),
              "siblings under node {} still share {} mutations",
              node.key(),
              shared.len(),
            );
          }
        }
      }
    }
  }

  // === Boundary unit tests ===

  #[test]
  fn test_merge_all_children_share_same_mutation() -> Result<(), Report> {
    // Every child shares the same mutation. One group = all children.
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.1,C:0.1,D:0.1,E:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let partition = helpers::make_partition_from_static(
      &graph,
      &names,
      100,
      &[
        ("root", "A", vec![sub_at(0)]),
        ("root", "B", vec![sub_at(0)]),
        ("root", "C", vec![sub_at(0)]),
        ("root", "D", vec![sub_at(0)]),
        ("root", "E", vec![sub_at(0)]),
      ],
    )?;
    let mut partitions = vec![partition];

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 1);
    graph.build()?;

    // Root has 1 child (the new node), which has all 5 original children.
    let root = graph.get_roots().next().expect("root");
    assert_eq!(root.degree_out(), 1);

    let internal_edge = root.outbound()[0];
    let internal_key = graph.get_target_node_key(internal_edge)?;
    let internal = graph.get_node(internal_key).expect("internal node");
    assert_eq!(internal.degree_out(), 5);

    Ok(())
  }

  #[test]
  fn test_merge_overlapping_groups_greedy_selection() -> Result<(), Report> {
    // A,B share {sub0}. B,C share {sub1}. B in both groups.
    // Greedy picks one, second group excluded for this round.
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.1,C:0.1,D:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let partition = helpers::make_partition_from_static(
      &graph,
      &names,
      100,
      &[
        ("root", "A", vec![sub_at(0)]),
        ("root", "B", vec![sub_at(0), sub_at(1)]),
        ("root", "C", vec![sub_at(1)]),
        ("root", "D", vec![sub_at(2)]),
      ],
    )?;
    let mut partitions = vec![partition];

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut partitions, &mut branch_lengths)?;
    assert!(merged >= 1);

    Ok(())
  }

  #[test]
  fn test_merge_polytomy_reduced_to_binary_stops() -> Result<(), Report> {
    // 3 children, 2 share. After merge: binary tree, loop stops.
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.1,C:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let partition = helpers::make_partition_from_static(
      &graph,
      &names,
      100,
      &[
        ("root", "A", vec![sub_at(0)]),
        ("root", "B", vec![sub_at(0)]),
        ("root", "C", vec![sub_at(1)]),
      ],
    )?;
    let mut partitions = vec![partition];

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 1);
    graph.build()?;

    let root = graph.get_roots().next().expect("root");
    assert_eq!(root.degree_out(), 2);

    Ok(())
  }

  #[test]
  fn test_merge_multiple_polytomies_in_one_tree() -> Result<(), Report> {
    // Two independent polytomies: root has {I, D, E, F}, I has {A, B, C}.
    // A,B share sub0 under I. D,E share sub1 under root.
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1,C:0.1)I:0.1,D:0.1,E:0.1,F:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let partition = helpers::make_partition_from_static(
      &graph,
      &names,
      100,
      &[
        ("root", "I", vec![]),
        ("I", "A", vec![sub_at(0)]),
        ("I", "B", vec![sub_at(0)]),
        ("I", "C", vec![sub_at(5)]),
        ("root", "D", vec![sub_at(1)]),
        ("root", "E", vec![sub_at(1)]),
        ("root", "F", vec![sub_at(6)]),
      ],
    )?;
    let mut partitions = vec![partition];

    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 2);

    Ok(())
  }

  #[test]
  fn test_merge_disjoint_sub_and_indel_groups_same_round() -> Result<(), Report> {
    // A,B share a sub. C,D share an indel. No overlap between groups.
    // Both groups should merge (possibly in one round).
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.1,C:0.1,D:0.1,E:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let mut partition = helpers::make_partition_from_static(
      &graph,
      &names,
      100,
      &[
        ("root", "A", vec![sub_at(0)]),
        ("root", "B", vec![sub_at(0)]),
        ("root", "C", vec![]),
        ("root", "D", vec![]),
        ("root", "E", vec![sub_at(5)]),
      ],
    )?;

    let shared_indel = InDel::del((10, 13), Seq::try_from_str("GTA").unwrap()).unwrap();
    {
      let p = &mut partition;
      let edge_c = find_edge_key(&graph, &names, "root", "C").expect("edge root->C");
      let edge_d = find_edge_key(&graph, &names, "root", "D").expect("edge root->D");
      p.obs_edges.get_mut(&edge_c).expect("partition edge C").indels = vec![shared_indel.clone()];
      p.obs_edges.get_mut(&edge_d).expect("partition edge D").indels = vec![shared_indel];
    }

    let mut partitions = vec![partition];
    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut partitions, &mut branch_lengths)?;
    assert_eq!(merged, 2);

    Ok(())
  }

  #[test]
  fn test_merge_multi_partition_asymmetric_sharing() -> Result<(), Report> {
    // Partition 1: A,B share sub at pos 0. Partition 2: A,C share sub at pos 50.
    // Total shared(A,B) = 1 (from p1). Total shared(A,C) = 1 (from p2).
    // Both groups have equal score. Greedy picks one.
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.1,C:0.1,D:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let p1 = helpers::make_partition_from_static(
      &graph,
      &names,
      100,
      &[
        ("root", "A", vec![sub_at(0)]),
        ("root", "B", vec![sub_at(0)]),
        ("root", "C", vec![sub_at(5)]),
        ("root", "D", vec![sub_at(6)]),
      ],
    )?;

    let p2 = helpers::make_partition_from_static_indexed(
      &graph,
      &names,
      1,
      200,
      &[
        ("root", "A", vec![Sub::new(c(b'C'), 50_usize, c(b'G')).unwrap()]),
        ("root", "B", vec![Sub::new(c(b'G'), 60_usize, c(b'T')).unwrap()]),
        ("root", "C", vec![Sub::new(c(b'C'), 50_usize, c(b'G')).unwrap()]),
        ("root", "D", vec![Sub::new(c(b'T'), 70_usize, c(b'A')).unwrap()]),
      ],
    )?;

    let mut partitions = vec![p1, p2];
    let mut branch_lengths = branch_lengths;
    let merged = merge_shared_mutation_branches(&mut graph, &mut partitions, &mut branch_lengths)?;
    assert!(merged >= 1);
    graph.build()?;

    // A must be in exactly one merged group, not both.
    let leaf_count = graph
      .get_leaves()
      .filter(|n| names.get(&n.key()).and_then(|x| x.as_ref()).is_some())
      .count();
    assert_eq!(leaf_count, 4);

    Ok(())
  }

  mod helpers {
    use super::*;

    type PolytomySetup = (
      Graph,
      BTreeMap<GraphNodeKey, Option<String>>,
      Vec<(String, String, Vec<Sub>)>,
      usize,
      BTreeMap<GraphEdgeKey, Option<f64>>,
    );

    pub fn build_polytomy(
      n_children: usize,
      n_shared: usize,
      unique_counts: &[usize],
      min_length: usize,
    ) -> PolytomySetup {
      let names: Vec<String> = (0..n_children).map(|i| format!("N{i}")).collect();
      let newick_children = names.iter().map(|n| format!("{n}:0.1")).join(",");
      let newick = format!("({newick_children})root;");
      let nwk_parsed = nwk_read_str(&newick).unwrap();
      let node_names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;

      let mut pos_counter = 0_usize;
      let shared: Vec<Sub> = (pos_counter..pos_counter + n_shared).map(sub_at).collect();
      pos_counter += n_shared;

      let mut edge_mutations: Vec<(String, String, Vec<Sub>)> = Vec::new();
      for (i, name) in names.iter().enumerate() {
        let n_unique = unique_counts[i];
        let mut subs = shared.clone();
        for _ in 0..n_unique {
          subs.push(sub_at(pos_counter));
          pos_counter += 1;
        }
        subs.sort();
        edge_mutations.push(("root".to_owned(), name.clone(), subs));
      }

      let length = min_length.max(pos_counter + 1);
      (graph, node_names, edge_mutations, length, branch_lengths)
    }

    pub fn make_partition(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      length: usize,
      edge_mutations: &[(String, String, Vec<Sub>)],
    ) -> PartitionMarginalSparse {
      let refs: Vec<(&str, &str, Vec<Sub>)> = edge_mutations
        .iter()
        .map(|(src, tgt, subs)| (src.as_str(), tgt.as_str(), subs.clone()))
        .collect();
      make_partition_from_static(graph, names, length, &refs).unwrap()
    }

    pub fn make_partition_from_static(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      length: usize,
      edge_mutations: &[(&str, &str, Vec<Sub>)],
    ) -> Result<PartitionMarginalSparse, Report> {
      make_partition_from_static_indexed(graph, names, 0, length, edge_mutations)
    }

    pub fn make_partition_from_static_indexed(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      index: usize,
      length: usize,
      edge_mutations: &[(&str, &str, Vec<Sub>)],
    ) -> Result<PartitionMarginalSparse, Report> {
      let alphabet = Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?;

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
          find_edge_key(graph, names, source, target).unwrap_or_else(|| panic!("edge {source}->{target} not found"));
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

      Ok(partition)
    }
  }
}

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;

  use crate::gtr::jc_distance::jukes_cantor_distance;
  use crate::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::payload::ancestral::GraphAncestral;
  use crate::partition::sparse::{SparseEdgePartition, SparseNodePartition};
  use crate::seq::indel::InDel;
  use crate::seq::mutation::Sub;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use approx::assert_relative_eq;
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::seq;
  use treetime_primitives::{AsciiChar, Seq};

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
    Sub::new(c(reff), pos, c(qry)).unwrap()
  }

  fn make_partition(
    graph: &GraphAncestral,
    length: usize,
    edge_mutations: &[(&str, &str, Vec<Sub>)],
  ) -> Result<Arc<RwLock<PartitionMarginalSparse>>, Report> {
    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length,
      nodes: btreemap! {},
      edges: btreemap! {},
      root_sequence: seq![],
    };

    // Build root reference sequence consistent with edge subs.
    // Set each position to the sub's ref character so that edge_subs()
    // produces the same mutations as the stored subs.
    let mut ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
    for (_, _, subs) in edge_mutations {
      for s in subs {
        if s.pos() < length {
          ref_seq[s.pos()] = s.reff();
        }
      }
    }

    partition.root_sequence = ref_seq.clone();

    // Populate node entries so edge_subs() can reconstruct states
    for node in graph.get_nodes() {
      let key = node.read_arc().key();
      let mut node_part = SparseNodePartition::empty(&partition.alphabet);
      node_part.seq.sequence = ref_seq.clone();
      partition.nodes.insert(key, node_part);
    }

    for (source, target, subs) in edge_mutations {
      let edge_key =
        find_edge_key(graph, source, target).unwrap_or_else(|| panic!("edge {source}->{target} not found in graph"));
      partition
        .edges
        .insert(edge_key, SparseEdgePartition::with_fitch_subs(subs.clone()));
    }

    Ok(Arc::new(RwLock::new(partition)))
  }

  /// Find the new internal node (unnamed, non-root, non-leaf).
  fn find_unnamed_internal_nodes(graph: &GraphAncestral) -> Vec<GraphNodeKey> {
    graph
      .get_nodes()
      .iter()
      .filter_map(|node| {
        let node = node.read_arc();
        let payload = node.payload().read_arc();
        let is_unnamed = payload.name.is_none();
        let is_internal = !node.is_leaf() && !node.is_root();
        (is_unnamed && is_internal).then_some(node.key())
      })
      .collect()
  }

  #[test]
  fn test_merge_no_polytomy() -> Result<(), Report> {
    // Binary tree: no polytomy, nothing to merge
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)internal:0.3)root;")?;
    let partition = make_partition(
      &graph,
      100,
      &[
        ("root", "internal", vec![sub(b'A', 0, b'T')]),
        ("internal", "A", vec![sub(b'A', 0, b'T')]),
        ("internal", "B", vec![sub(b'A', 0, b'T')]),
      ],
    )?;
    let partitions = vec![partition];

    let merged = merge_shared_mutation_branches(&mut graph, &partitions)?;
    assert_eq!(merged, 0);
    assert_eq!(graph.get_nodes().len(), 4); // root, internal, A, B
    Ok(())
  }

  #[test]
  fn test_merge_polytomy_no_shared_mutations() -> Result<(), Report> {
    // Polytomy with 3 children, no shared mutations
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2,C:0.3)root;")?;
    let partition = make_partition(
      &graph,
      100,
      &[
        ("root", "A", vec![sub(b'A', 0, b'T')]),
        ("root", "B", vec![sub(b'G', 5, b'C')]),
        ("root", "C", vec![sub(b'T', 10, b'A')]),
      ],
    )?;
    let partitions = vec![partition];

    let merged = merge_shared_mutation_branches(&mut graph, &partitions)?;
    assert_eq!(merged, 0);
    assert_eq!(graph.get_nodes().len(), 4); // root, A, B, C
    Ok(())
  }

  #[test]
  fn test_merge_polytomy_two_siblings_share_all_mutations() -> Result<(), Report> {
    // Polytomy: A and B share identical mutations, C is different
    //   root -> A (subs: A0T, G5C)
    //   root -> B (subs: A0T, G5C)
    //   root -> C (subs: T10A)
    // After merge: root -> N (subs: A0T, G5C), N -> A (subs: []), N -> B (subs: [])
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2,C:0.3)root;")?;
    let shared = vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')];
    let partition = make_partition(
      &graph,
      100,
      &[
        ("root", "A", shared.clone()),
        ("root", "B", shared),
        ("root", "C", vec![sub(b'T', 10, b'A')]),
      ],
    )?;
    let partitions = vec![partition];

    let merged = merge_shared_mutation_branches(&mut graph, &partitions)?;
    assert_eq!(merged, 1);

    graph.build()?;
    // New topology: root -> N, N -> A, N -> B, root -> C
    assert_eq!(graph.get_nodes().len(), 5); // root, N, A, B, C
    assert_eq!(graph.get_edges().len(), 4);

    // The new internal node should exist
    let unnamed = find_unnamed_internal_nodes(&graph);
    assert_eq!(unnamed.len(), 1);

    // Check partition data: new edge to N has shared mutations
    let p = partitions[0].read_arc();
    // A and B should have no remaining mutations
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      if let Some(edge_data) = p.edges.get(&edge.key()) {
        match target_name.as_deref() {
          Some("A" | "B") => assert_eq!(
            edge_data.fitch_subs().len(),
            0,
            "child should have no remaining mutations"
          ),
          Some("C") => assert_eq!(edge_data.fitch_subs().len(), 1),
          None => assert_eq!(
            edge_data.fitch_subs().len(),
            2,
            "new internal edge should carry shared mutations"
          ),
          _ => {},
        }
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_polytomy_partial_overlap() -> Result<(), Report> {
    // A has {A0T, G5C}, B has {A0T, G5C, T10A}, C has {C20G}
    // Shared(A,B) = {A0T, G5C} (2), so merge A and B
    // After: N -> A (subs: []), N -> B (subs: {T10A})
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2,C:0.3)root;")?;
    let partition = make_partition(
      &graph,
      100,
      &[
        ("root", "A", vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')]),
        (
          "root",
          "B",
          vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C'), sub(b'T', 10, b'A')],
        ),
        ("root", "C", vec![sub(b'C', 20, b'G')]),
      ],
    )?;
    let partitions = vec![partition];

    let merged = merge_shared_mutation_branches(&mut graph, &partitions)?;
    assert_eq!(merged, 1);

    graph.build()?;
    let p = partitions[0].read_arc();

    // Check that B retains only its unique mutation
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      if let Some(edge_data) = p.edges.get(&edge.key()) {
        match target_name.as_deref() {
          Some("A") => assert_eq!(edge_data.fitch_subs().len(), 0),
          Some("B") => {
            assert_eq!(edge_data.fitch_subs().len(), 1);
            assert_eq!(edge_data.fitch_subs()[0], sub(b'T', 10, b'A'));
          },
          Some("C") => assert_eq!(edge_data.fitch_subs().len(), 1),
          None => assert_eq!(
            edge_data.fitch_subs().len(),
            2,
            "internal edge carries shared mutations"
          ),
          _ => {},
        }
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_greedy_picks_best_pair() -> Result<(), Report> {
    // A has {A0T, G5C, T10A}, B has {A0T, G5C, T10A}, C has {C20G}, D has {T30A}
    // Shared(A,B) = 3, no other pair shares mutations
    // Greedy merges A and B
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2,C:0.3,D:0.4)root;")?;
    let partition = make_partition(
      &graph,
      100,
      &[
        (
          "root",
          "A",
          vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C'), sub(b'T', 10, b'A')],
        ),
        (
          "root",
          "B",
          vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C'), sub(b'T', 10, b'A')],
        ),
        ("root", "C", vec![sub(b'C', 20, b'G')]),
        ("root", "D", vec![sub(b'T', 30, b'A')]),
      ],
    )?;
    let partitions = vec![partition];

    let merged = merge_shared_mutation_branches(&mut graph, &partitions)?;
    assert_eq!(merged, 1);

    graph.build()?;
    // Topology: root -> N -> {A, B}, root -> C, root -> D
    assert_eq!(graph.get_nodes().len(), 6); // root, N, A, B, C, D

    let p = partitions[0].read_arc();
    // The new internal edge carries all 3 shared mutations
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      if target_name.is_none() {
        if let Some(edge_data) = p.edges.get(&edge.key()) {
          assert_eq!(edge_data.fitch_subs().len(), 3);
        }
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_branch_length_adjustment() -> Result<(), Report> {
    // A and B share 2 mutations out of length 100. Both have 0 unique mutations.
    // Parent edge: jc(2/100).
    // Child edges: jc(0/100) = 0.0 (no remaining mutations).
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2,C:0.3)root;")?;
    let shared = vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')];
    let partition = make_partition(
      &graph,
      100,
      &[
        ("root", "A", shared.clone()),
        ("root", "B", shared),
        ("root", "C", vec![sub(b'T', 10, b'A')]),
      ],
    )?;
    let partitions = vec![partition];

    merge_shared_mutation_branches(&mut graph, &partitions)?;
    graph.build()?;

    let d = -0.75 * f64::ln(1.0 - 4.0 * 0.02 / 3.0);
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      let bl = edge.payload().read_arc().branch_length;

      match target_name.as_deref() {
        None => assert_relative_eq!(bl.unwrap(), d, epsilon = 1e-15),
        Some("A") => assert_relative_eq!(bl.unwrap(), 0.0, epsilon = 1e-15),
        Some("B") => assert_relative_eq!(bl.unwrap(), 0.0, epsilon = 1e-15),
        Some("C") => assert_relative_eq!(bl.unwrap(), 0.3, epsilon = 1e-6),
        _ => {},
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_child_bl_zero_when_all_shared() -> Result<(), Report> {
    // 10 shared mutations out of length 100. Both children have only shared
    // mutations (0 remaining), so child BLs are jc(0) = 0.0.
    let mut graph: GraphAncestral = nwk_read_str("(A:0.05,B:0.2,C:0.3)root;")?;
    let shared: Vec<Sub> = (0..10).map(|i| sub(b'A', i, b'T')).collect();
    let partition = make_partition(
      &graph,
      100,
      &[
        ("root", "A", shared.clone()),
        ("root", "B", shared),
        ("root", "C", vec![sub(b'T', 50, b'A')]),
      ],
    )?;
    let partitions = vec![partition];

    merge_shared_mutation_branches(&mut graph, &partitions)?;
    graph.build()?;

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      let bl = edge.payload().read_arc().branch_length;

      if target_name.as_deref() == Some("A") {
        assert_relative_eq!(bl.unwrap(), 0.0, epsilon = 1e-15);
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_multiple_partitions() -> Result<(), Report> {
    // Two partitions: shared in both
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2,C:0.3)root;")?;

    let p1 = make_partition(
      &graph,
      100,
      &[
        ("root", "A", vec![sub(b'A', 0, b'T')]),
        ("root", "B", vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')]),
        ("root", "C", vec![sub(b'T', 10, b'A')]),
      ],
    )?;

    let mut p2_inner = PartitionMarginalSparse {
      index: 1,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length: 200,
      nodes: btreemap! {},
      edges: btreemap! {},
      root_sequence: seq![],
    };
    let mut p2_ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(200).collect();
    p2_ref_seq[50] = c(b'C'); // sub C50G uses ref='C'
    p2_inner.root_sequence = p2_ref_seq.clone();
    for node in graph.get_nodes() {
      let key = node.read_arc().key();
      let mut node_part = SparseNodePartition::empty(&p2_inner.alphabet);
      node_part.seq.sequence = p2_ref_seq.clone();
      p2_inner.nodes.insert(key, node_part);
    }

    // Partition 2: A and B share mutation at pos 50
    let edge_a = find_edge_key(&graph, "root", "A").unwrap();
    let edge_b = find_edge_key(&graph, "root", "B").unwrap();
    let edge_c = find_edge_key(&graph, "root", "C").unwrap();
    p2_inner
      .edges
      .insert(edge_a, SparseEdgePartition::with_fitch_subs(vec![sub(b'C', 50, b'G')]));
    p2_inner
      .edges
      .insert(edge_b, SparseEdgePartition::with_fitch_subs(vec![sub(b'C', 50, b'G')]));
    p2_inner.edges.insert(edge_c, SparseEdgePartition::default());
    let p2 = Arc::new(RwLock::new(p2_inner));

    let partitions = vec![p1, p2];

    let merged = merge_shared_mutation_branches(&mut graph, &partitions)?;
    assert_eq!(merged, 1);
    graph.build()?;

    // Total shared = 1 (from p1) + 1 (from p2) = 2 across total length 100 + 200 = 300
    // Pooled p-distance = 2/300 ≈ 0.006667
    // JC69 correction with k=4: d = -3/4 * ln(1 - 4*p/3)
    let p_pooled = 2.0 / 300.0;
    let d_expected = -0.75 * f64::ln(1.0 - 4.0 * p_pooled / 3.0);
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      let bl = edge.payload().read_arc().branch_length;

      if target_name.is_none() {
        assert_relative_eq!(bl.unwrap(), d_expected, epsilon = 1e-15);
      }
    }

    // Verify partition 1: internal edge has 1 shared sub, B has 1 remaining
    let p1 = partitions[0].read_arc();
    let p2 = partitions[1].read_arc();
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      if target_name.as_deref() == Some("B") {
        // P1: B had {A0T, G5C}, shared = {A0T}, remaining = {G5C}
        assert_eq!(p1.edges[&edge.key()].fitch_subs().len(), 1);
        assert_eq!(p1.edges[&edge.key()].fitch_subs()[0], sub(b'G', 5, b'C'));
        // P2: B had {C50G}, shared = {C50G}, remaining = {}
        assert_eq!(p2.edges[&edge.key()].fitch_subs().len(), 0);
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_repeated_until_exhausted() -> Result<(), Report> {
    // 5 children: A, B share {A0T}; C, D share {G5C}; E has nothing shared
    // First merge: pick pair with more shared (both have 1, so either pair)
    // After first merge: remaining polytomy of 4 children (N1, <one of C/D/E>, <other>, <other>)
    // If C and D still share mutations, they get merged too
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.1,C:0.1,D:0.1,E:0.1)root;")?;
    let partition = make_partition(
      &graph,
      100,
      &[
        ("root", "A", vec![sub(b'A', 0, b'T')]),
        ("root", "B", vec![sub(b'A', 0, b'T')]),
        ("root", "C", vec![sub(b'G', 5, b'C')]),
        ("root", "D", vec![sub(b'G', 5, b'C')]),
        ("root", "E", vec![sub(b'T', 10, b'A')]),
      ],
    )?;
    let partitions = vec![partition];

    let merged = merge_shared_mutation_branches(&mut graph, &partitions)?;
    // Two merge rounds: A+B, then C+D
    assert_eq!(merged, 2);

    graph.build()?;
    // root -> N1 -> {A, B}, root -> N2 -> {C, D}, root -> E
    assert_eq!(graph.get_nodes().len(), 8); // root, N1, N2, A, B, C, D, E
    assert_eq!(graph.get_edges().len(), 7);

    Ok(())
  }

  #[test]
  fn test_merge_empty_partitions() -> Result<(), Report> {
    // No partition data at all
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2,C:0.3)root;")?;
    let partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];

    let merged = merge_shared_mutation_branches(&mut graph, &partitions)?;
    assert_eq!(merged, 0);
    Ok(())
  }

  #[test]
  fn test_merge_preserves_tree_structure_for_non_polytomies() -> Result<(), Report> {
    // Mix of polytomy and binary nodes
    // root -> internal1 -> {A, B, C} (polytomy), root -> D
    // A and B share mutations, C is different
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1,C:0.1)internal1:0.1,D:0.2)root;")?;
    let partition = make_partition(
      &graph,
      100,
      &[
        ("root", "internal1", vec![sub(b'C', 20, b'G')]),
        ("internal1", "A", vec![sub(b'A', 0, b'T')]),
        ("internal1", "B", vec![sub(b'A', 0, b'T')]),
        ("internal1", "C", vec![sub(b'G', 5, b'C')]),
        ("root", "D", vec![sub(b'T', 10, b'A')]),
      ],
    )?;
    let partitions = vec![partition];

    let merged = merge_shared_mutation_branches(&mut graph, &partitions)?;
    assert_eq!(merged, 1);

    graph.build()?;
    // internal1 -> N -> {A, B}, internal1 -> C
    assert_eq!(graph.get_nodes().len(), 7); // root, internal1, N, A, B, C, D

    // D should be unaffected
    assert!(find_node_key_by_name(&graph, "D").is_some());
    assert!(find_edge_key(&graph, "root", "D").is_some());

    Ok(())
  }

  #[test]
  fn test_merge_single_mutation_shared() -> Result<(), Report> {
    // Minimal case: exactly 1 shared mutation
    let mut graph: GraphAncestral = nwk_read_str("(A:0.05,B:0.05,C:0.05)root;")?;
    let partition = make_partition(
      &graph,
      1000,
      &[
        ("root", "A", vec![sub(b'A', 42, b'T')]),
        ("root", "B", vec![sub(b'A', 42, b'T')]),
        ("root", "C", vec![]),
      ],
    )?;
    let partitions = vec![partition];

    let merged = merge_shared_mutation_branches(&mut graph, &partitions)?;
    assert_eq!(merged, 1);
    graph.build()?;

    // 1 shared mutation out of 1000 => p = 0.001
    // JC69 correction: d = -3/4 * ln(1 - 4*0.001/3) ≈ 0.001001
    // At small p the correction is a few parts per thousand above raw p.
    let d_expected = -0.75 * f64::ln(1.0 - 4.0 * 0.001 / 3.0);
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      let bl = edge.payload().read_arc().branch_length;

      if target_name.is_none() {
        assert_relative_eq!(bl.unwrap(), d_expected, epsilon = 1e-15);
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_branch_length_jc_correction_differs_from_raw() -> Result<(), Report> {
    // 10 shared mutations out of length 100 places the pooled p-distance at
    // 0.10, where JC69 correction differs from the raw ratio by about 7%.
    // Both children have 0 unique mutations, so child BLs are 0.0.
    let mut graph: GraphAncestral = nwk_read_str("(A:0.5,B:0.5,C:0.5)root;")?;
    let shared: Vec<Sub> = (0..10).map(|i| sub(b'A', i, b'T')).collect();
    let partition = make_partition(
      &graph,
      100,
      &[
        ("root", "A", shared.clone()),
        ("root", "B", shared),
        ("root", "C", vec![sub(b'T', 50, b'A')]),
      ],
    )?;
    let partitions = vec![partition];

    merge_shared_mutation_branches(&mut graph, &partitions)?;
    graph.build()?;

    let p = 0.10;
    let d = -0.75 * f64::ln(1.0 - 4.0 * p / 3.0);
    assert!(d > p * 1.05, "JC correction must exceed raw p by >5%: d={d} p={p}");
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      let bl = edge.payload().read_arc().branch_length;

      match target_name.as_deref() {
        None => assert_relative_eq!(bl.unwrap(), d, epsilon = 1e-15),
        Some("A" | "B") => assert_relative_eq!(bl.unwrap(), 0.0, epsilon = 1e-15),
        _ => {},
      }
    }

    Ok(())
  }

  // === Child branch length = JC69(remaining_mutations / alignment_length) ===

  #[rustfmt::skip]
  #[rstest]
  #[case::basic_remaining(       (2, 0, 1),  100, 0.2)]
  #[case::zero_remaining(        (2, 0, 0),  100, 0.2)]
  #[case::asymmetric_remaining(  (3, 2, 5), 1000, 0.5)]
  #[case::newick_bl_independent( (2, 0, 1),  100, 99.0)]
  #[trace]
  fn test_merge_child_bl(
    #[case] (n_shared, n_unique_a, n_unique_b): (usize, usize, usize),
    #[case] length: usize,
    #[case] newick_bl: f64,
  ) -> Result<(), Report> {
    let newick = format!("(A:{newick_bl},B:{newick_bl},C:{newick_bl})root;");
    let mut graph: GraphAncestral = nwk_read_str(&newick)?;

    let edge_subs = helpers::build_shared_unique_subs(n_shared, n_unique_a, n_unique_b);
    let partition = make_partition(&graph, length, &edge_subs)?;
    let partitions = vec![partition];

    merge_shared_mutation_branches(&mut graph, &partitions)?;
    graph.build()?;

    let bls = helpers::extract_branch_lengths(&graph);
    let expected_bl_a = helpers::expected_jc_bl(n_unique_a, length);
    let expected_bl_b = helpers::expected_jc_bl(n_unique_b, length);
    assert_relative_eq!(bls["A"], expected_bl_a, epsilon = 1e-15);
    assert_relative_eq!(bls["B"], expected_bl_b, epsilon = 1e-15);

    Ok(())
  }

  #[test]
  fn test_merge_child_bl_includes_indels_in_remaining() -> Result<(), Report> {
    // A has 2 shared subs + 1 indel on its edge (counted in remaining).
    // B has 2 shared subs only. bl_a = jc(1/100), bl_b = 0.0.
    let mut graph: GraphAncestral = nwk_read_str("(A:0.5,B:0.5,C:0.5)root;")?;
    let shared = vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')];
    let partition = make_partition(
      &graph,
      100,
      &[
        ("root", "A", shared.clone()),
        ("root", "B", shared),
        ("root", "C", vec![sub(b'T', 20, b'A')]),
      ],
    )?;

    {
      let mut p = partition.write_arc();
      let edge_a = find_edge_key(&graph, "root", "A").expect("edge root->A");
      p.edges.get_mut(&edge_a).expect("partition edge A").indels =
        vec![InDel::del((10, 13), Seq::try_from_str("GTA").unwrap())];
    }

    let partitions = vec![partition];
    merge_shared_mutation_branches(&mut graph, &partitions)?;
    graph.build()?;

    let bls = helpers::extract_branch_lengths(&graph);
    assert_relative_eq!(bls["A"], helpers::expected_jc_bl(1, 100), epsilon = 1e-15);
    assert_relative_eq!(bls["B"], helpers::expected_jc_bl(0, 100), epsilon = 1e-15);

    Ok(())
  }

  #[test]
  fn test_merge_child_bl_across_partitions() -> Result<(), Report> {
    // Two partitions (lengths 100 and 200). A has 1 unique sub in p1, 2 in p2.
    // Total remaining = 3, total length = 300. bl_a = jc(3/300).
    let mut graph: GraphAncestral = nwk_read_str("(A:0.5,B:0.5,C:0.5)root;")?;

    let p1 = make_partition(
      &graph,
      100,
      &[
        (
          "root",
          "A",
          vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C'), sub(b'T', 10, b'A')],
        ),
        ("root", "B", vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')]),
        ("root", "C", vec![sub(b'T', 20, b'A')]),
      ],
    )?;

    let p2 = helpers::make_second_partition(
      &graph,
      200,
      &[
        (
          ("root", "A"),
          vec![sub(b'C', 50, b'G'), sub(b'G', 60, b'T'), sub(b'T', 70, b'A')],
        ),
        (("root", "B"), vec![sub(b'C', 50, b'G')]),
        (("root", "C"), vec![]),
      ],
    )?;

    let partitions = vec![p1, p2];
    merge_shared_mutation_branches(&mut graph, &partitions)?;
    graph.build()?;

    let bls = helpers::extract_branch_lengths(&graph);
    assert_relative_eq!(bls["A"], helpers::expected_jc_bl(3, 300), epsilon = 1e-15);
    assert_relative_eq!(bls["B"], helpers::expected_jc_bl(0, 300), epsilon = 1e-15);

    Ok(())
  }

  // === Group merge (k >= 2 siblings) and indel sharing ===

  #[test]
  fn test_merge_group_three_siblings_same_mutation() -> Result<(), Report> {
    // A, B, C all share {A0T}. All 3 should land under one new internal node.
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.1,C:0.1,D:0.1)root;")?;
    let partition = make_partition(
      &graph,
      100,
      &[
        ("root", "A", vec![sub(b'A', 0, b'T')]),
        ("root", "B", vec![sub(b'A', 0, b'T')]),
        ("root", "C", vec![sub(b'A', 0, b'T')]),
        ("root", "D", vec![sub(b'G', 5, b'C')]),
      ],
    )?;
    let partitions = vec![partition];

    let merged = merge_shared_mutation_branches(&mut graph, &partitions)?;
    assert_eq!(merged, 1);
    graph.build()?;

    let unnamed = find_unnamed_internal_nodes(&graph);
    assert_eq!(unnamed.len(), 1);

    let new_node = graph.get_node(unnamed[0]).expect("new internal node");
    assert_eq!(new_node.read_arc().degree_out(), 3);

    Ok(())
  }

  #[test]
  fn test_merge_shared_indels_only() -> Result<(), Report> {
    // A and B share an identical deletion but no substitutions.
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.1,C:0.1)root;")?;
    let partition = make_partition(
      &graph,
      100,
      &[
        ("root", "A", vec![]),
        ("root", "B", vec![]),
        ("root", "C", vec![sub(b'T', 20, b'A')]),
      ],
    )?;

    let shared_indel = InDel::del((5, 8), Seq::try_from_str("GTA").unwrap());
    {
      let mut p = partition.write_arc();
      let edge_a = find_edge_key(&graph, "root", "A").expect("edge root->A");
      let edge_b = find_edge_key(&graph, "root", "B").expect("edge root->B");
      p.edges.get_mut(&edge_a).expect("partition edge A").indels = vec![shared_indel.clone()];
      p.edges.get_mut(&edge_b).expect("partition edge B").indels = vec![shared_indel];
    }

    let partitions = vec![partition];
    let merged = merge_shared_mutation_branches(&mut graph, &partitions)?;
    assert_eq!(merged, 1);

    Ok(())
  }

  #[test]
  fn test_merge_shared_subs_and_indels_split_correctly() -> Result<(), Report> {
    // A and B share 1 sub + 1 indel. A also has 1 unique sub.
    let mut graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.1,C:0.1)root;")?;
    let partition = make_partition(
      &graph,
      100,
      &[
        ("root", "A", vec![sub(b'A', 0, b'T'), sub(b'G', 15, b'C')]),
        ("root", "B", vec![sub(b'A', 0, b'T')]),
        ("root", "C", vec![sub(b'T', 20, b'A')]),
      ],
    )?;

    let shared_indel = InDel::del((5, 8), Seq::try_from_str("GTA").unwrap());
    {
      let mut p = partition.write_arc();
      let edge_a = find_edge_key(&graph, "root", "A").expect("edge root->A");
      let edge_b = find_edge_key(&graph, "root", "B").expect("edge root->B");
      p.edges.get_mut(&edge_a).expect("partition edge A").indels = vec![shared_indel.clone()];
      p.edges.get_mut(&edge_b).expect("partition edge B").indels = vec![shared_indel];
    }

    let partitions = vec![partition];
    merge_shared_mutation_branches(&mut graph, &partitions)?;
    graph.build()?;

    let edge_data = helpers::extract_edge_mutation_counts(&graph, &partitions[0]);
    assert_eq!(edge_data[&None], (1, 1), "parent edge: 1 shared sub, 1 shared indel");
    assert_eq!(edge_data[&Some("A")], (1, 0), "A: 1 unique sub, no indels");
    assert_eq!(edge_data[&Some("B")], (0, 0), "B: no remaining mutations");

    Ok(())
  }

  mod helpers {
    use super::*;

    pub fn expected_jc_bl(count: usize, length: usize) -> f64 {
      if length == 0 || count == 0 {
        return 0.0;
      }
      jukes_cantor_distance(count as f64 / length as f64, 4)
    }

    pub fn extract_branch_lengths(graph: &GraphAncestral) -> BTreeMap<String, f64> {
      graph
        .get_edges()
        .iter()
        .filter_map(|edge| {
          let edge = edge.read_arc();
          let target = graph.get_node(edge.target())?;
          let name = target.read_arc().payload().read_arc().name.clone()?;
          let bl = edge.payload().read_arc().branch_length?;
          Some((name, bl))
        })
        .collect()
    }

    pub fn extract_edge_mutation_counts<'a>(
      graph: &GraphAncestral,
      partition: &Arc<RwLock<PartitionMarginalSparse>>,
    ) -> BTreeMap<Option<&'a str>, (usize, usize)> {
      let p = partition.read_arc();
      let mut result = BTreeMap::new();
      for edge in graph.get_edges() {
        let edge = edge.read_arc();
        let target = graph.get_node(edge.target()).expect("target node");
        let target_name = target.read_arc().payload().read_arc().name.clone();
        if let Some(edge_data) = p.edges.get(&edge.key()) {
          let key: Option<&'a str> = match target_name.as_deref() {
            Some("A") => Some("A"),
            Some("B") => Some("B"),
            Some("C") => Some("C"),
            Some("D") => Some("D"),
            None => None,
            _ => continue,
          };
          result.insert(key, (edge_data.fitch_subs().len(), edge_data.indels.len()));
        }
      }
      result
    }

    pub fn build_shared_unique_subs<'a>(
      n_shared: usize,
      n_unique_a: usize,
      n_unique_b: usize,
    ) -> Vec<(&'a str, &'a str, Vec<Sub>)> {
      let nucs = [b'A', b'C', b'G', b'T'];
      let shared: Vec<Sub> = (0..n_shared).map(|i| sub(nucs[i % 4], i, nucs[(i + 1) % 4])).collect();

      let mut subs_a = shared.clone();
      for i in 0..n_unique_a {
        let pos = n_shared + i;
        subs_a.push(sub(nucs[pos % 4], pos, nucs[(pos + 1) % 4]));
      }
      subs_a.sort();

      let mut subs_b = shared;
      for i in 0..n_unique_b {
        let pos = n_shared + n_unique_a + i;
        subs_b.push(sub(nucs[pos % 4], pos, nucs[(pos + 1) % 4]));
      }
      subs_b.sort();

      let c_pos = n_shared + n_unique_a + n_unique_b;
      vec![
        ("root", "A", subs_a),
        ("root", "B", subs_b),
        ("root", "C", vec![sub(nucs[c_pos % 4], c_pos, nucs[(c_pos + 1) % 4])]),
      ]
    }

    pub fn make_second_partition(
      graph: &GraphAncestral,
      length: usize,
      edge_subs: &[((&str, &str), Vec<Sub>)],
    ) -> Result<Arc<RwLock<PartitionMarginalSparse>>, Report> {
      let mut partition = PartitionMarginalSparse {
        index: 1,
        gtr: jc69(JC69Params::default())?,
        alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
        length,
        nodes: btreemap! {},
        edges: btreemap! {},
        root_sequence: seq![],
      };

      let mut ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
      for (_, subs) in edge_subs {
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

      for ((source, target), subs) in edge_subs {
        let edge_key =
          find_edge_key(graph, source, target).unwrap_or_else(|| panic!("edge {source}->{target} not found in graph"));
        partition
          .edges
          .insert(edge_key, SparseEdgePartition::with_fitch_subs(subs.clone()));
      }

      Ok(Arc::new(RwLock::new(partition)))
    }
  }
}

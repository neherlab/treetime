#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::commands::prune::run::merge_shared_mutation_branches;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::representation::payload::sparse::{SparseEdgePartition, SparseNodePartition};
  use crate::seq::mutation::Sub;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use approx::assert_relative_eq;
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AsciiChar;

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
    };

    // Build root reference sequence consistent with edge subs.
    // Set each position to the sub's ref character so that edge_subs_from_graph()
    // produces the same mutations as the stored subs.
    let mut ref_seq = treetime_primitives::Seq::from_iter((0..length).map(|_| c(b'A')));
    for (_, _, subs) in edge_mutations {
      for s in subs {
        if s.pos() < length {
          ref_seq[s.pos()] = s.reff();
        }
      }
    }

    // Populate node entries so edge_subs_from_graph() can reconstruct states
    for node in graph.get_nodes() {
      let key = node.read_arc().key();
      let mut node_part = SparseNodePartition::empty(&partition.alphabet);
      node_part.seq.sequence = ref_seq.clone();
      partition.nodes.insert(key, node_part);
    }

    for (source, target, subs) in edge_mutations {
      let edge_key =
        find_edge_key(graph, source, target).unwrap_or_else(|| panic!("edge {source}->{target} not found in graph"));
      partition.edges.insert(
        edge_key,
        SparseEdgePartition {
          subs: subs.clone(),
          ..SparseEdgePartition::default()
        },
      );
    }

    Ok(Arc::new(RwLock::new(partition)))
  }

  /// Find the new internal node (unnamed, non-root, non-leaf).
  fn find_unnamed_internal_nodes(graph: &GraphAncestral) -> Vec<treetime_graph::node::GraphNodeKey> {
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
          Some("A" | "B") => assert_eq!(edge_data.subs.len(), 0, "child should have no remaining mutations"),
          Some("C") => assert_eq!(edge_data.subs.len(), 1),
          None => assert_eq!(
            edge_data.subs.len(),
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
          Some("A") => assert_eq!(edge_data.subs.len(), 0),
          Some("B") => {
            assert_eq!(edge_data.subs.len(), 1);
            assert_eq!(edge_data.subs[0], sub(b'T', 10, b'A'));
          },
          Some("C") => assert_eq!(edge_data.subs.len(), 1),
          None => assert_eq!(edge_data.subs.len(), 2, "internal edge carries shared mutations"),
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
          assert_eq!(edge_data.subs.len(), 3);
        }
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_branch_length_adjustment() -> Result<(), Report> {
    // A and B share 2 mutations out of length 100 => new_bl = 2/100 = 0.02
    // A original bl = 0.1, adjusted = max(0, 0.1 - 0.02) = 0.08
    // B original bl = 0.2, adjusted = max(0, 0.2 - 0.02) = 0.18
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

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      let bl = edge.payload().read_arc().branch_length;

      match target_name.as_deref() {
        None => assert_relative_eq!(bl.unwrap(), 0.02, epsilon = 1e-10),
        Some("A") => assert_relative_eq!(bl.unwrap(), 0.08, epsilon = 1e-6),
        Some("B") => assert_relative_eq!(bl.unwrap(), 0.18, epsilon = 1e-6),
        Some("C") => assert_relative_eq!(bl.unwrap(), 0.3, epsilon = 1e-6),
        _ => {},
      }
    }

    Ok(())
  }

  #[test]
  fn test_merge_branch_length_clamp_to_zero() -> Result<(), Report> {
    // new_bl = 10/100 = 0.1, A original bl = 0.05 => clamped to 0.0
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
        assert_relative_eq!(bl.unwrap(), 0.0, epsilon = 1e-10);
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
    };
    // Populate nodes for p2 with a reference sequence matching the sub ref chars
    let mut p2_ref_seq = treetime_primitives::Seq::from_iter((0..200).map(|_| c(b'A')));
    p2_ref_seq[50] = c(b'C'); // sub C50G uses ref='C'
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
    p2_inner.edges.insert(
      edge_a,
      SparseEdgePartition {
        subs: vec![sub(b'C', 50, b'G')],
        ..SparseEdgePartition::default()
      },
    );
    p2_inner.edges.insert(
      edge_b,
      SparseEdgePartition {
        subs: vec![sub(b'C', 50, b'G')],
        ..SparseEdgePartition::default()
      },
    );
    p2_inner.edges.insert(edge_c, SparseEdgePartition::default());
    let p2 = Arc::new(RwLock::new(p2_inner));

    let partitions = vec![p1, p2];

    let merged = merge_shared_mutation_branches(&mut graph, &partitions)?;
    assert_eq!(merged, 1);
    graph.build()?;

    // Total shared = 1 (from p1) + 1 (from p2) = 2
    // new_bl = 2 / (100 + 200) = 2/300
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      let bl = edge.payload().read_arc().branch_length;

      if target_name.is_none() {
        assert_relative_eq!(bl.unwrap(), 2.0 / 300.0, epsilon = 1e-10);
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
        assert_eq!(p1.edges[&edge.key()].subs.len(), 1);
        assert_eq!(p1.edges[&edge.key()].subs[0], sub(b'G', 5, b'C'));
        // P2: B had {C50G}, shared = {C50G}, remaining = {}
        assert_eq!(p2.edges[&edge.key()].subs.len(), 0);
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

    // New edge bl = 1/1000 = 0.001
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      let bl = edge.payload().read_arc().branch_length;

      if target_name.is_none() {
        assert_relative_eq!(bl.unwrap(), 0.001, epsilon = 1e-10);
      }
    }

    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};

  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::update_marginal;
  use crate::clock::clock_regression::{ClockParams, clock_regression_backward, clock_regression_forward};
  use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::o;
  use crate::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::partition::sparse::{SparseEdgePartition, SparseNodePartition, SparseSeqDistribution};
  use crate::partition::timetree::GraphTimetree;
  use crate::partition::traits::{PartitionRerootOps, PartitionTimetreeAll};
  use crate::payload::timetree::EdgeTimetree;
  use crate::payload::timetree::NodeTimetree;
  use crate::pretty_assert_ulps_eq;
  use crate::seq::indel::InDel;
  use crate::seq::mutation::Sub;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::optimization::reroot::reroot_tree;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::node::{Named, TimeConstraint};
  use treetime_graph::reroot::RerootChanges;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::{AsciiChar, Seq, seq};
  use treetime_utils::make_report;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  fn setup_dates(graph: &GraphTimetree) {
    let dates = btreemap! {
      o!("A") => 2013.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2005.0,
    };

    for n in graph.get_leaves() {
      let name = n.read_arc().payload().read_arc().name().map(|s| s.as_ref().to_owned());
      if let Some(name) = name {
        let date = dates[&name];
        let dist = Arc::new(Distribution::point(date, 1.0));
        n.write_arc().payload().write_arc().set_time_distribution(Some(dist));
      }
    }
  }

  fn gap_free_alignment() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      indoc! {r#"
        >A
        ACGTACGTACGTACGT
        >B
        ACGTACGTACGTACGA
        >C
        ACGTACGTACGTACGG
        >D
        ACGTACGTACGTACGC
      "#},
      &alphabet,
    )
  }

  #[test]
  fn test_reroot_tree_sparse_with_edge_split() -> Result<(), Report> {
    // Test that reroot works correctly with sparse partitions when edge split is enabled
    let aln = gap_free_alignment()?;
    let mut graph: GraphTimetree = nwk_read_str(TREE_NEWICK)?;
    setup_dates(&graph);

    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let fitch = create_fitch_partition(&graph, 0, alphabet, &aln)?;
    let sparse_partition = Arc::new(RwLock::new(fitch.into_marginal_sparse(gtr, &graph)?));

    let clock_params = ClockParams::default();
    clock_regression_backward(&graph, &clock_params, None);
    clock_regression_forward(&graph, &clock_params, None);

    let partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> = sparse_partition;
    let partitions = vec![partition];

    // Record initial state
    let initial_leaf_count = graph.get_leaves().len();
    let initial_node_count = graph.get_nodes().len();

    // Should complete without error - edge split and trivial root removal are now always enabled
    let clock_model = reroot_tree(
      &mut graph,
      &partitions,
      &clock_params,
      None,
      &BranchPointOptimizationParams::default(),
      &RerootSpec::default(),
      true,
    )?;

    // Verify we still have a valid tree with exactly one root
    let root = graph.get_exactly_one_root()?;

    // Root should have no parent edges (it's the true root)
    assert!(
      root.read_arc().inbound().is_empty(),
      "Root should have no inbound edges"
    );

    // Leaf count must be preserved
    assert_eq!(
      graph.get_leaves().len(),
      initial_leaf_count,
      "Leaf count should be unchanged"
    );

    // Node count may increase by 1 if edge was split, but never decrease
    assert!(
      graph.get_nodes().len() >= initial_node_count,
      "Node count should not decrease after reroot"
    );

    // Clock model should have reasonable R² (r_val² > 0.5 for this well-structured tree)
    let r_val = clock_model.r_val().expect("Clock model should have r_val");
    let r_squared = r_val * r_val;
    assert!(r_squared > 0.5, "R² should be > 0.5 for this tree, got {r_squared}");

    // Clock model chisq should be finite and positive
    let chisq = clock_model.chisq().expect("Clock model should have chisq");
    assert!(
      chisq.is_finite() && chisq >= 0.0,
      "Chisq should be finite and non-negative"
    );

    Ok(())
  }

  #[test]
  fn test_sparse_reroot_inverts_subs_and_indels_on_path() -> Result<(), Report> {
    // Tree: (A:0.1,B:0.2)root;
    // After reroot to A, edge direction inverts
    let graph: GraphTimetree = nwk_read_str("(A:0.1,B:0.2)root;")?;

    // Set dates so clock regression works
    for n in graph.get_leaves() {
      let name = n.read_arc().payload().read_arc().name().map(|s| s.as_ref().to_owned());
      if name.as_deref() == Some("A") {
        n.write_arc().payload().write_arc().time = Some(2020.0);
      } else if name.as_deref() == Some("B") {
        n.write_arc().payload().write_arc().time = Some(2010.0);
      }
    }

    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params::default())?;

    let root_key = find_node_key_by_name(&graph, "root").ok_or_else(|| make_report!("root not found"))?;
    let a_key = find_node_key_by_name(&graph, "A").ok_or_else(|| make_report!("A not found"))?;

    // Find edge from root to A
    let edge_to_a_key = graph
      .get_edges()
      .iter()
      .find(|e| {
        let e = e.read_arc();
        let src = e.source();
        let tgt = e.target();
        (src == root_key && tgt == a_key) || (src == a_key && tgt == root_key)
      })
      .map(|e| e.read_arc().key())
      .ok_or_else(|| make_report!("Edge to A not found"))?;

    // Create sparse partition with manually seeded edge data
    let sub_original = Sub::new(c(b'A'), 5_usize, c(b'G'))?; // A5G: ref=A, qry=G
    let indel_original = InDel::del(
      (10, 12),
      seq![
        AsciiChar::from_byte_unchecked(b'A'),
        AsciiChar::from_byte_unchecked(b'C')
      ],
    ); // deletion=true

    let mut sparse_partition = PartitionMarginalSparse {
      index: 0,
      gtr,
      alphabet: alphabet.clone(),
      length: 16,
      root_sequence: seq![AsciiChar::from_byte_unchecked(b'A'); 16],
      nodes: btreemap! {
        root_key => SparseNodePartition::new(&seq![AsciiChar::from_byte_unchecked(b'A'); 16], &alphabet)?,
        a_key => SparseNodePartition::new(&seq![AsciiChar::from_byte_unchecked(b'A'); 16], &alphabet)?,
      },
      edges: btreemap! {
        edge_to_a_key => {
          let mut edge = SparseEdgePartition::with_fitch_subs_and_indels(vec![sub_original], vec![indel_original]);
          edge.msg_from_child = SparseSeqDistribution {
            log_lh: 1.0, // non-default to verify it gets cleared
            ..SparseSeqDistribution::default()
          };
          edge
        },
      },
    };

    // Build RerootChanges with inverted edge keys (simulating reroot from root to A)
    let changes = RerootChanges {
      inverted_edge_keys: vec![edge_to_a_key],
      ..RerootChanges::default()
    };

    // Call apply_reroot directly
    sparse_partition.apply_reroot(&changes)?;

    // Verify substitution is inverted
    let edge_data = &sparse_partition.edges[&edge_to_a_key];
    let sub_after = &edge_data.fitch_subs()[0];
    assert_eq!(sub_after.reff(), c(b'G'), "Sub ref should be swapped to G");
    assert_eq!(sub_after.qry(), c(b'A'), "Sub qry should be swapped to A");

    // Verify indel is inverted
    let indel_after = &edge_data.indels[0];
    assert_eq!(indel_after.deletion, false, "Indel deletion flag should be toggled");

    // Verify msg_from_child is cleared
    pretty_assert_ulps_eq!(edge_data.msg_from_child.log_lh, 0.0, max_ulps = 5);

    Ok(())
  }

  #[test]
  fn test_sparse_reroot_inverts_edge_mutations() -> Result<(), Report> {
    // Test that update_partition_after_reroot inverts edge mutations on the reroot path.
    // Note: This method does NOT compute node sequences - that's done by the subsequent
    // marginal update pass (process_node_backward + process_node_forward).
    //
    // Tree: (A:0.1,B:0.2)root;
    let graph: GraphTimetree = nwk_read_str("(A:0.1,B:0.2)root;")?;

    // Set dates
    for n in graph.get_leaves() {
      let name = n.read_arc().payload().read_arc().name().map(|s| s.as_ref().to_owned());
      if name.as_deref() == Some("A") {
        n.write_arc().payload().write_arc().time = Some(2020.0);
      } else if name.as_deref() == Some("B") {
        n.write_arc().payload().write_arc().time = Some(2010.0);
      }
    }

    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params::default())?;

    let root_key = find_node_key_by_name(&graph, "root").ok_or_else(|| make_report!("root not found"))?;
    let a_key = find_node_key_by_name(&graph, "A").ok_or_else(|| make_report!("A not found"))?;

    // Find edge from root to A
    let edge_to_a_key = graph
      .get_edges()
      .iter()
      .find(|e| {
        let e = e.read_arc();
        let src = e.source();
        let tgt = e.target();
        (src == root_key && tgt == a_key) || (src == a_key && tgt == root_key)
      })
      .map(|e| e.read_arc().key())
      .ok_or_else(|| make_report!("Edge to A not found"))?;

    // Create root sequence with specific characters
    let root_seq = Seq::try_from_slice(b"ACGTACGT")?;

    // Edge has substitution at position 2: root has G, child has T (G->T in parent->child direction)
    let sub = Sub::new(c(b'G'), 2_usize, c(b'T'))?;

    let mut sparse_partition = PartitionMarginalSparse {
      index: 0,
      gtr,
      alphabet: alphabet.clone(),
      length: 8,
      root_sequence: root_seq.clone(),
      nodes: btreemap! {
        root_key => SparseNodePartition::new(&root_seq, &alphabet)?,
        a_key => SparseNodePartition::new(&seq![AsciiChar::from_byte_unchecked(b'A'); 8], &alphabet)?,
      },
      edges: btreemap! {
        edge_to_a_key => SparseEdgePartition::with_fitch_subs(vec![sub]),
      },
    };

    // Build RerootChanges with inverted edge keys (simulating reroot from root to A)
    let changes = RerootChanges {
      inverted_edge_keys: vec![edge_to_a_key],
      ..RerootChanges::default()
    };

    // Call apply_reroot
    sparse_partition.apply_reroot(&changes)?;

    // Verify edge mutation is inverted: was G->T, now should be T->G
    let edge_data = &sparse_partition.edges[&edge_to_a_key];
    assert_eq!(edge_data.fitch_subs().len(), 1);
    let inverted_sub = &edge_data.fitch_subs()[0];
    assert_eq!(
      inverted_sub.reff(),
      AsciiChar::from_byte_unchecked(b'T'),
      "After inversion, reff should be T (was qry)"
    );
    assert_eq!(
      inverted_sub.qry(),
      AsciiChar::from_byte_unchecked(b'G'),
      "After inversion, qry should be G (was reff)"
    );
    assert_eq!(inverted_sub.pos(), 2, "Position should remain unchanged");

    // Verify root_sequence is updated: original root had G at pos 2,
    // child A had T. After reroot to A, new root_sequence should have T at pos 2.
    let expected_new_root_seq = {
      let mut s = root_seq;
      s[2] = c(b'T');
      s
    };
    assert_eq!(
      sparse_partition.root_sequence, expected_new_root_seq,
      "root_sequence should reflect the new root's state after edge inversion"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_root_sequence_updated_with_indel() -> Result<(), Report> {
    let graph: GraphTimetree = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params::default())?;

    let root_key = find_node_key_by_name(&graph, "root").ok_or_else(|| make_report!("root not found"))?;
    let a_key = find_node_key_by_name(&graph, "A").ok_or_else(|| make_report!("A not found"))?;
    let edge_to_a_key = graph
      .get_edges()
      .iter()
      .find(|e| {
        let e = e.read_arc();
        e.source() == root_key && e.target() == a_key
      })
      .map(|e| e.read_arc().key())
      .ok_or_else(|| make_report!("Edge to A not found"))?;

    let root_seq = Seq::try_from_slice(b"ACGTACGT")?;
    let indel = InDel::del((2, 4), seq![c(b'G'), c(b'T')]);

    let mut sparse_partition = PartitionMarginalSparse {
      index: 0,
      gtr,
      alphabet: alphabet.clone(),
      length: 8,
      root_sequence: root_seq.clone(),
      nodes: btreemap! {
        root_key => SparseNodePartition::new(&root_seq, &alphabet)?,
        a_key => SparseNodePartition::new(&seq![c(b'A'); 8], &alphabet)?,
      },
      edges: btreemap! {
        edge_to_a_key => SparseEdgePartition::with_fitch_subs_and_indels(vec![], vec![indel]),
      },
    };

    let changes = RerootChanges {
      inverted_edge_keys: vec![edge_to_a_key],
      ..RerootChanges::default()
    };

    sparse_partition.apply_reroot(&changes)?;

    // Original: root has "ACGTACGT", edge to A has deletion at [2,4) (G,T -> gap).
    // After inversion the indel becomes an insertion. Going from old root to new
    // root in original direction, the child had gaps at [2,4).
    let mut expected = root_seq;
    expected[2] = alphabet.gap();
    expected[3] = alphabet.gap();
    assert_eq!(
      sparse_partition.root_sequence, expected,
      "root_sequence should have gaps at positions 2-3 after indel-based reroot"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_root_sequence_multi_hop() -> Result<(), Report> {
    let graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.2)AB:0.1,C:0.3)root;")?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params::default())?;

    let root_key = find_node_key_by_name(&graph, "root").ok_or_else(|| make_report!("root not found"))?;
    let ab_key = find_node_key_by_name(&graph, "AB").ok_or_else(|| make_report!("AB not found"))?;
    let a_key = find_node_key_by_name(&graph, "A").ok_or_else(|| make_report!("A not found"))?;

    let edge_root_ab = graph
      .get_edges()
      .iter()
      .find(|e| {
        let e = e.read_arc();
        e.source() == root_key && e.target() == ab_key
      })
      .map(|e| e.read_arc().key())
      .ok_or_else(|| make_report!("Edge root->AB not found"))?;
    let edge_ab_a = graph
      .get_edges()
      .iter()
      .find(|e| {
        let e = e.read_arc();
        e.source() == ab_key && e.target() == a_key
      })
      .map(|e| e.read_arc().key())
      .ok_or_else(|| make_report!("Edge AB->A not found"))?;

    let root_seq = Seq::try_from_slice(b"ACGTACGT")?;
    let sub1 = Sub::new(c(b'A'), 0_usize, c(b'G'))?; // root->AB: A0G
    let sub2 = Sub::new(c(b'G'), 0_usize, c(b'T'))?; // AB->A: G0T (cumulative at pos 0)
    let sub3 = Sub::new(c(b'C'), 1_usize, c(b'A'))?; // AB->A: C1A

    let mut sparse_partition = PartitionMarginalSparse {
      index: 0,
      gtr,
      alphabet: alphabet.clone(),
      length: 8,
      root_sequence: root_seq.clone(),
      nodes: btreemap! {
        root_key => SparseNodePartition::new(&root_seq, &alphabet)?,
        ab_key => SparseNodePartition::new(&seq![c(b'A'); 8], &alphabet)?,
        a_key => SparseNodePartition::new(&seq![c(b'A'); 8], &alphabet)?,
      },
      edges: btreemap! {
        edge_root_ab => SparseEdgePartition::with_fitch_subs(vec![sub1]),
        edge_ab_a => SparseEdgePartition::with_fitch_subs(vec![sub2, sub3]),
      },
    };

    // Reroot from root through AB to A: two inverted edges
    let changes = RerootChanges {
      inverted_edge_keys: vec![edge_root_ab, edge_ab_a],
      ..RerootChanges::default()
    };

    sparse_partition.apply_reroot(&changes)?;

    // Original path: root(ACGTACGT) -> AB (pos0: A->G) -> A (pos0: G->T, pos1: C->A)
    // New root = A: pos0 = T, pos1 = A, rest unchanged from root
    let mut expected = root_seq;
    expected[0] = c(b'T');
    expected[1] = c(b'A');
    assert_eq!(
      sparse_partition.root_sequence, expected,
      "root_sequence should reflect cumulative subs across multi-hop reroot"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_tree_sparse_flow_does_not_panic() -> Result<(), Report> {
    // Regression test: verify reroot_tree completes without panicking
    // when keep_root=false (reroot enabled) with sparse partitions
    let aln = gap_free_alignment()?;
    let mut graph: GraphTimetree = nwk_read_str(TREE_NEWICK)?;
    setup_dates(&graph);

    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let fitch = create_fitch_partition(&graph, 0, alphabet, &aln)?;
    let sparse_partition = Arc::new(RwLock::new(fitch.into_marginal_sparse(gtr, &graph)?));

    let clock_params = ClockParams::default();
    clock_regression_backward(&graph, &clock_params, None);
    clock_regression_forward(&graph, &clock_params, None);

    let partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> = sparse_partition;
    let partitions = vec![partition];

    // Record initial state
    let initial_leaf_count = graph.get_leaves().len();

    // Initialize marginal for the sparse partition
    update_marginal(&graph, &partitions)?;

    // First reroot call (simulating keep_root=false flow)
    let clock_model_1 = reroot_tree(
      &mut graph,
      &partitions,
      &clock_params,
      None,
      &BranchPointOptimizationParams::default(),
      &RerootSpec::default(),
      true,
    )?;

    // Verify tree validity after first reroot
    drop(graph.get_exactly_one_root()?);
    assert_eq!(
      graph.get_leaves().len(),
      initial_leaf_count,
      "Leaf count should be unchanged after first reroot"
    );

    let r_squared_1 = clock_model_1.r_val().map(|r| r * r);

    // Second reroot call (simulating refinement iteration)
    let clock_model_2 = reroot_tree(
      &mut graph,
      &partitions,
      &clock_params,
      Some(clock_model_1.clock_rate()),
      &BranchPointOptimizationParams::default(),
      &RerootSpec::default(),
      true,
    )?;

    // Verify tree validity after second reroot
    drop(graph.get_exactly_one_root()?);
    assert_eq!(
      graph.get_leaves().len(),
      initial_leaf_count,
      "Leaf count should be unchanged after second reroot"
    );

    // Second reroot with fixed rate should maintain or improve R²
    // (or have no r_val if rate was fixed)
    if let (Some(r2_1), Some(r2_2)) = (r_squared_1, clock_model_2.r_val().map(|r| r * r)) {
      // Allow small tolerance for floating point
      assert!(
        r2_2 >= r2_1 - 1e-6,
        "Second reroot R² ({r2_2}) should be >= first R² ({r2_1})"
      );
    }

    // Both clock models should have finite chisq
    let chisq_1 = clock_model_1.chisq().expect("First clock model should have chisq");
    assert!(
      chisq_1.is_finite() && chisq_1 >= 0.0,
      "First chisq should be finite and non-negative"
    );

    Ok(())
  }
}

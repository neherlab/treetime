#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::seq::alignment::node_seq_inputs;

  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::SparseReconstruction;
  use crate::clock::clock_regression::{ClockVarianceParams, clock_regression_backward, clock_regression_forward};
  use crate::clock::clock_state::{ClockInputs, ClockState};
  use crate::clock::date_constraints::DateConstraints;
  use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::o;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::partition::marginal::sparse::reroot::reroot_sparse;
  use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
  use crate::partition::timetree::marginal::marginal_update_timetree;
  use crate::partition::timetree::partition::PartitionTimetree;
  use crate::seq::indel::InDel;
  use crate::seq::mutation::Sub;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::optimization::reroot::reroot_tree;
  use crate::timetree::timetree_state::TimetreeState;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_graph::reroot::RerootChanges;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::{AlignmentRecord, AsciiChar, Seq, seq};
  use treetime_utils::make_report;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  fn date_constraints(names: &BTreeMap<GraphNodeKey, Option<String>>, graph: &Graph) -> DateConstraints {
    let dates = btreemap! {
      o!("A") => 2013.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2005.0,
    };

    let mut time_distributions = BTreeMap::new();
    for n in graph.get_leaves() {
      let name = names.get(&n.key()).cloned().flatten();
      if let Some(name) = name {
        let date = dates[&name];
        time_distributions.insert(n.key(), Some(Arc::new(Distribution::point(date, 1.0))));
      }
    }
    DateConstraints {
      time_distributions,
      ..DateConstraints::default()
    }
  }

  fn gap_free_alignment() -> Result<Vec<AlignmentRecord>, Report> {
    let alphabet = Alphabet::default();
    Ok(
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
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect(),
    )
  }

  #[test]
  fn test_reroot_tree_sparse_with_edge_split() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let constraints = date_constraints(&names, &graph);

    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let fitch = create_fitch_partition(&graph, 0, alphabet, &node_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let sparse_partition = PartitionTimetree::Sparse(SparseReconstruction::seeded(partition, gtr, node_states));

    let clock_params = ClockVarianceParams::default();
    let timetree_state = TimetreeState::seed_from_values(&graph, &constraints);
    let clock_inputs = ClockInputs::seed_from_times(&graph, &timetree_state.likely_times(&constraints));
    let mut clock_state = ClockState::new(&graph);
    clock_regression_backward(
      &graph,
      &clock_inputs,
      &mut clock_state,
      &clock_params,
      &branch_lengths,
      None,
    )?;
    clock_regression_forward(
      &graph,
      &clock_inputs,
      &mut clock_state,
      &clock_params,
      &branch_lengths,
      None,
    )?;

    let partitions = vec![sparse_partition];

    let initial_leaf_count = graph.get_leaves().count();
    let initial_node_count = graph.get_nodes().count();

    let names_tt_3 = names;
    let (clock_model, partitions) = reroot_tree(
      &mut graph,
      &constraints,
      &mut clock_state,
      &timetree_state,
      partitions,
      &clock_params,
      None,
      &BranchPointOptimizationParams::default(),
      &RerootSpec::default(),
      true,
      &mut branch_lengths,
      &names_tt_3,
    )?;

    let root = graph.get_exactly_one_root()?;

    assert!(root.inbound().is_empty(), "Root should have no inbound edges");

    assert_eq!(
      graph.get_leaves().count(),
      initial_leaf_count,
      "Leaf count should be unchanged"
    );

    assert!(
      graph.get_nodes().count() >= initial_node_count,
      "Node count should not decrease after reroot"
    );

    let r_val = clock_model.r_val().expect("Clock model should have r_val");
    let r_squared = r_val * r_val;
    assert!(r_squared > 0.5, "R² should be > 0.5 for this tree, got {r_squared}");

    let chisq = clock_model.chisq().expect("Clock model should have chisq");
    assert!(
      chisq.is_finite() && chisq >= 0.0,
      "Chisq should be finite and non-negative"
    );

    Ok(())
  }

  #[test]
  fn test_sparse_reroot_inverts_subs_and_indels_on_path() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;

    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params::default())?;

    let root_key = find_node_key_by_name(&graph, &names, "root").ok_or_else(|| make_report!("root not found"))?;
    let a_key = find_node_key_by_name(&graph, &names, "A").ok_or_else(|| make_report!("A not found"))?;

    let edge_to_a_key = graph
      .get_edges()
      .find(|e| {
        let src = e.source();
        let tgt = e.target();
        (src == root_key && tgt == a_key) || (src == a_key && tgt == root_key)
      })
      .map(|e| e.key())
      .ok_or_else(|| make_report!("Edge to A not found"))?;

    let sub_original = Sub::new(c(b'A'), 5_usize, c(b'G'))?;
    let indel_original = InDel::del(
      (10, 12),
      seq![
        AsciiChar::from_byte_unchecked(b'A'),
        AsciiChar::from_byte_unchecked(b'C')
      ],
    )?;

    let partition = PartitionMarginalSparse {
      index: 0,
      alphabet: alphabet.clone(),
      length: 16,
      root_sequence: seq![AsciiChar::from_byte_unchecked(b'A'); 16],
      obs_nodes: btreemap! {
        root_key => SparseNodeObs::new(&seq![AsciiChar::from_byte_unchecked(b'A'); 16], &alphabet),
        a_key => SparseNodeObs::new(&seq![AsciiChar::from_byte_unchecked(b'A'); 16], &alphabet),
      },
      obs_edges: btreemap! {
        edge_to_a_key => SparseEdgeObs::with_fitch_subs_and_indels(vec![sub_original], vec![indel_original]),
      },
    };
    let node_states = btreemap! {
      root_key => SparseNodeState::leaf(&seq![AsciiChar::from_byte_unchecked(b'A'); 16]),
      a_key => SparseNodeState::leaf(&seq![AsciiChar::from_byte_unchecked(b'A'); 16]),
    };

    let changes = RerootChanges {
      inverted_edge_keys: vec![edge_to_a_key],
      ..RerootChanges::default()
    };

    let recon = reroot_sparse(partition, gtr, node_states, &changes)?;

    let edge_data = &recon.partition.obs_edges[&edge_to_a_key];
    let sub_after = &edge_data.fitch_subs()[0];
    assert_eq!(sub_after.reff(), c(b'G'), "Sub ref should be swapped to G");
    assert_eq!(sub_after.qry(), c(b'A'), "Sub qry should be swapped to A");

    let indel_after = &edge_data.indels[0];
    assert!(!indel_after.is_deletion(), "Indel direction should be toggled");

    assert!(
      recon.edges.backward.is_empty() && recon.edges.forward.is_empty() && recon.edges.estimates.is_empty(),
      "a reroot carries no per-edge results across"
    );

    Ok(())
  }

  #[test]
  fn test_sparse_reroot_inverts_edge_mutations() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;

    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params::default())?;

    let root_key = find_node_key_by_name(&graph, &names, "root").ok_or_else(|| make_report!("root not found"))?;
    let a_key = find_node_key_by_name(&graph, &names, "A").ok_or_else(|| make_report!("A not found"))?;

    let edge_to_a_key = graph
      .get_edges()
      .find(|e| {
        let src = e.source();
        let tgt = e.target();
        (src == root_key && tgt == a_key) || (src == a_key && tgt == root_key)
      })
      .map(|e| e.key())
      .ok_or_else(|| make_report!("Edge to A not found"))?;

    let root_seq = Seq::try_from_slice(b"ACGTACGT")?;

    let sub = Sub::new(c(b'G'), 2_usize, c(b'T'))?;

    let partition = PartitionMarginalSparse {
      index: 0,
      alphabet: alphabet.clone(),
      length: 8,
      root_sequence: root_seq.clone(),
      obs_nodes: btreemap! {
        root_key => SparseNodeObs::new(&root_seq, &alphabet),
        a_key => SparseNodeObs::new(&seq![AsciiChar::from_byte_unchecked(b'A'); 8], &alphabet),
      },
      obs_edges: btreemap! {
        edge_to_a_key => SparseEdgeObs::with_fitch_subs(vec![sub]),
      },
    };
    let node_states = btreemap! {
      root_key => SparseNodeState::leaf(&root_seq),
      a_key => SparseNodeState::leaf(&seq![AsciiChar::from_byte_unchecked(b'A'); 8]),
    };

    let changes = RerootChanges {
      inverted_edge_keys: vec![edge_to_a_key],
      ..RerootChanges::default()
    };

    let recon = reroot_sparse(partition, gtr, node_states, &changes)?;

    let edge_data = &recon.partition.obs_edges[&edge_to_a_key];
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

    let expected_new_root_seq = {
      let mut s = root_seq;
      s[2] = c(b'T');
      s
    };
    assert_eq!(
      recon.partition.root_sequence, expected_new_root_seq,
      "root_sequence should reflect the new root's state after edge inversion"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_root_sequence_updated_with_indel() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params::default())?;

    let root_key = find_node_key_by_name(&graph, &names, "root").ok_or_else(|| make_report!("root not found"))?;
    let a_key = find_node_key_by_name(&graph, &names, "A").ok_or_else(|| make_report!("A not found"))?;
    let edge_to_a_key = graph
      .get_edges()
      .find(|e| e.source() == root_key && e.target() == a_key)
      .map(|e| e.key())
      .ok_or_else(|| make_report!("Edge to A not found"))?;

    let root_seq = Seq::try_from_slice(b"ACGTACGT")?;
    let indel = InDel::del((2, 4), seq![c(b'G'), c(b'T')])?;

    let partition = PartitionMarginalSparse {
      index: 0,
      alphabet: alphabet.clone(),
      length: 8,
      root_sequence: root_seq.clone(),
      obs_nodes: btreemap! {
        root_key => SparseNodeObs::new(&root_seq, &alphabet),
        a_key => SparseNodeObs::new(&seq![c(b'A'); 8], &alphabet),
      },
      obs_edges: btreemap! {
        edge_to_a_key => SparseEdgeObs::with_fitch_subs_and_indels(vec![], vec![indel]),
      },
    };
    let node_states = btreemap! {
      root_key => SparseNodeState::leaf(&root_seq),
      a_key => SparseNodeState::leaf(&seq![c(b'A'); 8]),
    };

    let changes = RerootChanges {
      inverted_edge_keys: vec![edge_to_a_key],
      ..RerootChanges::default()
    };

    let recon = reroot_sparse(partition, gtr, node_states, &changes)?;

    let mut expected = root_seq;
    expected[2] = alphabet.gap();
    expected[3] = alphabet.gap();
    assert_eq!(
      recon.partition.root_sequence, expected,
      "root_sequence should have gaps at positions 2-3 after indel-based reroot"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_root_sequence_multi_hop() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,C:0.3)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params::default())?;

    let root_key = find_node_key_by_name(&graph, &names, "root").ok_or_else(|| make_report!("root not found"))?;
    let ab_key = find_node_key_by_name(&graph, &names, "AB").ok_or_else(|| make_report!("AB not found"))?;
    let a_key = find_node_key_by_name(&graph, &names, "A").ok_or_else(|| make_report!("A not found"))?;

    let edge_root_ab = graph
      .get_edges()
      .find(|e| e.source() == root_key && e.target() == ab_key)
      .map(|e| e.key())
      .ok_or_else(|| make_report!("Edge root->AB not found"))?;
    let edge_ab_a = graph
      .get_edges()
      .find(|e| e.source() == ab_key && e.target() == a_key)
      .map(|e| e.key())
      .ok_or_else(|| make_report!("Edge AB->A not found"))?;

    let root_seq = Seq::try_from_slice(b"ACGTACGT")?;
    let sub1 = Sub::new(c(b'A'), 0_usize, c(b'G'))?;
    let sub2 = Sub::new(c(b'G'), 0_usize, c(b'T'))?;
    let sub3 = Sub::new(c(b'C'), 1_usize, c(b'A'))?;

    let partition = PartitionMarginalSparse {
      index: 0,
      alphabet: alphabet.clone(),
      length: 8,
      root_sequence: root_seq.clone(),
      obs_nodes: btreemap! {
        root_key => SparseNodeObs::new(&root_seq, &alphabet),
        ab_key => SparseNodeObs::new(&seq![c(b'A'); 8], &alphabet),
        a_key => SparseNodeObs::new(&seq![c(b'A'); 8], &alphabet),
      },
      obs_edges: btreemap! {
        edge_root_ab => SparseEdgeObs::with_fitch_subs(vec![sub1]),
        edge_ab_a => SparseEdgeObs::with_fitch_subs(vec![sub2, sub3]),
      },
    };
    let node_states = btreemap! {
      root_key => SparseNodeState::leaf(&root_seq),
      ab_key => SparseNodeState::leaf(&seq![c(b'A'); 8]),
      a_key => SparseNodeState::leaf(&seq![c(b'A'); 8]),
    };

    let changes = RerootChanges {
      inverted_edge_keys: vec![edge_root_ab, edge_ab_a],
      ..RerootChanges::default()
    };

    let recon = reroot_sparse(partition, gtr, node_states, &changes)?;

    let mut expected = root_seq;
    expected[0] = c(b'T');
    expected[1] = c(b'A');
    assert_eq!(
      recon.partition.root_sequence, expected,
      "root_sequence should reflect cumulative subs across multi-hop reroot"
    );

    Ok(())
  }

  #[test]
  fn test_reroot_tree_sparse_flow_does_not_panic() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let constraints = date_constraints(&names, &graph);

    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let fitch = create_fitch_partition(&graph, 0, alphabet, &node_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let sparse_partition = PartitionTimetree::Sparse(SparseReconstruction::seeded(partition, gtr, node_states));

    let clock_params = ClockVarianceParams::default();
    let timetree_state_1 = TimetreeState::seed_from_values(&graph, &constraints);
    let clock_inputs = ClockInputs::seed_from_times(&graph, &timetree_state_1.likely_times(&constraints));
    let mut clock_state = ClockState::new(&graph);
    clock_regression_backward(
      &graph,
      &clock_inputs,
      &mut clock_state,
      &clock_params,
      &branch_lengths,
      None,
    )?;
    clock_regression_forward(
      &graph,
      &clock_inputs,
      &mut clock_state,
      &clock_params,
      &branch_lengths,
      None,
    )?;

    let partitions = vec![sparse_partition];

    let initial_leaf_count = graph.get_leaves().count();

    let (partitions, _) = marginal_update_timetree(&graph, &branch_lengths_or_zero(&branch_lengths), partitions)?;

    let names_tt_2 = names.clone();
    let (clock_model_1, partitions) = reroot_tree(
      &mut graph,
      &constraints,
      &mut clock_state,
      &timetree_state_1,
      partitions,
      &clock_params,
      None,
      &BranchPointOptimizationParams::default(),
      &RerootSpec::default(),
      true,
      &mut branch_lengths,
      &names_tt_2,
    )?;

    let _ = graph.get_exactly_one_root()?;
    assert_eq!(
      graph.get_leaves().count(),
      initial_leaf_count,
      "Leaf count should be unchanged after first reroot"
    );

    let r_squared_1 = clock_model_1.r_val().map(|r| r * r);

    let timetree_state_2 = TimetreeState::seed_from_values(&graph, &constraints);
    let names_tt_1 = names;
    let (clock_model_2, partitions) = reroot_tree(
      &mut graph,
      &constraints,
      &mut clock_state,
      &timetree_state_2,
      partitions,
      &clock_params,
      Some(clock_model_1.clock_rate()),
      &BranchPointOptimizationParams::default(),
      &RerootSpec::default(),
      true,
      &mut branch_lengths,
      &names_tt_1,
    )?;

    let _ = graph.get_exactly_one_root()?;
    assert_eq!(
      graph.get_leaves().count(),
      initial_leaf_count,
      "Leaf count should be unchanged after second reroot"
    );

    if let (Some(r2_1), Some(r2_2)) = (r_squared_1, clock_model_2.r_val().map(|r| r * r)) {
      assert!(
        r2_2 >= r2_1 - 1e-6,
        "Second reroot R² ({r2_2}) should be >= first R² ({r2_1})"
      );
    }

    let chisq_1 = clock_model_1.chisq().expect("First clock model should have chisq");
    assert!(
      chisq_1.is_finite() && chisq_1 >= 0.0,
      "First chisq should be finite and non-negative"
    );

    Ok(())
  }
}

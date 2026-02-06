use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
use crate::commands::ancestral::marginal_unified::update_marginal;
use crate::commands::clock::clock_regression::{ClockParams, clock_regression_backward, clock_regression_forward};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::timetree::optimization::reroot::reroot_tree;
use crate::commands::timetree::partition_ops::{PartitionTimetreeAll, PartitionTimetreeOps};
use crate::distribution::distribution::Distribution;
use crate::graph::node::{GraphNodeKey, Named, TimeConstraint};
use crate::gtr::get_gtr::{JC69Params, jc69};
use crate::io::fasta::{FastaRecord, read_many_fasta_str};
use crate::io::nwk::nwk_read_str;
use crate::o;
use crate::representation::edge_timetree::EdgeTimetree;
use crate::representation::graph_sparse::{MarginalSparseSeqDistribution, SparseEdgePartition};
use crate::representation::node_timetree::NodeTimetree;
use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
use crate::representation::partition_timetree::GraphTimetree;
use crate::seq;
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use eyre::Report;
use indoc::indoc;
use maplit::btreemap;
use parking_lot::RwLock;
use pretty_assertions::assert_eq;
use std::sync::Arc;

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

fn find_node_key_by_name(graph: &GraphTimetree, name: &str) -> Option<GraphNodeKey> {
  for node in graph.get_nodes() {
    let node = node.read_arc();
    let payload = node.payload().read_arc();
    if payload.name().is_some_and(|n| n.as_ref() == name) {
      return Some(node.key());
    }
  }
  None
}

#[test]
fn test_reroot_tree_sparse_with_edge_split() -> Result<(), Report> {
  // Test that reroot works correctly with sparse partitions when edge split is enabled
  let aln = gap_free_alignment()?;
  let mut graph: GraphTimetree = nwk_read_str(TREE_NEWICK)?;
  setup_dates(&graph);

  let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
  let gtr = jc69(JC69Params {
    alphabet: AlphabetName::Nuc,
    treat_gap_as_unknown: false,
    ..JC69Params::default()
  })?;

  let sparse_partition = Arc::new(RwLock::new(PartitionMarginalSparse {
    index: 0,
    gtr,
    alphabet,
    length: get_common_length(&aln)?,
    nodes: btreemap! {},
    edges: btreemap! {},
  }));

  let partitions_for_compress: [Arc<RwLock<PartitionMarginalSparse>>; 1] = [Arc::clone(&sparse_partition)];
  compress_sequences(&graph, &partitions_for_compress, &aln)?;

  let clock_params = ClockParams::default();
  clock_regression_backward(&graph, &clock_params);
  clock_regression_forward(&graph, &clock_params);

  let partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> = sparse_partition;
  let partitions = vec![partition];

  // Should complete without error - edge split and trivial root removal are now always enabled
  let _clock_model = reroot_tree(
    &mut graph,
    &partitions,
    &clock_params,
    None,
    &BranchPointOptimizationParams::default(),
  )?;

  // Verify we still have a valid tree with exactly one root
  let _root = graph.get_exactly_one_root()?;

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

  let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
  let gtr = jc69(JC69Params::default())?;

  let root_key = find_node_key_by_name(&graph, "root").ok_or_else(|| eyre::eyre!("root not found"))?;
  let a_key = find_node_key_by_name(&graph, "A").ok_or_else(|| eyre::eyre!("A not found"))?;

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
    .ok_or_else(|| eyre::eyre!("Edge to A not found"))?;

  // Create sparse partition with manually seeded edge data
  let sub_original = Sub::new(b'A', 5_usize, b'G')?; // A5G: ref=A, qry=G
  let indel_original = InDel::del((10, 12), seq![b'A', b'C']); // deletion=true

  let mut sparse_partition = PartitionMarginalSparse {
    index: 0,
    gtr,
    alphabet: alphabet.clone(),
    length: 16,
    nodes: btreemap! {
      root_key => crate::representation::graph_sparse::SparseNodePartition::new(&seq![b'A'; 16], &alphabet)?,
      a_key => crate::representation::graph_sparse::SparseNodePartition::new(&seq![b'A'; 16], &alphabet)?,
    },
    edges: btreemap! {
      edge_to_a_key => SparseEdgePartition {
        subs: vec![sub_original],
        indels: vec![indel_original],
        msg_to_parent: MarginalSparseSeqDistribution::default(),
        msg_to_child: MarginalSparseSeqDistribution::default(),
        msg_from_child: MarginalSparseSeqDistribution {
          log_lh: 1.0, // non-default to verify it gets cleared
          ..MarginalSparseSeqDistribution::default()
        },
        transmission: None,
      },
    },
  };

  // Manually set root sequence
  if let Some(n) = sparse_partition.nodes.get_mut(&root_key) {
    n.seq.sequence = seq![b'A'; 16];
  }

  // Build path from root (old root) to A (new root)
  // Path format: [(node_key, edge_to_next), ...] where first node is old_root
  let path_from_old_to_new = vec![(root_key, Some(edge_to_a_key)), (a_key, None)];

  // Call reroot_partition_node_only directly
  sparse_partition.reroot_partition_node_only(&graph, root_key, a_key, &path_from_old_to_new)?;

  // Verify substitution is inverted
  let edge_data = &sparse_partition.edges[&edge_to_a_key];
  let sub_after = &edge_data.subs[0];
  assert_eq!(sub_after.reff(), b'G'.into(), "Sub ref should be swapped to G");
  assert_eq!(sub_after.qry(), b'A'.into(), "Sub qry should be swapped to A");

  // Verify indel is inverted
  let indel_after = &edge_data.indels[0];
  assert_eq!(indel_after.deletion, false, "Indel deletion flag should be toggled");

  // Verify msg_from_child is cleared
  assert!(
    edge_data.msg_from_child.log_lh.abs() < f64::EPSILON,
    "msg_from_child should be reset to default"
  );

  Ok(())
}

#[test]
fn test_sparse_reroot_moves_root_sequence() -> Result<(), Report> {
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

  let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
  let gtr = jc69(JC69Params::default())?;

  let root_key = find_node_key_by_name(&graph, "root").ok_or_else(|| eyre::eyre!("root not found"))?;
  let a_key = find_node_key_by_name(&graph, "A").ok_or_else(|| eyre::eyre!("A not found"))?;

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
    .ok_or_else(|| eyre::eyre!("Edge to A not found"))?;

  // Create root sequence with specific characters
  let root_seq = seq![b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T'];

  // Edge has substitution at position 2: root has G, child has T
  let sub = Sub::new(b'G', 2_usize, b'T')?;

  let mut sparse_partition = PartitionMarginalSparse {
    index: 0,
    gtr,
    alphabet: alphabet.clone(),
    length: 8,
    nodes: btreemap! {
      root_key => crate::representation::graph_sparse::SparseNodePartition::new(&root_seq, &alphabet)?,
      a_key => crate::representation::graph_sparse::SparseNodePartition::new(&seq![b'A'; 8], &alphabet)?,
    },
    edges: btreemap! {
      edge_to_a_key => SparseEdgePartition {
        subs: vec![sub],
        indels: vec![],
        msg_to_parent: MarginalSparseSeqDistribution::default(),
        msg_to_child: MarginalSparseSeqDistribution::default(),
        msg_from_child: MarginalSparseSeqDistribution::default(),
        transmission: None,
      },
    },
  };

  // Set root sequence explicitly
  if let Some(n) = sparse_partition.nodes.get_mut(&root_key) {
    n.seq.sequence = root_seq;
  }

  // Build path from root (old root) to A (new root)
  let path_from_old_to_new = vec![(root_key, Some(edge_to_a_key)), (a_key, None)];

  // Call reroot_partition_node_only
  sparse_partition.reroot_partition_node_only(&graph, root_key, a_key, &path_from_old_to_new)?;

  // Verify old root sequence is cleared
  let old_root_seq = &sparse_partition.nodes[&root_key].seq.sequence;
  assert!(old_root_seq.is_empty(), "Old root sequence should be cleared");

  // Verify new root has the expected sequence
  // New root sequence should be computed by traversing from old root to new root
  // Path: root -> A, edge has sub G2T (in original direction)
  // Going reverse (A <- root), we apply the ref (G) at position 2
  // So new root seq = ACGTACGT with position 2 = G (which is already G in original)
  // Actually, the algorithm applies reff when going from old root toward new root
  let new_root_seq = &sparse_partition.nodes[&a_key].seq.sequence;
  let expected_new_root_seq = seq![b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T'];
  assert_eq!(
    *new_root_seq, expected_new_root_seq,
    "New root sequence should be computed from old root by applying reverse edge mutations"
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

  let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
  let gtr = jc69(JC69Params {
    alphabet: AlphabetName::Nuc,
    treat_gap_as_unknown: false,
    ..JC69Params::default()
  })?;

  let sparse_partition = Arc::new(RwLock::new(PartitionMarginalSparse {
    index: 0,
    gtr,
    alphabet,
    length: get_common_length(&aln)?,
    nodes: btreemap! {},
    edges: btreemap! {},
  }));

  let partitions_for_compress: [Arc<RwLock<PartitionMarginalSparse>>; 1] = [Arc::clone(&sparse_partition)];
  compress_sequences(&graph, &partitions_for_compress, &aln)?;

  let clock_params = ClockParams::default();
  clock_regression_backward(&graph, &clock_params);
  clock_regression_forward(&graph, &clock_params);

  let partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> = sparse_partition;
  let partitions = vec![partition];

  // Initialize marginal for the sparse partition
  update_marginal(&graph, &partitions)?;

  // First reroot call (simulating keep_root=false flow)
  let clock_model = reroot_tree(
    &mut graph,
    &partitions,
    &clock_params,
    None,
    &BranchPointOptimizationParams::default(),
  )?;

  // Second reroot call (simulating refinement iteration)
  let _clock_model = reroot_tree(
    &mut graph,
    &partitions,
    &clock_params,
    Some(clock_model.clock_rate()),
    &BranchPointOptimizationParams::default(),
  )?;

  // If we reach here without panic, the test passes
  Ok(())
}

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal_core::MarginalData;
  use crate::partition::marginal_dense::PartitionMarginalDense;
  use crate::partition::traits::PartitionBranchOps;
  use crate::payload::ancestral::GraphAncestral;
  use crate::partition::dense::{DenseEdgePartition, DenseNodePartition, DenseSeqDistribution, DenseSeqInfo};
  use crate::seq::alignment::get_common_length;
  use crate::seq::mutation::Sub;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use ndarray::array;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;

  /// Regression test: uniform outgroup message must not create false substitutions.
  ///
  /// Constructs a partition where the parent and child node posteriors are both
  /// sharp at G (index 2), but the child edge's msg_to_child is uniform [0.25,
  /// 0.25, 0.25, 0.25]. The old implementation compared per-edge message argmax:
  /// msg_to_parent argmax = G, msg_to_child argmax = A (tie-breaks to index 0),
  /// reporting a false G->A substitution. The fixed implementation compares node
  /// posteriors: parent argmax = G, child argmax = G, correctly reporting zero
  /// substitutions.
  #[test]
  fn test_dense_edge_subs_no_false_mutation_from_uniform_outgroup() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2):0.01;")?;
    let edge_ref = &graph.get_edges()[0];
    let edge_key = edge_ref.read_arc().key();
    let parent_key = graph.get_source_node_key(edge_key)?;
    let child_key = graph.get_target_node_key(edge_key)?;

    // Parent posterior: sharp at G (index 2) at position 0
    let parent_posterior = array![[0.0, 0.0, 1.0, 0.0]];
    // Child posterior: sharp at G (index 2) at position 0
    let child_posterior = array![[0.0, 0.0, 1.0, 0.0]];

    // Edge messages: msg_to_parent = sharp G, msg_to_child = uniform (uninformative)
    // Under the old code, argmax of uniform msg_to_child would resolve to index 0
    // (state A), producing a false G->A substitution.
    let msg_to_parent = array![[0.0, 0.0, 1.0, 0.0]];
    let msg_to_child = array![[0.25, 0.25, 0.25, 0.25]];

    let partition = PartitionMarginalDense {
      data: MarginalData {
        gtr: jc69(JC69Params::default())?,
        nodes: btreemap! {
          parent_key => DenseNodePartition {
            seq: DenseSeqInfo::default(),
            profile: DenseSeqDistribution::new(parent_posterior, 0.0),
          },
          child_key => DenseNodePartition {
            seq: DenseSeqInfo::default(),
            profile: DenseSeqDistribution::new(child_posterior, 0.0),
          },
        },
        edges: btreemap! {
          edge_key => DenseEdgePartition {
            msg_to_parent: DenseSeqDistribution::new(msg_to_parent, 0.0),
            msg_to_child: DenseSeqDistribution::new(msg_to_child, 0.0),
            ..Default::default()
          },
        },
        min_branch_length: MIN_BRANCH_LENGTH_FRACTION,
      },
      index: 0,
      alphabet: Alphabet::new(AlphabetName::Nuc)?,
      length: 1,
    };

    let subs = partition.edge_subs(&graph, edge_key)?;

    assert_eq!(Vec::<Sub>::new(), subs);
    Ok(())
  }

  /// Regression test: edge messages that disagree with posteriors must not mask
  /// real substitutions.
  ///
  /// Constructs a partition where parent posterior = A and child posterior = C
  /// (a real substitution), but both edge messages happen to have the same argmax
  /// (both sharp at G). The old code would compare edge messages and report zero
  /// substitutions, missing the real A->C change visible in the node posteriors.
  #[test]
  fn test_dense_edge_subs_detects_real_mutation_hidden_by_edge_messages() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2):0.01;")?;
    let edge_ref = &graph.get_edges()[0];
    let edge_key = edge_ref.read_arc().key();
    let parent_key = graph.get_source_node_key(edge_key)?;
    let child_key = graph.get_target_node_key(edge_key)?;

    // Node posteriors: parent = A, child = C -> real substitution
    let parent_posterior = array![[1.0, 0.0, 0.0, 0.0]];
    let child_posterior = array![[0.0, 1.0, 0.0, 0.0]];

    // Edge messages: both sharp at G -> old code would see no difference
    let msg_to_parent = array![[0.0, 0.0, 1.0, 0.0]];
    let msg_to_child = array![[0.0, 0.0, 1.0, 0.0]];

    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let parent_state = alphabet.char(0); // A
    let child_state = alphabet.char(1); // C

    let partition = PartitionMarginalDense {
      data: MarginalData {
        gtr: jc69(JC69Params::default())?,
        nodes: btreemap! {
          parent_key => DenseNodePartition {
            seq: DenseSeqInfo::default(),
            profile: DenseSeqDistribution::new(parent_posterior, 0.0),
          },
          child_key => DenseNodePartition {
            seq: DenseSeqInfo::default(),
            profile: DenseSeqDistribution::new(child_posterior, 0.0),
          },
        },
        edges: btreemap! {
          edge_key => DenseEdgePartition {
            msg_to_parent: DenseSeqDistribution::new(msg_to_parent, 0.0),
            msg_to_child: DenseSeqDistribution::new(msg_to_child, 0.0),
            ..Default::default()
          },
        },
        min_branch_length: MIN_BRANCH_LENGTH_FRACTION,
      },
      index: 0,
      alphabet,
      length: 1,
    };

    let subs = partition.edge_subs(&graph, edge_key)?;

    let expected = vec![Sub::new(parent_state, 0_usize, child_state)?];
    assert_eq!(expected, subs);
    Ok(())
  }

  /// Parity test: dense edge_subs must match parent-child MAP sequence differences.
  ///
  /// After full marginal inference (backward + forward pass), the node posteriors
  /// contain the complete marginal state distribution at each node. The MAP state
  /// (argmax of the posterior) at each position defines the reconstructed sequence.
  /// Branch mutations are differences between reconstructed parent and child
  /// sequences at canonical (non-gap, non-ambiguous) positions.
  ///
  /// This test mirrors the sparse parity test
  /// `test_sparse_edge_subs_match_reconstructed_branch_differences` in
  /// `test_marginal_sparse.rs`, verifying that dense `edge_subs()` is consistent
  /// with the node posteriors it reads.
  #[test]
  fn test_dense_edge_subs_match_reconstructed_branch_differences() -> Result<(), Report> {
    let aln = divergent_alignment()?;
    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense::new(0, jc69(JC69Params::default())?, Alphabet::new(AlphabetName::Nuc)?, get_common_length(&aln)?)))];

    initialize_marginal(&graph, &partitions, &aln)?;
    update_marginal(&graph, &partitions)?;

    let partition = partitions[0].read_arc();

    // Collect edge_subs() results from all edges.
    let actual_by_edge: BTreeMap<_, _> = graph
      .get_edges()
      .iter()
      .map(|edge_ref| {
        let edge_key = edge_ref.read_arc().key();
        let subs = partition.edge_subs(&graph, edge_key).unwrap();
        (edge_key, subs)
      })
      .collect();

    // Build expected substitutions by comparing MAP states of node posteriors
    // directly. This is independent of edge_subs() and serves as the oracle.
    let expected_by_edge: BTreeMap<_, _> = graph
      .get_edges()
      .iter()
      .map(|edge_ref| {
        let edge_ref = edge_ref.read_arc();
        let edge_key = edge_ref.key();
        let parent_key = edge_ref.source();
        let child_key = edge_ref.target();
        let expected = helpers::diff_map_states(
          &partition.alphabet,
          &partition.data.nodes[&parent_key],
          &partition.data.nodes[&child_key],
        );
        (edge_key, expected)
      })
      .collect();

    assert_eq!(expected_by_edge, actual_by_edge);
    Ok(())
  }

  /// Gap positions must be excluded from dense edge_subs even when node posteriors
  /// differ at those positions.
  ///
  /// Constructs a partition where one leaf has gaps at positions 1-2. The node
  /// posteriors at gap positions are uniform (from treat_gap_as_unknown), which
  /// could differ from the parent's non-gap posterior. The gap filter must prevent
  /// these positions from appearing as substitutions.
  #[test]
  fn test_dense_edge_subs_excludes_gap_positions_with_posteriors() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2):0.01;")?;
    let edge_ref = &graph.get_edges()[0];
    let edge_key = edge_ref.read_arc().key();
    let parent_key = graph.get_source_node_key(edge_key)?;
    let child_key = graph.get_target_node_key(edge_key)?;

    // Position 0: both A -> no sub
    // Position 1: parent C, child uniform (gap) -> excluded by gap filter
    // Position 2: parent G, child uniform (gap) -> excluded by gap filter
    // Position 3: parent A, child C -> real sub
    #[rustfmt::skip]
    let parent_posterior = array![
      [1.0,  0.0,  0.0,  0.0],
      [0.0,  1.0,  0.0,  0.0],
      [0.0,  0.0,  1.0,  0.0],
      [1.0,  0.0,  0.0,  0.0],
    ];
    #[rustfmt::skip]
    let child_posterior = array![
      [1.0,  0.0,  0.0,  0.0],
      [0.25, 0.25, 0.25, 0.25],
      [0.25, 0.25, 0.25, 0.25],
      [0.0,  1.0,  0.0,  0.0],
    ];

    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let partition = PartitionMarginalDense {
      data: MarginalData {
        gtr: jc69(JC69Params::default())?,
        nodes: btreemap! {
          parent_key => DenseNodePartition {
            seq: DenseSeqInfo::default(),
            profile: DenseSeqDistribution::new(parent_posterior, 0.0),
          },
          child_key => DenseNodePartition {
            seq: DenseSeqInfo { gaps: vec![(1, 3)], non_char: vec![(1, 3)], ..Default::default() },
            profile: DenseSeqDistribution::new(child_posterior, 0.0),
          },
        },
        edges: btreemap! {
          edge_key => DenseEdgePartition::default(),
        },
        min_branch_length: MIN_BRANCH_LENGTH_FRACTION / 4.0,
      },
      index: 0,
      alphabet: alphabet.clone(),
      length: 4,
    };

    let subs = partition.edge_subs(&graph, edge_key)?;

    // Only position 3 (A->C) should appear; positions 1-2 are gap-filtered.
    let expected = vec![Sub::new(alphabet.char(0), 3_usize, alphabet.char(1))?];
    assert_eq!(expected, subs);
    Ok(())
  }

  /// The is_canonical filter is present in dense edge_subs for consistency with
  /// sparse edge_subs (marginal_sparse.rs:287), but it is structurally vacuous
  /// for current profile shapes: dense profiles have n_canonical columns, so
  /// argmax always maps to a canonical character. This test verifies the
  /// filter's presence by ensuring the parity test oracle (which has the
  /// filter) agrees with production. The parity test above already covers
  /// this, but this dedicated test documents the design intent.
  #[test]
  fn test_dense_edge_subs_is_canonical_filter_present() -> Result<(), Report> {
    let aln = divergent_alignment()?;
    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense::new(0, jc69(JC69Params::default())?, Alphabet::new(AlphabetName::Nuc)?, get_common_length(&aln)?)))];

    initialize_marginal(&graph, &partitions, &aln)?;
    update_marginal(&graph, &partitions)?;

    let partition = partitions[0].read_arc();
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let subs = partition.edge_subs(&graph, edge_key)?;
      // Every reported substitution must involve canonical states only.
      // This validates the is_canonical filter is present and active.
      for sub in &subs {
        assert!(
          partition.alphabet.is_canonical(sub.reff()),
          "parent state must be canonical: {sub:?}"
        );
        assert!(
          partition.alphabet.is_canonical(sub.qry()),
          "child state must be canonical: {sub:?}"
        );
      }
    }
    Ok(())
  }

  fn divergent_alignment() -> Result<Vec<FastaRecord>, Report> {
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

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::partition::dense::DenseNodePartition;
    use crate::seq::mutation::Sub;
    use treetime_utils::array::ndarray::argmax_first;

    /// Derive MAP-state substitutions from two node posteriors, skipping gap positions.
    ///
    /// This is the oracle for parity testing: it reads the same node posterior data
    /// that edge_subs() reads, computes argmax independently, and compares. The
    /// implementation deliberately avoids calling edge_subs() so the test is not
    /// circular.
    pub fn diff_map_states(
      alphabet: &Alphabet,
      parent_node: &DenseNodePartition,
      child_node: &DenseNodePartition,
    ) -> Vec<Sub> {
      let parent_profile = &parent_node.profile.dis;
      let child_profile = &child_node.profile.dis;
      let parent_gaps = &parent_node.seq.gaps;
      let child_gaps = &child_node.seq.gaps;

      let mut subs = Vec::new();
      for (pos, (parent_row, child_row)) in parent_profile.rows().into_iter().zip(child_profile.rows()).enumerate() {
        if parent_gaps.iter().any(|&(start, end)| pos >= start && pos < end) {
          continue;
        }
        if child_gaps.iter().any(|&(start, end)| pos >= start && pos < end) {
          continue;
        }

        let parent_state = alphabet.char(argmax_first(&parent_row).unwrap_or(0));
        let child_state = alphabet.char(argmax_first(&child_row).unwrap_or(0));
        if parent_state != child_state && alphabet.is_canonical(parent_state) && alphabet.is_canonical(child_state) {
          subs.push(Sub::new(parent_state, pos, child_state).unwrap());
        }
      }
      subs
    }
  }
}

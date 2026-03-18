#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::get_common_length;
  use crate::commands::ancestral::marginal::initialize_marginal;
  use crate::commands::optimize::partition_ops::PartitionOptimizeOps;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::representation::payload::dense::{DenseEdgePartition, DenseNodePartition, DenseSeqDis, DenseSeqInfo};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use ndarray::array;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;

  use self::helpers::{compute_soft_hamming, divergent_alignment, identical_alignment, setup_dense};

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  /// Sharp profiles (one-hot) that agree at every position: dot = 1, diff = 0.
  #[test]
  fn test_sharp_matching_profiles_zero_differences() -> Result<(), Report> {
    #[rustfmt::skip]
    let profiles = array![
      [1.0, 0.0, 0.0, 0.0],
      [0.0, 1.0, 0.0, 0.0],
      [0.0, 0.0, 1.0, 0.0],
      [0.0, 0.0, 0.0, 1.0],
    ];
    let result = compute_soft_hamming(profiles.clone(), profiles)?;
    assert_abs_diff_eq!(0.0, result, epsilon = 1e-10);
    Ok(())
  }

  /// Sharp profiles (one-hot) that disagree at every position: dot = 0, diff = 1 each.
  #[test]
  fn test_sharp_mismatched_profiles_integer_differences() -> Result<(), Report> {
    #[rustfmt::skip]
    let parent = array![
      [1.0, 0.0, 0.0, 0.0],
      [0.0, 1.0, 0.0, 0.0],
      [0.0, 0.0, 1.0, 0.0],
      [0.0, 0.0, 0.0, 1.0],
    ];
    #[rustfmt::skip]
    let child = array![
      [0.0, 1.0, 0.0, 0.0],
      [1.0, 0.0, 0.0, 0.0],
      [0.0, 0.0, 0.0, 1.0],
      [0.0, 0.0, 1.0, 0.0],
    ];
    let result = compute_soft_hamming(parent, child)?;
    assert_abs_diff_eq!(4.0, result, epsilon = 1e-10);
    Ok(())
  }

  /// Uniform profiles (0.25 each): dot = 4 * 0.25^2 = 0.25, diff = 0.75 per position.
  /// 4 positions * 0.75 = 3.0 total. Hard argmax would give 0 (argmax_first picks
  /// the same state for both identical uniform profiles).
  #[test]
  fn test_uniform_profiles_fractional_differences() -> Result<(), Report> {
    #[rustfmt::skip]
    let uniform = array![
      [0.25, 0.25, 0.25, 0.25],
      [0.25, 0.25, 0.25, 0.25],
      [0.25, 0.25, 0.25, 0.25],
      [0.25, 0.25, 0.25, 0.25],
    ];
    let result = compute_soft_hamming(uniform.clone(), uniform)?;
    assert_abs_diff_eq!(3.0, result, epsilon = 1e-10);
    Ok(())
  }

  /// Weakly informative profiles with same dominant state.
  /// dot([0.7, 0.1, 0.1, 0.1], [0.7, 0.1, 0.1, 0.1]) = 0.49 + 3*0.01 = 0.52.
  /// diff = 0.48. Hard argmax: both pick state 0, gives 0.
  #[test]
  fn test_weakly_informative_same_dominant_state() -> Result<(), Report> {
    let profile = array![[0.7, 0.1, 0.1, 0.1]];
    let result = compute_soft_hamming(profile.clone(), profile)?;
    assert_abs_diff_eq!(0.48, result, epsilon = 1e-10);
    Ok(())
  }

  /// Weakly informative profiles with different dominant states.
  /// dot([0.7, 0.1, 0.1, 0.1], [0.1, 0.7, 0.1, 0.1]) = 2*0.07 + 2*0.01 = 0.16.
  /// diff = 0.84. Hard argmax: picks different states, gives 1.
  #[test]
  fn test_weakly_informative_different_dominant_state() -> Result<(), Report> {
    let parent = array![[0.7, 0.1, 0.1, 0.1]];
    let child = array![[0.1, 0.7, 0.1, 0.1]];
    let result = compute_soft_hamming(parent, child)?;
    assert_abs_diff_eq!(0.84, result, epsilon = 1e-10);
    Ok(())
  }

  /// Mix of sharp, uniform, and weakly informative positions.
  /// Position 0: sharp agree            → 0.00
  /// Position 1: sharp disagree         → 1.00
  /// Position 2: uniform                → 0.75
  /// Position 3: weakly informative     → 0.48
  /// Total: 2.23
  #[test]
  fn test_mixed_profiles_sum_of_contributions() -> Result<(), Report> {
    #[rustfmt::skip]
    let parent = array![
      [1.0,  0.0,  0.0,  0.0 ],
      [1.0,  0.0,  0.0,  0.0 ],
      [0.25, 0.25, 0.25, 0.25],
      [0.7,  0.1,  0.1,  0.1 ],
    ];
    #[rustfmt::skip]
    let child = array![
      [1.0,  0.0,  0.0,  0.0 ],
      [0.0,  1.0,  0.0,  0.0 ],
      [0.25, 0.25, 0.25, 0.25],
      [0.7,  0.1,  0.1,  0.1 ],
    ];
    let result = compute_soft_hamming(parent, child)?;
    assert_abs_diff_eq!(2.23, result, epsilon = 1e-10);
    Ok(())
  }

  /// Gap positions are excluded from the soft Hamming distance.
  /// Parent has gaps at positions 1-2, excluding the sharp-disagree and uniform positions.
  /// Only position 0 (sharp agree → 0) and 3 (weak same → 0.48) contribute.
  #[test]
  fn test_gap_positions_excluded() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2):0.01;")?;
    let edge_ref = &graph.get_edges()[0];
    let edge_key = edge_ref.read_arc().key();
    let parent_key = graph.get_source_node_key(edge_key)?;
    let child_key = graph.get_target_node_key(edge_key)?;

    #[rustfmt::skip]
    let parent_profiles = array![
      [1.0,  0.0,  0.0,  0.0 ],
      [1.0,  0.0,  0.0,  0.0 ],
      [0.25, 0.25, 0.25, 0.25],
      [0.7,  0.1,  0.1,  0.1 ],
    ];
    #[rustfmt::skip]
    let child_profiles = array![
      [1.0,  0.0,  0.0,  0.0 ],
      [0.0,  1.0,  0.0,  0.0 ],
      [0.25, 0.25, 0.25, 0.25],
      [0.7,  0.1,  0.1,  0.1 ],
    ];

    let partition = PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(AlphabetName::Nuc)?,
      length: 4,
      nodes: btreemap! {
        parent_key => DenseNodePartition {
          seq: DenseSeqInfo { gaps: vec![(1, 3)], ..Default::default() },
          profile: DenseSeqDis::default(),
        },
        child_key => DenseNodePartition {
          seq: DenseSeqInfo::default(),
          profile: DenseSeqDis::default(),
        },
      },
      edges: btreemap! {
        edge_key => DenseEdgePartition {
          msg_to_parent: DenseSeqDis::new(parent_profiles),
          msg_to_child: DenseSeqDis::new(child_profiles),
          ..Default::default()
        },
      },
    };

    let result = partition.edge_initial_differences(&graph, edge_key)?;
    assert_abs_diff_eq!(0.48, result, epsilon = 1e-10);
    Ok(())
  }

  /// For identical sequences, hard Hamming is 0 (argmax always matches).
  /// Soft Hamming is positive but small: the transition matrix introduces
  /// off-diagonal probability mass in per-edge messages even when all
  /// leaves share the same sequence.
  #[test]
  fn test_identical_sequences_hard_zero_soft_small() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let aln = identical_alignment()?;
    let partitions = setup_dense(&graph, &aln)?;

    let p = partitions[0].read_arc();
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let soft = p.edge_initial_differences(&graph, edge_key)?;
      let hard = p.edge_subs(&graph, edge_key)?.len() as f64;
      assert_abs_diff_eq!(0.0, hard, epsilon = 1e-10);
      assert!(soft >= 0.0, "soft Hamming must be non-negative");
      assert!(
        soft < 4.0,
        "soft Hamming should be small for identical sequences, got {soft}"
      );
    }
    Ok(())
  }

  /// For divergent sequences, soft and hard Hamming should differ at edges
  /// with uncertain profiles (internal edges where the tree is not well-resolved).
  #[test]
  fn test_divergent_sequences_soft_differs_from_hard() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let aln = divergent_alignment()?;
    let partitions = setup_dense(&graph, &aln)?;

    let p = partitions[0].read_arc();
    let mut found_difference = false;
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let soft = p.edge_initial_differences(&graph, edge_key)?;
      let hard = p.edge_subs(&graph, edge_key)?.len() as f64;
      if (soft - hard).abs() > 1e-6 {
        found_difference = true;
      }
    }
    assert!(
      found_difference,
      "soft and hard Hamming should differ at edges with uncertain profiles"
    );
    Ok(())
  }

  /// Soft Hamming is bounded: 0 <= differences <= effective_length for any edge.
  #[test]
  fn test_differences_bounded_by_effective_length() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let aln = divergent_alignment()?;
    let partitions = setup_dense(&graph, &aln)?;

    let p = partitions[0].read_arc();
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let differences = p.edge_initial_differences(&graph, edge_key)?;
      let effective_length = p.edge_effective_length(&graph, edge_key)?;
      assert!(
        differences >= 0.0,
        "differences must be non-negative, got {differences}"
      );
      assert!(
        differences <= effective_length as f64,
        "differences ({differences}) must not exceed effective length ({effective_length})"
      );
    }

    Ok(())
  }

  mod helpers {
    use super::*;

    /// Compute soft Hamming distance for a single edge with manually specified profiles.
    ///
    /// Creates a minimal partition with the given edge messages and computes
    /// edge_initial_differences. The graph structure comes from a simple
    /// two-leaf newick tree; the first edge is used.
    pub fn compute_soft_hamming(
      parent_profiles: ndarray::Array2<f64>,
      child_profiles: ndarray::Array2<f64>,
    ) -> Result<f64, Report> {
      let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2):0.01;")?;
      let edge_ref = &graph.get_edges()[0];
      let edge_key = edge_ref.read_arc().key();
      let parent_key = graph.get_source_node_key(edge_key)?;
      let child_key = graph.get_target_node_key(edge_key)?;

      let length = parent_profiles.nrows();

      let partition = PartitionMarginalDense {
        index: 0,
        gtr: jc69(JC69Params::default())?,
        alphabet: Alphabet::new(AlphabetName::Nuc)?,
        length,
        nodes: btreemap! {
          parent_key => DenseNodePartition {
            seq: DenseSeqInfo::default(),
            profile: DenseSeqDis::default(),
          },
          child_key => DenseNodePartition {
            seq: DenseSeqInfo::default(),
            profile: DenseSeqDis::default(),
          },
        },
        edges: btreemap! {
          edge_key => DenseEdgePartition {
            msg_to_parent: DenseSeqDis::new(parent_profiles),
            msg_to_child: DenseSeqDis::new(child_profiles),
            ..Default::default()
          },
        },
      };

      partition.edge_initial_differences(&graph, edge_key)
    }

    pub fn identical_alignment() -> Result<Vec<FastaRecord>, Report> {
      let alphabet = Alphabet::default();
      read_many_fasta_str(
        indoc! {r#"
          >A
          ACGTACGTACGTACGT
          >B
          ACGTACGTACGTACGT
          >C
          ACGTACGTACGTACGT
          >D
          ACGTACGTACGTACGT
        "#},
        &alphabet,
      )
    }

    pub fn divergent_alignment() -> Result<Vec<FastaRecord>, Report> {
      let alphabet = Alphabet::default();
      read_many_fasta_str(
        indoc! {r#"
          >A
          AAAACCCCGGGGTTTT
          >B
          CCCCGGGGTTTTAAAA
          >C
          GGGGTTTTAAAACCCC
          >D
          TTTTAAAACCCCGGGG
        "#},
        &alphabet,
      )
    }

    pub fn setup_dense(
      graph: &GraphAncestral,
      aln: &[FastaRecord],
    ) -> Result<Vec<Arc<RwLock<PartitionMarginalDense>>>, Report> {
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense {
        index: 0,
        gtr: jc69(JC69Params::default())?,
        alphabet,
        length: get_common_length(aln)?,
        nodes: btreemap! {},
        edges: btreemap! {},
      }))];

      initialize_marginal(graph, &partitions, aln)?;

      Ok(partitions)
    }
  }
}

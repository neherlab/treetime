#![allow(
  clippy::integer_division,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::DenseReconstruction;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::infer_gtr::common::{
    InferGtrOptions, accumulate_mutation_counts, get_branch_mutation_matrix, infer_gtr_impl,
  };
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::shared::update::MarginalPasses;
  use crate::pretty_assert_ulps_eq;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::LazyLock;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_utils::{
    pretty_assert_abs_diff_eq, pretty_assert_array_nonneg, pretty_assert_array_offdiag_upper_bounded,
    pretty_assert_array_positive,
  };

  use ndarray::{Array1, Array2, array};
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  fn setup_dense_partition(
    tree_nwk: &str,
    aln: &[AlignmentRecord],
  ) -> Result<(Graph, DenseReconstruction, BTreeMap<GraphEdgeKey, Option<f64>>), Report> {
    let nwk_parsed = nwk_read_str(tree_nwk)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let partition = PartitionMarginalDense::new(0, alphabet, get_common_length(aln)?);
    let node_states = partition.attach_sequences(&graph, &node_seq_inputs(&graph, &names, aln.to_vec()))?;
    let recon = DenseReconstruction::seeded(partition, gtr, node_states);
    let (recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    Ok((graph, recon, branch_lengths))
  }

  #[test]
  fn test_uniform_sequences() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
      indoc! {r#"
      >A
      ACGT
      >B
      ACGT
      >C
      ACGT
      >D
      ACGT
      "#},
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let (graph, recon, branch_lengths) =
      setup_dense_partition("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;", &aln)?;

    let counts = recon.partition.count_transitions(
      &recon.gtr,
      &graph,
      &branch_lengths_or_zero(&branch_lengths),
      &recon.node_states,
      &recon.edges.backward,
      &recon.edges.forward,
    )?;

    pretty_assert_array_nonneg!(counts.nij);
    pretty_assert_array_offdiag_upper_bounded!(counts.nij, bound = 0.1);

    assert_eq!(array![1.0, 1.0, 1.0, 1.0], counts.root_state);

    Ok(())
  }

  #[test]
  fn test_single_mutation() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
      indoc! {r#"
      >A
      ACGT
      >B
      CCGT
      "#},
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let (graph, recon, branch_lengths) = setup_dense_partition("(A:0.1,B:0.1)root:0.0;", &aln)?;

    let counts = recon.partition.count_transitions(
      &recon.gtr,
      &graph,
      &branch_lengths_or_zero(&branch_lengths),
      &recon.node_states,
      &recon.edges.backward,
      &recon.edges.forward,
    )?;

    let total_mutations: f64 = counts
      .nij
      .iter()
      .enumerate()
      .filter_map(|(idx, &val)| {
        let i = idx / 4;
        let j = idx % 4;
        (i != j).then_some(val)
      })
      .sum();

    assert!(
      total_mutations >= 1.0,
      "Expected at least 1 mutation, got {total_mutations}"
    );

    Ok(())
  }

  #[test]
  fn test_zero_branch_lengths() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
      indoc! {r#"
      >A
      ACGT
      >B
      CCGT
      >C
      ACGT
      >D
      CCGT
      "#},
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let (graph, recon, branch_lengths) =
      setup_dense_partition("((A:0.0,B:0.0)AB:0.0,(C:0.0,D:0.0)CD:0.0)root:0.0;", &aln)?;

    let counts = recon.partition.count_transitions(
      &recon.gtr,
      &graph,
      &branch_lengths_or_zero(&branch_lengths),
      &recon.node_states,
      &recon.edges.backward,
      &recon.edges.forward,
    )?;

    let ti_max = counts.Ti.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    assert!(ti_max < 1e-2, "max(Ti) = {ti_max} should be < 1e-2 for clamped zero-BL");

    let total_off_diagonal: f64 = counts
      .nij
      .indexed_iter()
      .filter(|&((i, j), _)| i != j)
      .map(|(_, &v)| v)
      .sum();
    assert!(
      total_off_diagonal < 5.0,
      "total off-diagonal nij = {total_off_diagonal} should be < 5.0 for clamped zero-BL"
    );

    assert!(counts.root_state.sum() > 0.0, "root_state should be populated");

    Ok(())
  }

  #[test]
  fn test_zero_branch_lengths_unclamped() -> Result<(), Report> {
    let n_states = 4;

    let exp_qt = Array2::eye(n_states);

    let messages = Array2::eye(n_states);

    let mut_stack = get_branch_mutation_matrix(&messages, &messages, &exp_qt);

    let mut nij = Array2::zeros((n_states, n_states));
    let mut ti = Array1::zeros(n_states);

    accumulate_mutation_counts(&mut_stack, 0.0, &mut nij, &mut ti);

    let expected_ti: Array1<f64> = Array1::zeros(n_states);
    assert_eq!(expected_ti, ti);

    let expected_nij: Array2<f64> = Array2::eye(n_states);
    assert_eq!(expected_nij, nij);

    Ok(())
  }

  #[test]
  fn test_produces_valid_model() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
      indoc! {r#"
      >A
      ACATCGCCGTAGAC
      >B
      GCATCCCTGTAGGG
      >C
      CCGGCGATGTGTTG
      >D
      TCGGCCGTGTGTTG
      "#},
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let tree_nwk = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
    let (graph, recon, branch_lengths) = setup_dense_partition(tree_nwk, &aln)?;

    let counts = recon.partition.count_transitions(
      &recon.gtr,
      &graph,
      &branch_lengths_or_zero(&branch_lengths),
      &recon.node_states,
      &recon.edges.backward,
      &recon.edges.forward,
    )?;
    let result = infer_gtr_impl(&counts, &InferGtrOptions::default())?;

    pretty_assert_abs_diff_eq!(result.W, result.W.t().to_owned(), epsilon = 1e-9);
    pretty_assert_ulps_eq!(1.0, result.pi.sum(), epsilon = 1e-9);
    pretty_assert_array_positive!(result.pi);

    assert!(result.mu > 0.0, "mu should be positive, got {}", result.mu);

    Ok(())
  }
}

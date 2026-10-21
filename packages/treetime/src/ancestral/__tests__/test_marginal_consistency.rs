#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, Mu};
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction, SparseReconstruction};
  use crate::partition::marginal::sample::SampleMode;
  use crate::pretty_assert_ulps_eq;
  use crate::seq::mutation::Sub;
  use crate::test_utils::get_common_length;
  use crate::test_utils::leaf_seq_inputs;
  use crate::test_utils::{
    dense_reconstruction, emitted_sequences_by_name, find_node_key_by_name, internal_node_keys,
    marginal_update_reconstruction, sparse_reconstruction,
  };
  use eyre::Report;
  use indoc::indoc;
  use itertools::Itertools;
  use ndarray::{Array1, array, aview1};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::sync::LazyLock;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node_names::NodeNames;

  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::{AlignmentRecord, Seq};
  use treetime_utils::sync::random::get_random_number_generator;
  use treetime_utils::{assert_error, make_report};

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  fn gap_free_alignment() -> Result<Vec<AlignmentRecord>, Report> {
    Ok(
      fasta_read(
        indoc! {r#"
      >A
      ACGTACGTACGTACGT
      >B
      ACGTACGTACGTACGA
      >C
      ACGTACGTACGTACGG
      >D
      ACGTACGTACGTACGC
    "#}
        .as_bytes(),
        &*NUC_ALPHABET,
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect(),
    )
  }

  fn ambiguous_r_in_g_clade_alignment() -> Result<Vec<AlignmentRecord>, Report> {
    Ok(
      fasta_read(
        indoc! {r#"
      >A
      RCGTACGT
      >B
      GCGTACGT
      >C
      GCGTACGT
      >D
      GCGTACGT
    "#}
        .as_bytes(),
        &*NUC_ALPHABET,
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect(),
    )
  }

  fn run_dense_marginal(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    names: &NodeNames,
    aln: &[AlignmentRecord],
    gtr: GTR,
  ) -> Result<(f64, MarginalReconstruction), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let partition = PartitionMarginalDense::new(alphabet, graph, &leaf_seq_inputs(graph, names, aln.to_vec()))?;
    let recon = MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, gtr));
    let (recon, log_lh) = marginal_update_reconstruction(recon, graph, &branch_lengths_or_zero(branch_lengths))?;
    let log_lh = log_lh.value();
    Ok((log_lh, recon))
  }

  fn run_sparse_marginal(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    names: &NodeNames,
    aln: &[AlignmentRecord],
    gtr: GTR,
  ) -> Result<(f64, MarginalReconstruction), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let fitch = create_fitch_partition(graph, alphabet, leaf_seq_inputs(graph, names, aln.to_vec()))?;
    let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
    let recon = MarginalReconstruction::Sparse(SparseReconstruction::seeded(partition, gtr, node_states));
    let (recon, log_lh) = marginal_update_reconstruction(recon, graph, &branch_lengths_or_zero(branch_lengths))?;
    let log_lh = log_lh.value();
    Ok((log_lh, recon))
  }

  #[test]
  fn test_marginal_dense_sparse_log_lh_consistency_gap_free() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let gtr_dense = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let gtr_sparse = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let (log_lh_dense, _) = run_dense_marginal(&graph, &branch_lengths, &names, &aln, gtr_dense)?;
    let (log_lh_sparse, _) = run_sparse_marginal(&graph, &branch_lengths, &names, &aln, gtr_sparse)?;

    pretty_assert_ulps_eq!(log_lh_dense, log_lh_sparse, epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_marginal_sparse_varpos_matches_dense_profile_gap_free() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let gtr_dense = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let gtr_sparse = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let (_, dense_partition) = run_dense_marginal(&graph, &branch_lengths, &names, &aln, gtr_dense)?;
    let (_, sparse_partition) = run_sparse_marginal(&graph, &branch_lengths, &names, &aln, gtr_sparse)?;

    let root_key = find_node_key_by_name(&graph, &names, "root").ok_or_else(|| make_report!("Root node not found"))?;
    let ab_key = find_node_key_by_name(&graph, &names, "AB").ok_or_else(|| make_report!("AB node not found"))?;

    let dense = dense_reconstruction(&dense_partition);
    let sparse = sparse_reconstruction(&sparse_partition);

    for node_key in [root_key, ab_key] {
      let dense_node = &dense.node_states[&node_key];
      let sparse_node = &sparse.node_states[&node_key];

      for (&pos, var_pos) in &sparse_node.profile.variable {
        let dense_row = dense_node.profile.dis.row(pos).to_owned();
        pretty_assert_ulps_eq!(dense_row, aview1(var_pos.dis).to_owned(), epsilon = 1e-6);
      }
    }

    Ok(())
  }

  #[test]
  fn test_marginal_dense_sparse_ambiguous_character_expectations_documented() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = fasta_read(
      indoc! {r#"
      >A
      ACGTACGTACGTACNT
      >B
      ACGTACGTACGTACRA
      >C
      ACGTACGTACGTACGG
      >D
      ACGTACGTACGTACGC
    "#}
      .as_bytes(),
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;

    let gtr_dense = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let gtr_sparse = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let (log_lh_dense, dense_partition) = run_dense_marginal(&graph, &branch_lengths, &names, &aln, gtr_dense)?;
    let (log_lh_sparse, sparse_partition) = run_sparse_marginal(&graph, &branch_lengths, &names, &aln, gtr_sparse)?;

    pretty_assert_ulps_eq!(log_lh_dense, log_lh_sparse, epsilon = 1e-10);

    let dense = dense_reconstruction(&dense_partition);
    let sparse = sparse_reconstruction(&sparse_partition);

    for node_data in dense.node_states.values() {
      if !node_data.profile.dis.is_empty() {
        for row in node_data.profile.dis.rows() {
          let sum: f64 = row.sum();
          assert!(sum.is_finite(), "Dense node profile row sum is not finite: {sum}");
          pretty_assert_ulps_eq!(sum, 1.0, epsilon = 1e-6);
        }
      }
    }

    for node_data in sparse.node_states.values() {
      assert!(
        node_data.profile.log_lh.value().is_finite(),
        "Sparse node profile log_lh is not finite: {}",
        node_data.profile.log_lh.value()
      );

      for (pos, var_pos) in &node_data.profile.variable {
        let sum: f64 = aview1(var_pos.dis).sum();
        assert!(
          sum.is_finite(),
          "Sparse variable position {pos} sum is not finite: {sum}"
        );
        pretty_assert_ulps_eq!(sum, 1.0, epsilon = 1e-6);
      }

      for (char_key, fixed_dis) in node_data.profile.fixed.iter().map(|fixed| (fixed.state, &fixed.dis)) {
        let sum: f64 = aview1(fixed_dis).sum();
        assert!(
          sum.is_finite(),
          "Sparse fixed distribution for char {char_key:?} sum is not finite: {sum}"
        );
        pretty_assert_ulps_eq!(sum, 1.0, epsilon = 1e-6);
      }
    }

    Ok(())
  }

  #[test]
  fn test_marginal_dense_sparse_ambiguous_r_reference_state_consistency() -> Result<(), Report> {
    let aln = ambiguous_r_in_g_clade_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let gtr_dense = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let gtr_sparse = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let (log_lh_dense, dense_partition) = run_dense_marginal(&graph, &branch_lengths, &names, &aln, gtr_dense)?;
    let (log_lh_sparse, sparse_partition) = run_sparse_marginal(&graph, &branch_lengths, &names, &aln, gtr_sparse)?;

    pretty_assert_ulps_eq!(log_lh_dense, log_lh_sparse, epsilon = 1e-10);

    let dense_sequences = reconstruct_named_sequences(&graph, &names, &dense_partition)?;
    let sparse_sequences = reconstruct_named_sequences(&graph, &names, &sparse_partition)?;
    assert_eq!(dense_sequences, sparse_sequences);

    let dense_branch_subs =
      edge_subs_by_edge_name(&graph, &names, |key| dense_partition.view().edge_subs(&graph, key))?;
    let sparse_branch_subs =
      edge_subs_by_edge_name(&graph, &names, |key| sparse_partition.view().edge_subs(&graph, key))?;
    assert_eq!(dense_branch_subs, sparse_branch_subs);

    Ok(())
  }

  #[test]
  fn test_marginal_consistency_leaf_unknown_under_gapped_ancestor_stays_unknown() -> Result<(), Report> {
    let aln = fasta_read(
      indoc! {r#"
      >L1
      ACNNGT
      >L2
      ACNNGT
      >L3
      AC--GT
    "#}
      .as_bytes(),
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect_vec();
    let nwk_parsed = nwk_read(b"((L1:0.1,L2:0.1)X:0.1,L3:0.1)root:0.01;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph: Graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let x_key = find_node_key_by_name(&graph, &names, "X").expect("node X must exist");

    let (_, dense) = run_dense_marginal(&graph, &branch_lengths, &names, &aln, jc69(JC69Params::default())?)?;
    let (dense, _) = marginal_update_reconstruction(dense, &graph, &branch_lengths_or_zero(&branch_lengths))?;
    let (_, sparse) = run_sparse_marginal(&graph, &branch_lengths, &names, &aln, jc69(JC69Params::default())?)?;

    let expected = "ACNNGT";
    let dense_x = dense.node_sequence(&graph, false, x_key)?;
    let sparse_x = sparse.node_sequence(&graph, false, x_key)?;
    assert_eq!(
      expected,
      dense_x.to_string(),
      "X has only unknown children at columns 2-3, so a pass that reads observed leaf gaps keeps them unknown"
    );
    assert_eq!(
      expected,
      sparse_x.to_string(),
      "sparse reads observed leaf gaps from its single Fitch pass"
    );
    for (leaf, expected_leaf) in [("L1", "ACNNGT"), ("L2", "ACNNGT"), ("L3", "AC--GT")] {
      let leaf_key = find_node_key_by_name(&graph, &names, leaf).expect("leaf must exist");
      assert_eq!(
        expected_leaf,
        dense.node_sequence(&graph, false, leaf_key)?.to_string(),
        "dense leaf {leaf} keeps its observed unknown and gap states"
      );
      assert_eq!(
        expected_leaf,
        sparse.node_sequence(&graph, false, leaf_key)?.to_string(),
        "sparse leaf {leaf} keeps its observed unknown and gap states"
      );
    }
    Ok(())
  }

  fn reconstruct_named_sequences(
    graph: &Graph,
    names: &NodeNames,
    recon: &MarginalReconstruction,
  ) -> Result<BTreeMap<String, Seq>, Report> {
    let sampled = recon.sample_sequences(graph, SampleMode::Argmax, &mut get_random_number_generator(0))?;
    emitted_sequences_by_name(names, &internal_node_keys(graph), &sampled, |key| {
      recon.node_sequence(graph, false, key)
    })
  }

  fn edge_subs_by_edge_name(
    graph: &Graph,
    names: &NodeNames,
    edge_subs: impl Fn(GraphEdgeKey) -> Result<Vec<Sub>, Report>,
  ) -> Result<BTreeMap<String, Vec<Sub>>, Report> {
    graph
      .get_edges()
      .map(|edge_ref| {
        let edge = edge_ref;
        let parent_name = names[edge.source()].to_owned();
        let child_name = names[edge.target()].to_owned();
        let edge_name = format!("{parent_name}->{child_name}");
        let subs = edge_subs(edge.key())?;
        Ok((edge_name, subs))
      })
      .collect()
  }

  #[test]
  fn test_marginal_posteriors_sum_to_one_skewed_gtr() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let pi = array![0.9, 0.06, 0.02, 0.02];

    let W = array![
      [0.0, 1.0, 1.0, 1.0],
      [1.0, 0.0, 1.0, 1.0],
      [1.0, 1.0, 0.0, 1.0],
      [1.0, 1.0, 1.0, 0.0],
    ];

    let gtr = GTR::builder()
      .n_states(alphabet.n_canonical())
      .mu(Mu::FromMatrix(1.0))
      .W(W)
      .pi(pi)
      .build()?;

    let tree_newick = "((A:0.601,B:0.301):0.1,C:0.2):0.001;";
    let nwk_parsed = nwk_read(tree_newick.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let aln: Vec<AlignmentRecord> = fasta_read(
      indoc! {r#"
        >A
        AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT
        >B
        AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
        >C
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
      "#}
      .as_bytes(),
      &alphabet,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let partition = PartitionMarginalDense::new(alphabet, &graph, &leaf_seq_inputs(&graph, &names, aln))?;
    let recon = DenseReconstruction::seeded(partition, gtr);

    let (recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    for (node_key, node_data) in &recon.node_states {
      if node_data.profile.dis.is_empty() {
        continue;
      }
      for (pos, row) in node_data.profile.dis.rows().into_iter().enumerate() {
        let sum: f64 = row.sum();
        assert!(
          sum.is_finite(),
          "Node {node_key:?} position {pos}: posterior sum is not finite: {sum}"
        );
        pretty_assert_ulps_eq!(sum, 1.0, epsilon = 1e-6);
      }
    }

    Ok(())
  }

  #[test]
  fn test_marginal_sparse_uniform_site_rates_matches_scalar() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let seq_len = get_common_length(&aln)?;

    let gtr_scalar = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let (log_lh_scalar, _) = run_sparse_marginal(&graph, &branch_lengths, &names, &aln, gtr_scalar)?;

    let gtr_uniform = GTR::builder()
      .n_states(4)
      .mu(Mu::Normalized(1.0))
      .pi(Array1::ones(4))
      .site_rates(Array1::ones(seq_len))
      .build()?;
    let (log_lh_uniform, _) = run_sparse_marginal(&graph, &branch_lengths, &names, &aln, gtr_uniform)?;

    pretty_assert_ulps_eq!(log_lh_scalar, log_lh_uniform, epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_marginal_dense_uniform_site_rates_matches_scalar() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let seq_len = get_common_length(&aln)?;

    let gtr_scalar = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let (log_lh_scalar, _) = run_dense_marginal(&graph, &branch_lengths, &names, &aln, gtr_scalar)?;

    let gtr_uniform = GTR::builder()
      .n_states(4)
      .mu(Mu::Normalized(1.0))
      .pi(Array1::ones(4))
      .site_rates(Array1::ones(seq_len))
      .build()?;
    let (log_lh_uniform, _) = run_dense_marginal(&graph, &branch_lengths, &names, &aln, gtr_uniform)?;

    pretty_assert_ulps_eq!(log_lh_scalar, log_lh_uniform, epsilon = 1e-10);

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dense( run_dense_marginal)]
  #[case::sparse(run_sparse_marginal)]
  #[trace]
  fn test_marginal_site_rates_rejected_by_branch_length_contribution(#[case] #[notrace] run: MarginalRunner) -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let seq_len = get_common_length(&aln)?;

    let gtr = GTR::builder()
      .n_states(4)
      .mu(Mu::Normalized(1.0))
      .pi(Array1::ones(4))
      .site_rates(Array1::ones(seq_len))
      .build()?;
    let (_, reconstruction) = run(&graph, &branch_lengths, &names, &aln, gtr)?;
    let edge_key = graph.get_edges().next().expect("tree has an edge").key();

    assert_error!(
      reconstruction.view().create_edge_contribution(edge_key),
      "Branch length likelihood does not support per-site rates: the substitution model has site rates, but the branch length objective applies one rate to all sites"
    );

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dense_longer( run_dense_marginal,  17, "Number of site rates (17) differs from the number of sites in the partition (16)")]
  #[case::dense_single( run_dense_marginal,  1,  "Number of site rates (1) differs from the number of sites in the partition (16)")]
  #[case::sparse_longer(run_sparse_marginal, 17, "Number of site rates (17) differs from the number of sites in the partition (16)")]
  #[case::sparse_shorter(run_sparse_marginal, 15, "Number of site rates (15) differs from the number of sites in the partition (16)")]
  #[trace]
  fn test_marginal_site_rates_length_must_match_alignment(
    #[case] #[notrace] run: MarginalRunner,
    #[case] n_site_rates: usize,
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let gtr = GTR::builder()
      .n_states(4)
      .mu(Mu::Normalized(1.0))
      .pi(Array1::ones(4))
      .site_rates(Array1::ones(n_site_rates))
      .build()?;

    assert_error!(run(&graph, &branch_lengths, &names, &aln, gtr), expected);

    Ok(())
  }

  type MarginalRunner = fn(
    &Graph,
    &BTreeMap<GraphEdgeKey, Option<f64>>,
    &NodeNames,
    &[AlignmentRecord],
    GTR,
  ) -> Result<(f64, MarginalReconstruction), Report>;
}

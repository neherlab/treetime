#![allow(
  clippy::disallowed_methods,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::{ancestral_reconstruction, branch_lengths_or_zero};
  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
  use crate::ancestral::sample::SampleMode;
  use crate::ancestral::tip_states::TipStates;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::GTR;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::pretty_assert_ulps_eq;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use crate::seq::mutation::Sub;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use indoc::indoc;
  use treetime_graph::graph::Graph;

  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::LazyLock;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

  use treetime_utils::make_report;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  fn gap_free_alignment() -> Result<Vec<AlignmentRecord>, Report> {
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
        &*NUC_ALPHABET,
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect(),
    )
  }

  fn ambiguous_r_in_g_clade_alignment() -> Result<Vec<AlignmentRecord>, Report> {
    Ok(
      read_many_fasta_str(
        indoc! {r#"
      >A
      RCGTACGT
      >B
      GCGTACGT
      >C
      GCGTACGT
      >D
      GCGTACGT
    "#},
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
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    gtr: GTR,
  ) -> Result<(f64, DenseReconstruction), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let partition = PartitionMarginalDense::new(0, alphabet, get_common_length(aln)?);
    let node_states = partition.attach_sequences(graph, &node_seq_inputs(graph, names, aln.to_vec()))?;
    let recon = DenseReconstruction::seeded(partition, gtr, node_states);
    let (recon, log_lh) = recon.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    let log_lh = log_lh.value();
    Ok((log_lh, recon))
  }

  fn run_sparse_marginal(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    gtr: GTR,
  ) -> Result<(f64, SparseReconstruction), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let fitch = create_fitch_partition(graph, 0, alphabet, &node_seq_inputs(graph, names, aln.to_vec()))?;
    let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
    let recon = SparseReconstruction::seeded(partition, gtr, node_states);
    let (recon, log_lh) = recon.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    let log_lh = log_lh.value();
    Ok((log_lh, recon))
  }

  #[test]
  fn test_marginal_dense_sparse_log_lh_consistency_gap_free() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
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
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
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

    let dense = &dense_partition;
    let sparse = &sparse_partition;

    for node_key in [root_key, ab_key] {
      let dense_node = &dense.node_states[&node_key];
      let sparse_node = &sparse.node_states[&node_key];

      for (&pos, var_pos) in &sparse_node.profile.variable {
        let dense_row = dense_node.profile.dis.row(pos).to_owned();
        pretty_assert_ulps_eq!(dense_row, var_pos.dis.clone(), epsilon = 1e-6);
      }
    }

    Ok(())
  }

  #[test]
  fn test_marginal_dense_sparse_ambiguous_character_expectations_documented() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
      indoc! {r#"
      >A
      ACGTACGTACGTACNT
      >B
      ACGTACGTACGTACRA
      >C
      ACGTACGTACGTACGG
      >D
      ACGTACGTACGTACGC
    "#},
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
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

    let dense = &dense_partition;
    let sparse = &sparse_partition;

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
        let sum: f64 = var_pos.dis.sum();
        assert!(
          sum.is_finite(),
          "Sparse variable position {pos} sum is not finite: {sum}"
        );
        pretty_assert_ulps_eq!(sum, 1.0, epsilon = 1e-6);
      }

      for (char_key, fixed_dis) in &node_data.profile.fixed {
        let sum: f64 = fixed_dis.sum();
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
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
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

    let (log_lh_dense, mut dense_partition) = run_dense_marginal(&graph, &branch_lengths, &names, &aln, gtr_dense)?;
    let (log_lh_sparse, mut sparse_partition) = run_sparse_marginal(&graph, &branch_lengths, &names, &aln, gtr_sparse)?;

    pretty_assert_ulps_eq!(log_lh_dense, log_lh_sparse, epsilon = 1e-10);

    let dense_sequences = reconstruct_named_sequences_dense(&graph, &names, &mut dense_partition)?;
    let sparse_sequences = reconstruct_named_sequences_sparse(&graph, &names, &mut sparse_partition)?;
    assert_eq!(dense_sequences, sparse_sequences);

    let dense_branch_subs = edge_subs_by_edge_name(&graph, &names, |key| dense_partition.edge_subs(&graph, key))?;
    let sparse_branch_subs = edge_subs_by_edge_name(&graph, &names, |key| sparse_partition.edge_subs(key))?;
    assert_eq!(dense_branch_subs, sparse_branch_subs);

    Ok(())
  }

  fn reconstruct_named_sequences_dense(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    recon: &mut DenseReconstruction,
  ) -> Result<BTreeMap<String, String>, Report> {
    let mut actual = BTreeMap::new();
    let DenseReconstruction {
      partition, node_states, ..
    } = recon;
    let mut rng = rand::thread_rng();
    ancestral_reconstruction(graph, |node| {
      let Some(seq) = partition.reconstruct_node_sequence(
        node_states,
        node,
        TipStates {
          include_leaves: false,
          impute: false,
        },
        SampleMode::Argmax,
        &mut rng,
      ) else {
        return Ok(false);
      };
      actual.insert(
        names[&node.key].clone().expect("all test nodes are named"),
        seq.to_string(),
      );
      Ok(true)
    })?;
    Ok(actual)
  }

  fn reconstruct_named_sequences_sparse(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    recon: &mut SparseReconstruction,
  ) -> Result<BTreeMap<String, String>, Report> {
    let mut actual = BTreeMap::new();
    let SparseReconstruction {
      partition,
      node_states,
      edges,
      ..
    } = recon;
    let mut rng = rand::thread_rng();
    ancestral_reconstruction(graph, |node| {
      let Some(seq) = partition.reconstruct_node_sequence(
        node_states,
        &edges.forward,
        node,
        TipStates {
          include_leaves: false,
          impute: false,
        },
        SampleMode::Argmax,
        &mut rng,
      )?
      else {
        return Ok(false);
      };
      actual.insert(
        names[&node.key].clone().expect("all test nodes are named"),
        seq.to_string(),
      );
      Ok(true)
    })?;
    Ok(actual)
  }

  fn edge_subs_by_edge_name(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    edge_subs: impl Fn(GraphEdgeKey) -> Result<Vec<Sub>, Report>,
  ) -> Result<BTreeMap<String, Vec<Sub>>, Report> {
    graph
      .get_edges()
      .map(|edge_ref| {
        let edge = edge_ref;
        let parent_name = names.get(&edge.source()).cloned().flatten().expect("named parent");
        let child_name = names.get(&edge.target()).cloned().flatten().expect("named child");
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
      .mu(1.0)
      .W(W)
      .pi(pi)
      .build()?;

    let tree_newick = "((A:0.601,B:0.301):0.1,C:0.2):0.001;";
    let nwk_parsed = nwk_read_str(tree_newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
      indoc! {r#"
        >A
        AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT
        >B
        AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
        >C
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
      "#},
      &alphabet,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let partition = PartitionMarginalDense::new(0, alphabet, get_common_length(&aln)?);
    let node_states = partition.attach_sequences(&graph, &node_seq_inputs(&graph, &names, aln))?;
    let recon = DenseReconstruction::seeded(partition, gtr, node_states);

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
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let seq_len = get_common_length(&aln)?;

    let gtr_scalar = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let (log_lh_scalar, _) = run_sparse_marginal(&graph, &branch_lengths, &names, &aln, gtr_scalar)?;

    let mut gtr_uniform = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    gtr_uniform.site_rates = Some(Array1::ones(seq_len));
    let (log_lh_uniform, _) = run_sparse_marginal(&graph, &branch_lengths, &names, &aln, gtr_uniform)?;

    pretty_assert_ulps_eq!(log_lh_scalar, log_lh_uniform, epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_marginal_dense_uniform_site_rates_matches_scalar() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let seq_len = get_common_length(&aln)?;

    let gtr_scalar = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let (log_lh_scalar, _) = run_dense_marginal(&graph, &branch_lengths, &names, &aln, gtr_scalar)?;

    let mut gtr_uniform = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    gtr_uniform.site_rates = Some(Array1::ones(seq_len));
    let (log_lh_uniform, _) = run_dense_marginal(&graph, &branch_lengths, &names, &aln, gtr_uniform)?;

    pretty_assert_ulps_eq!(log_lh_scalar, log_lh_uniform, epsilon = 1e-10);

    Ok(())
  }
}

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
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::pretty_assert_ulps_eq;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;

  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;

  use ndarray::array;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::path::PathBuf;
  use std::sync::LazyLock;
  use treetime_io::fasta::{read_many_fasta_path, read_many_fasta_str};
  use treetime_io::nwk::{nwk_read_file, nwk_read_str};
  use treetime_primitives::AlignmentRecord;

  use treetime_utils::make_report;

  fn build_dense_recon(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    alphabet: Alphabet,
    index: usize,
    gtr: GTR,
  ) -> Result<DenseReconstruction, Report> {
    let partition = PartitionMarginalDense::new(index, alphabet, get_common_length(aln)?);
    let node_states = partition.attach_sequences(graph, &node_seq_inputs(graph, names, aln.to_vec()))?;
    let recon = DenseReconstruction::seeded(partition, gtr, node_states);
    let (recon, _) = recon.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    Ok(recon)
  }

  fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .map(PathBuf::from)
      .expect("project has workspace root")
  }

  #[test]
  fn test_root_sequence_matches_python_h3n2_na_20() -> Result<(), Report> {
    let root = project_root();
    let tree_path = root.join("data/flu/h3n2/20/tree.nwk");
    let aln_path = root.join("data/flu/h3n2/20/aln.fasta.xz");

    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln: Vec<AlignmentRecord> = read_many_fasta_path(&[aln_path], &alphabet)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();

    let nwk_parsed = nwk_read_file(&tree_path)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;

    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let mut recon = build_dense_recon(&graph, &branch_lengths, &names, &aln, alphabet, 0, gtr)?;

    let mut root_seq = String::new();
    {
      let DenseReconstruction {
        partition, node_states, ..
      } = &mut recon;
      let mut rng = rand::thread_rng();
      ancestral_reconstruction(&graph, |node| {
        let seq = partition.reconstruct_node_sequence(
          node_states,
          node,
          TipStates {
            include_leaves: false,
            impute: false,
          },
          SampleMode::Argmax,
          &mut rng,
        )?;
        if names[&node.key].as_deref() == Some("NODE_0000000") {
          root_seq = seq.to_string();
        }
        Some(())
      })?;
    }

    let expected = "ATGAATCCAAATCAAAAGATAATAACGATTGGCTCTGTTTCTCTCACCATTTCCACAATATGCTTCTTCATGCAAATTGCCATCTTGATAACTACTGTAACATTGCATTTCAAGCAATATGAATTCAACTCCCCCCCAAACAACCAAGTGATGCTGTGTGAACCAACAATAATAGAAAGAAACATAACAGAGATAGTGTATCTGACCAACACCACCATAGAGAAGGAAATATGCCCCAAACCAGCAGAATACAGAAATTGGTCAAAACCGCAATGTGGCATTACAGGATTTGCACCTTTCTCTAAGGACAATTCGATTAGGCTTTCCGCTGGTGGGGACATCTGGGTGACAAGAGAACCTTATGTGTCATGCGATCCTGACAAGTGTTATCAATTTGCCCTTGGACAGGGAACAACACTAAACAACGTGCATTCAAATAACACAGTACGTGATAGGACCCCTTATCGGACTCTATTGATGAATGAGTTGGGTGTTCCTTTTCATCTGGGGACCAAGCAAGTGTGCATAGCATGGTCCAGCTCAAGTTGTCACGATGGAAAAGCATGGCTGCATGTTTGTATAACGGGGGATGATAAAAATGCAACTGCTAGCTTCATTTACAATGGGAGGCTTGTAGATAGTGTTGTTTCATGGTCCAAAGAAATTCTCAGGACCCAGGAGTCAGAATGCGTTTGTATCAATGGAACTTGTACAGTAGTAATGACTGATGGAAGTGCTTCAGGAAAAGCTGATACTAAAATACTATTCATTGAGGAGGGGAAAATCGTTCATACTAGCACATTGTCAGGAAGTGCTCAGCATGTCGAAGAGTGCTCTTGCTATCCTCGATATCCTGGTGTCAGATGTGTCTGCAGAGACAACTGGAAAGGCTCCAATCGGCCCATCGTAGATATAAACATAAAGGATCATAGCATTGTTTCCAGTTATGTGTGTTCAGGACTTGTTGGAGACACACCCAGAAAAAACGACAGCTCCAGCAGTAGCCATTGTTTGGATCCTAACAATGAAGAAGGTGGTCATGGAGTGAAAGGCTGGGCCTTTGATGATGGAAATGACGTGTGGATGGGAAGAACAATCAACGAGACGTCACGCTTAGGGTATGAAACCTTCAAAGTCATTGAAGGCTGGTCCAACCCTAAGTCCAAATTGCAGATAAATAGGCAAGTCATAGTTGACAGAGGTGATAGGTCCGGTTATTCTGGTATTTTCTCTGTTGAAGGCAAAAGCTGCATCAATCGGTGCTTTTATGTGGAGTTGATTAGGGGAAGAAAAGAGGAAACTGAAGTCTTGTGGACCTCAAACAGTATTGTTGTGTTTTGTGGCACCTCAGGTACATATGGAACAGGCTCATGGCCTGATGGGGCGGACCTCAATCTCATGCCTATA";

    assert_eq!(expected, root_seq, "Root sequence mismatch with Python v0");

    Ok(())
  }

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  fn make_python_reference_gtr() -> Result<GTR, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    GTR::new(GTRParams {
      n_states,
      mu: 1.0,
      W: None,
      pi: array![0.2, 0.3, 0.15, 0.35],
    })
  }

  fn make_python_reference_gtr2() -> Result<GTR, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    GTR::new(GTRParams {
      n_states,
      mu: 1.0,
      W: None,
      pi: array![0.4, 0.15, 0.12, 0.33],
    })
  }

  const PYTHON_TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  const PYTHON_ALN: &str = ">A\nACATCGCCNNA--GAC\n>B\nGCATCCCTGTA-NG--\n>C\nCCGGCGATGTRTTG--\n>D\nTCGGCCGTGTRTTG--\n";

  #[test]
  fn test_internal_node_ab_profile_matches_python() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(PYTHON_TREE)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let gtr = make_python_reference_gtr()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let recon = build_dense_recon(&graph, &branch_lengths, &names, &aln, alphabet, 0, gtr)?;

    let ab_key = find_node_key_by_name(&graph, &names, "AB").ok_or_else(|| make_report!("Node AB not found"))?;
    let ab_profile = &recon.node_states[&ab_key].profile.dis;

    let pos0_profile = ab_profile.row(0);

    pretty_assert_ulps_eq!(pos0_profile[0], 0.51275208, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_profile[1], 0.09128506, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_profile[2], 0.24647255, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_profile[3], 0.14949031, epsilon = 1e-6);

    Ok(())
  }

  #[test]
  fn test_root_profile_matches_python() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(PYTHON_TREE)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let gtr = make_python_reference_gtr()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let recon = build_dense_recon(&graph, &branch_lengths, &names, &aln, alphabet, 0, gtr)?;

    let root_key = find_node_key_by_name(&graph, &names, "root").ok_or_else(|| make_report!("Node root not found"))?;
    let root_profile = &recon.node_states[&root_key].profile.dis;

    let pos0_profile = root_profile.row(0);

    pretty_assert_ulps_eq!(pos0_profile[0], 0.28212327, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_profile[1], 0.21643546, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_profile[2], 0.13800802, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_profile[3], 0.36343326, epsilon = 1e-6);

    Ok(())
  }

  #[test]
  fn test_internal_node_cd_profile_valid() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(PYTHON_TREE)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let gtr = make_python_reference_gtr()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let recon = build_dense_recon(&graph, &branch_lengths, &names, &aln, alphabet, 0, gtr)?;

    let cd_key = find_node_key_by_name(&graph, &names, "CD").ok_or_else(|| make_report!("Node CD not found"))?;
    let cd_profile = &recon.node_states[&cd_key].profile.dis;

    for (pos, row) in cd_profile.rows().into_iter().enumerate() {
      let sum: f64 = row.sum();
      assert!(
        (sum - 1.0).abs() < 1e-10,
        "CD profile row {pos} not normalized: sum={sum}"
      );
      for (state, &val) in row.iter().enumerate() {
        assert!(val >= 0.0, "CD profile row {pos} state {state} negative: {val}");
        assert!(val.is_finite(), "CD profile row {pos} state {state} not finite: {val}");
      }
    }

    Ok(())
  }

  #[test]
  fn test_all_internal_nodes_normalized() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(PYTHON_TREE)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let gtr = make_python_reference_gtr()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let recon = build_dense_recon(&graph, &branch_lengths, &names, &aln, alphabet, 0, gtr)?;

    for (key, node_data) in &recon.node_states {
      let profile = &node_data.profile.dis;
      for (pos, row) in profile.rows().into_iter().enumerate() {
        let sum: f64 = row.sum();
        assert!(
          (sum - 1.0).abs() < 1e-10,
          "Node {key:?} row {pos} not normalized: sum={sum}"
        );
      }
    }

    Ok(())
  }

  #[test]
  fn test_multi_partition_independent_computation() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(PYTHON_TREE)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();

    let gtr1 = make_python_reference_gtr()?;
    let gtr2 = make_python_reference_gtr2()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let recon1 = build_dense_recon(&graph, &branch_lengths, &names, &aln, alphabet.clone(), 0, gtr1)?;
    let recon2 = build_dense_recon(&graph, &branch_lengths, &names, &aln, alphabet, 1, gtr2)?;

    let root_key = find_node_key_by_name(&graph, &names, "root").ok_or_else(|| make_report!("Node root not found"))?;

    let root1 = &recon1.node_states[&root_key].profile.dis;
    let root2 = &recon2.node_states[&root_key].profile.dis;

    let pos0_p1 = root1.row(0);
    let pos0_p2 = root2.row(0);

    pretty_assert_ulps_eq!(pos0_p1[0], 0.28212327, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p1[1], 0.21643546, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p1[2], 0.13800802, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p1[3], 0.36343326, epsilon = 1e-6);

    pretty_assert_ulps_eq!(pos0_p2[0], 0.29664652, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p2[1], 0.20554302, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p2[2], 0.13582953, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_p2[3], 0.36198093, epsilon = 1e-6);

    Ok(())
  }

  #[test]
  fn test_multi_partition_internal_node_ab() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(PYTHON_TREE)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(PYTHON_ALN, &*NUC_ALPHABET)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();

    let gtr1 = make_python_reference_gtr()?;
    let gtr2 = make_python_reference_gtr2()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let recon1 = build_dense_recon(&graph, &branch_lengths, &names, &aln, alphabet.clone(), 0, gtr1)?;
    let recon2 = build_dense_recon(&graph, &branch_lengths, &names, &aln, alphabet, 1, gtr2)?;

    let ab_key = find_node_key_by_name(&graph, &names, "AB").ok_or_else(|| make_report!("Node AB not found"))?;

    let ab1 = &recon1.node_states[&ab_key].profile.dis;
    let pos0_ab1 = ab1.row(0);

    pretty_assert_ulps_eq!(pos0_ab1[0], 0.51275208, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab1[1], 0.09128506, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab1[2], 0.24647255, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab1[3], 0.14949031, epsilon = 1e-6);

    let ab2 = &recon2.node_states[&ab_key].profile.dis;
    let pos0_ab2 = ab2.row(0);

    pretty_assert_ulps_eq!(pos0_ab2[0], 0.52331521, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab2[1], 0.08336271, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab2[2], 0.24488808, epsilon = 1e-6);
    pretty_assert_ulps_eq!(pos0_ab2[3], 0.148434, epsilon = 1e-6);

    Ok(())
  }

  #[test]
  fn test_multi_partition_sparse_dense_consistency() -> Result<(), Report> {
    let simple_aln = ">A\nACATCGCCTTACGGAC\n>B\nGCATCCCTGTACTGAC\n>C\nCCGGCGATGTATTGAC\n>D\nTCGGCCGTGTATTGAC\n";

    let nwk_parsed = nwk_read_str(PYTHON_TREE)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(simple_aln, &*NUC_ALPHABET)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();

    let gtr = make_python_reference_gtr()?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let length = get_common_length(&aln)?;

    let dense_partition = PartitionMarginalDense::new(0, alphabet.clone(), length);
    let dense_node_states = dense_partition.attach_sequences(&graph, &node_seq_inputs(&graph, &names, aln.clone()))?;
    let dense_recon = DenseReconstruction::seeded(dense_partition, gtr.clone(), dense_node_states);
    let (dense_recon, dense_log_lh) = dense_recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let dense_log_lh = dense_log_lh.value();

    let fitch = create_fitch_partition(&graph, 0, alphabet, &node_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let sparse_recon = SparseReconstruction::seeded(partition, gtr, node_states);
    let (sparse_recon, sparse_log_lh) =
      sparse_recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let sparse_log_lh = sparse_log_lh.value();

    pretty_assert_ulps_eq!(dense_log_lh, sparse_log_lh, epsilon = 1e-10);

    Ok(())
  }
}

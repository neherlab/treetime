#[cfg(test)]
mod tests {
  //! Golden master tests for infer_gtr_dense.
  //!
  //! Compares Rust v1 `infer_gtr_dense` against Python v0 `TreeAnc.infer_gtr(marginal=True)`.
  //!
  //! Golden outputs captured via `gm_infer_gtr_dense_capture` script.

  use crate::seq::alignment::node_seq_inputs;

  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::DenseReconstruction;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::infer_gtr::common::{InferGtrOptions, InferGtrResult, infer_gtr_impl};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::shared::update::MarginalPasses;
  use crate::pretty_assert_ulps_eq;
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use lazy_static::lazy_static;
  use treetime_graph::graph::Graph;

  use rstest::rstest;
  use serde::Deserialize;
  use std::collections::BTreeMap;
  use std::fs;
  use std::path::{Path, PathBuf};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::fasta::{read_many_fasta_path, read_many_fasta_str};
  use treetime_io::nwk::{nwk_read_file, nwk_read_str};
  use treetime_primitives::AlignmentRecord;

  #[rstest]
  #[case::simple_4taxa("simple_4taxa")]
  #[case::star_topology("star_topology")]
  #[case::caterpillar("caterpillar")]
  #[case::multiple_mutations("multiple_mutations")]
  #[case::varying_branch_lengths("varying_branch_lengths")]
  #[case::large_branchy_uneven("large_branchy_uneven")]
  fn test_gm_infer_gtr_dense_synthetic(#[case] case_name: &str) -> Result<(), Report> {
    let case = &INPUTS.synthetic[case_name];
    let expected = &OUTPUTS.synthetic[case_name];

    let fasta_str = alignment_to_fasta(&case.alignment);
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(&fasta_str, &*NUC_ALPHABET)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let (graph, recon, branch_lengths) = setup_dense_partition(&case.tree, &aln)?;

    let counts = recon.partition.count_transitions(
      &recon.gtr,
      &graph,
      &branch_lengths_or_zero(&branch_lengths),
      &recon.node_states,
      &recon.edges.backward,
      &recon.edges.forward,
    )?;
    let actual = infer_gtr_impl(&counts, &InferGtrOptions::default())?;

    // Short synthetic sequences: limited floating-point accumulation, tight tolerance
    pretty_assert_ulps_eq!(&expected.W, &actual.W, epsilon = 1e-8);
    pretty_assert_ulps_eq!(&expected.pi, &actual.pi, epsilon = 1e-8);
    pretty_assert_ulps_eq!(expected.mu, actual.mu, epsilon = 1e-8);

    Ok(())
  }

  #[rstest]
  #[case::dengue_20("dengue_20")]
  #[case::ebola_20("ebola_20")]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  #[case::rsv_a_20("rsv_a_20")]
  #[case::lassa_L_50("lassa_L_50")]
  // #[case::mpox_clade_ii_20("mpox_clade_ii_20")] // Slow
  #[case::tb_20("tb_20")]
  fn test_gm_infer_gtr_dense_real(#[case] case_name: &str) -> Result<(), Report> {
    let case = &INPUTS.real[case_name];
    let expected = &OUTPUTS.real[case_name];

    let (graph, recon, branch_lengths) = setup_dense_partition_from_files(&case.tree_path, &case.alignment_path)?;
    let counts = recon.partition.count_transitions(
      &recon.gtr,
      &graph,
      &branch_lengths_or_zero(&branch_lengths),
      &recon.node_states,
      &recon.edges.backward,
      &recon.edges.forward,
    )?;
    let actual = infer_gtr_impl(&counts, &InferGtrOptions::default())?;

    // BLAS drift between NumPy and ndarray scales with sequence length. mpox_clade_ii_20
    // (~200k positions) shows max diff ~2.3e-7. Tightest passing: 1e-6.
    // Measured at commit bba8c177 by running all cases and collecting max abs diff.
    pretty_assert_ulps_eq!(&expected.W, &actual.W, epsilon = 1e-6);
    pretty_assert_ulps_eq!(&expected.pi, &actual.pi, epsilon = 1e-6);
    pretty_assert_ulps_eq!(expected.mu, actual.mu, epsilon = 1e-6);

    Ok(())
  }

  const FIXTURES_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/src/gtr/infer_gtr/__tests__/__fixtures__");

  #[derive(Debug, Deserialize)]
  struct Inputs {
    synthetic: BTreeMap<String, SyntheticCase>,
    real: BTreeMap<String, RealCase>,
  }

  #[derive(Debug, Deserialize)]
  struct SyntheticCase {
    tree: String,
    alignment: BTreeMap<String, String>,
  }

  #[derive(Debug, Deserialize)]
  struct RealCase {
    tree_path: String,
    alignment_path: String,
  }

  #[derive(Debug, Deserialize)]
  struct Outputs {
    synthetic: BTreeMap<String, InferGtrResult>,
    real: BTreeMap<String, InferGtrResult>,
  }

  lazy_static! {
    static ref INPUTS: Inputs = {
      let path = Path::new(FIXTURES_DIR).join("gm_infer_gtr_dense_inputs.json");
      let content = fs::read_to_string(&path).expect("Failed to read inputs JSON");
      serde_json::from_str(&content).expect("Failed to parse inputs JSON")
    };
    static ref OUTPUTS: Outputs = {
      let path = Path::new(FIXTURES_DIR).join("gm_infer_gtr_dense_outputs.json");
      let content = fs::read_to_string(&path).expect("Failed to read outputs JSON");
      serde_json::from_str(&content).expect("Failed to parse outputs JSON")
    };
    static ref NUC_ALPHABET: Alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    static ref PROJECT_ROOT: PathBuf = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .expect("Failed to find project root")
      .to_path_buf();
  }

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

  fn setup_dense_partition_from_files(
    tree_path: impl AsRef<Path>,
    alignment_path: impl AsRef<Path>,
  ) -> Result<(Graph, DenseReconstruction, BTreeMap<GraphEdgeKey, Option<f64>>), Report> {
    let tree_path = PROJECT_ROOT.join(tree_path);
    let alignment_path = PROJECT_ROOT.join(alignment_path);

    let nwk_parsed = nwk_read_file(&tree_path)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;
    let aln: Vec<AlignmentRecord> = read_many_fasta_path(&[&alignment_path], &*NUC_ALPHABET)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();

    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let partition = PartitionMarginalDense::new(0, NUC_ALPHABET.clone(), get_common_length(&aln)?);
    let node_states = partition.attach_sequences(&graph, &node_seq_inputs(&graph, &names, aln))?;
    let recon = DenseReconstruction::seeded(partition, gtr, node_states);
    let (recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    Ok((graph, recon, branch_lengths))
  }

  fn alignment_to_fasta(aln: &BTreeMap<String, String>) -> String {
    aln
      .iter()
      .map(|(name, seq)| format!(">{name}\n{seq}"))
      .collect::<Vec<_>>()
      .join("\n")
  }
}

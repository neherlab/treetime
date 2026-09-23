#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::gtr_inference::infer_gtr_fitch;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::DenseReconstruction;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::GTR;
  use crate::gtr::infer_gtr::common::{InferGtrOptions, InferGtrResult, infer_gtr_impl};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::shared::update::MarginalPasses;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;

  use eyre::Report;
  use std::sync::LazyLock;
  use treetime_graph::graph::Graph;

  use ndarray::{Array1, Array2};
  use rstest::rstest;

  use std::path::PathBuf;
  use treetime_io::fasta::read_many_fasta_path;
  use treetime_io::nwk::nwk_read_file;
  use treetime_primitives::AlignmentRecord;

  #[rustfmt::skip]
  #[rstest]
  #[case::flu_h3n2_20(      "data/flu/h3n2/20/tree.nwk",      "data/flu/h3n2/20/aln.fasta.xz")]
  #[case::ebola_20(          "data/ebola/20/tree.nwk",          "data/ebola/20/aln.fasta.xz")]
  #[case::rsv_a_20(          "data/rsv/a/20/tree.nwk",          "data/rsv/a/20/aln.fasta.xz")]
  #[case::dengue_20(         "data/dengue/20/tree.nwk",         "data/dengue/20/aln.fasta.xz")]
  #[case::tb_20(             "data/tb/20/tree.nwk",             "data/tb/20/aln.fasta.xz")]
  #[case::lassa_L_50(        "data/lassa/L/50/tree.nwk",        "data/lassa/L/50/aln.fasta.xz")]
  #[trace]
  fn test_contract_dense_sparse_real_gtr(
    #[case] tree_path: &str,
    #[case] alignment_path: &str,
  ) -> Result<(), Report> {
    let result = infer_both(tree_path, alignment_path)?;

    let pi_cos = cosine_similarity(&result.dense.pi, &result.sparse.pi);
    let w_rel = relative_frobenius(&result.dense.W, &result.sparse.W);
    let mu_rel = relative_diff(result.dense.mu, result.sparse.mu);

    assert!(
      pi_cos > 0.997,
      "pi cosine similarity too low: {pi_cos} (threshold 0.997)\n  dense:  {:?}\n  sparse: {:?}",
      result.dense.pi.as_slice().unwrap(),
      result.sparse.pi.as_slice().unwrap(),
    );

    assert!(
      w_rel < 0.21,
      "W relative Frobenius norm too large: {w_rel} (threshold 0.21)\n  dense W:\n{}\n  sparse W:\n{}",
      result.dense.W,
      result.sparse.W,
    );

    assert!(
      mu_rel < 0.23,
      "mu relative difference too large: {mu_rel} (threshold 0.23)\n  dense:  {}\n  sparse: {}",
      result.dense.mu,
      result.sparse.mu,
    );

    Ok(())
  }

  static DENSE_NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(|| Alphabet::new(AlphabetName::Nuc).unwrap());
  static SPARSE_NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(|| Alphabet::new(AlphabetName::Nuc).unwrap());
  static PROJECT_ROOT: LazyLock<PathBuf> = LazyLock::new(|| {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .expect("Failed to find project root")
      .to_path_buf()
  });

  struct DenseSparseGtr {
    dense: GTR,
    sparse: GTR,
  }

  fn infer_both(tree_path: &str, alignment_path: &str) -> Result<DenseSparseGtr, Report> {
    let tree_path = PROJECT_ROOT.join(tree_path);
    let alignment_path = PROJECT_ROOT.join(alignment_path);

    let aln: Vec<AlignmentRecord> = read_many_fasta_path(&[&alignment_path], &*DENSE_NUC_ALPHABET)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();

    let dense = {
      let nwk_parsed = nwk_read_file(&tree_path)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;
      let partition = PartitionMarginalDense::new(0, DENSE_NUC_ALPHABET.clone(), get_common_length(&aln)?);
      let node_states = partition.attach_sequences(&graph, &node_seq_inputs(&graph, &names, aln.clone()))?;
      let recon = DenseReconstruction::seeded(
        partition,
        jc69(JC69Params {
          alphabet: AlphabetName::Nuc,
          ..JC69Params::default()
        })?,
        node_states,
      );
      let (recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
      let counts = recon.partition.count_transitions(
        &recon.gtr,
        &graph,
        &branch_lengths_or_zero(&branch_lengths),
        &recon.node_states,
        &recon.edges.backward,
        &recon.edges.forward,
      )?;
      let InferGtrResult { W, pi, mu } = infer_gtr_impl(&counts, &InferGtrOptions::default())?;
      let n_states = recon.partition.alphabet.n_canonical();
      GTR::builder().n_states(n_states).mu(mu).W(W).pi(pi).build()?
    };

    let sparse = {
      let nwk_parsed = nwk_read_file(&tree_path)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;
      let fitch = create_fitch_partition(
        &graph,
        0,
        SPARSE_NUC_ALPHABET.clone(),
        &node_seq_inputs(&graph, &names, aln),
      )?;
      infer_gtr_fitch(&fitch, &graph, &branch_lengths_or_zero(&branch_lengths))?
    };

    Ok(DenseSparseGtr { dense, sparse })
  }

  fn cosine_similarity(a: &Array1<f64>, b: &Array1<f64>) -> f64 {
    let dot = a.dot(b);
    let norm_a = a.dot(a).sqrt();
    let norm_b = b.dot(b).sqrt();
    dot / (norm_a * norm_b)
  }

  fn relative_frobenius(a: &Array2<f64>, b: &Array2<f64>) -> f64 {
    let diff = a - b;
    let diff_norm = diff.iter().map(|x| x * x).sum::<f64>().sqrt();
    let a_norm = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    diff_norm / a_norm.max(1e-15)
  }

  fn relative_diff(a: f64, b: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs()).max(1e-15)
  }
}

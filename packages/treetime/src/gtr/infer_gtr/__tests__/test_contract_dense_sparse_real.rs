//! Dense-sparse GTR cross-validation on real datasets.
//!
//! Dense (marginal posterior fractional counts) and sparse (Fitch parsimony integer counts)
//! produce structurally different mutation counts. The inferred GTR models (W, pi, mu) should
//! agree qualitatively on real data.
//!
//! Measured differences (2025-03):
//!
//! | dataset          | pi_cosine    | W_rel_frob | mu_rel_diff |
//! |------------------|--------------|------------|-------------|
//! | flu_h3n2_20      | 0.999998594  | 0.0173     | 0.0040      |
//! | ebola_20         | 0.999999994  | 0.0180     | 0.0236      |
//! | rsv_a_20         | 0.999991113  | 0.0070     | 0.0229      |
//! | dengue_20        | 0.999786491  | 0.0080     | 0.0413      |
//! | tb_20            | 0.999970946  | 0.0109     | 0.0247      |
//! | lassa_L_50       | 0.998686724  | 0.1015     | 0.0196      |
//! | mpox_clade_ii_20 | 0.999999989  | 0.0858     | 0.1126      |
//!
//! Outlier analysis:
//! - lassa_L_50: highest W_rel_frob (0.1015), lowest pi_cosine (0.9987). 50-sequence dataset
//!   with more internal branches amplifies parsimony vs marginal divergence on W shape.
//! - mpox_clade_ii_20: highest mu_rel_diff (0.1126). Very long genome (~200k positions) with
//!   few mutations between closely related sequences - parsimony and marginal diverge on
//!   overall rate scale when mutation signal is sparse relative to genome size.
//!
//! Thresholds set at ~2x headroom over worst observed:
//! - pi cosine > 0.997   (worst 0.9987, distance from 1.0 doubled)
//! - W rel frob < 0.21   (worst 0.1015, doubled)
//! - mu rel diff < 0.23  (worst 0.1126, doubled)

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::GTR;
  use crate::gtr::infer_gtr::dense::infer_gtr_dense;
  use crate::gtr::infer_gtr::sparse::infer_gtr_sparse;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use ndarray::{Array1, Array2};
  use parking_lot::RwLock;
  use rstest::rstest;
  use std::path::PathBuf;
  use std::sync::Arc;
  use treetime_io::fasta::read_many_fasta;
  use treetime_io::nwk::nwk_read_file;

  #[rustfmt::skip]
  #[rstest]
  #[case::flu_h3n2_20(      "data/flu/h3n2/20/tree.nwk",      "data/flu/h3n2/20/aln.fasta.xz")]
  #[case::ebola_20(          "data/ebola/20/tree.nwk",          "data/ebola/20/aln.fasta.xz")]
  #[case::rsv_a_20(          "data/rsv/a/20/tree.nwk",          "data/rsv/a/20/aln.fasta.xz")]
  #[case::dengue_20(         "data/dengue/20/tree.nwk",         "data/dengue/20/aln.fasta.xz")]
  #[case::tb_20(             "data/tb/20/tree.nwk",             "data/tb/20/aln.fasta.xz")]
  #[case::lassa_L_50(        "data/lassa/L/50/tree.nwk",        "data/lassa/L/50/aln.fasta.xz")]
  #[case::mpox_clade_ii_20(  "data/mpox/clade-ii/20/tree.nwk",  "data/mpox/clade-ii/20/aln.fasta.xz")]
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

  lazy_static! {
    static ref DENSE_NUC_ALPHABET: Alphabet = Alphabet::new(AlphabetName::Nuc, true).unwrap();
    static ref SPARSE_NUC_ALPHABET: Alphabet = Alphabet::new(AlphabetName::Nuc, false).unwrap();
    static ref PROJECT_ROOT: PathBuf = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .expect("Failed to find project root")
      .to_path_buf();
  }

  struct DenseSparseGtr {
    dense: GTR,
    sparse: GTR,
  }

  fn infer_both(tree_path: &str, alignment_path: &str) -> Result<DenseSparseGtr, Report> {
    let tree_path = PROJECT_ROOT.join(tree_path);
    let alignment_path = PROJECT_ROOT.join(alignment_path);

    let aln = read_many_fasta(&[&alignment_path], &*DENSE_NUC_ALPHABET)?;

    let dense = {
      let graph: GraphAncestral = nwk_read_file(&tree_path)?;
      let partition = Arc::new(RwLock::new(PartitionMarginalDense {
        index: 0,
        gtr: jc69(JC69Params {
          alphabet: AlphabetName::Nuc,
          treat_gap_as_unknown: true,
          ..JC69Params::default()
        })?,
        alphabet: DENSE_NUC_ALPHABET.clone(),
        length: get_common_length(&aln)?,
        nodes: btreemap! {},
        edges: btreemap! {},
      }));
      initialize_marginal(&graph, std::slice::from_ref(&partition), &aln)?;
      infer_gtr_dense(&partition, &graph)?
    };

    let sparse = {
      let graph: GraphAncestral = nwk_read_file(&tree_path)?;
      let partition = Arc::new(RwLock::new(PartitionMarginalSparse {
        index: 0,
        gtr: jc69(JC69Params {
          alphabet: AlphabetName::Nuc,
          treat_gap_as_unknown: false,
          ..JC69Params::default()
        })?,
        alphabet: SPARSE_NUC_ALPHABET.clone(),
        length: get_common_length(&aln)?,
        nodes: btreemap! {},
        edges: btreemap! {},
      }));
      compress_sequences(&graph, std::slice::from_ref(&partition), &aln)?;
      update_marginal(&graph, std::slice::from_ref(&partition))?;
      infer_gtr_sparse(&partition, &graph)?
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

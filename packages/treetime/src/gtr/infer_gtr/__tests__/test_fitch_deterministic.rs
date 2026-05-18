#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::gtr_inference::infer_gtr_fitch;
  use crate::payload::ancestral::GraphAncestral;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use lazy_static::lazy_static;
  use rstest::rstest;
  use std::path::PathBuf;
  use treetime_io::fasta::read_many_fasta;
  use treetime_io::nwk::nwk_read_file;

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    static ref PROJECT_ROOT: PathBuf = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .expect("Failed to find project root")
      .to_path_buf();
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::flu_h3n2_20(      "data/flu/h3n2/20/tree.nwk",      "data/flu/h3n2/20/aln.fasta.xz")]
  #[case::ebola_20(          "data/ebola/20/tree.nwk",          "data/ebola/20/aln.fasta.xz")]
  #[case::rsv_a_20(          "data/rsv/a/20/tree.nwk",          "data/rsv/a/20/aln.fasta.xz")]
  #[case::dengue_20(         "data/dengue/20/tree.nwk",         "data/dengue/20/aln.fasta.xz")]
  #[case::tb_20(             "data/tb/20/tree.nwk",             "data/tb/20/aln.fasta.xz")]
  #[trace]
  fn test_fitch_gtr_deterministic(
    #[case] tree_path: &str,
    #[case] alignment_path: &str,
  ) -> Result<(), Report> {
    let tree_path = PROJECT_ROOT.join(tree_path);
    let alignment_path = PROJECT_ROOT.join(alignment_path);
    let aln = read_many_fasta(&[&alignment_path], &*NUC_ALPHABET)?;

    let gtr_a = {
      let graph: GraphAncestral = nwk_read_file(&tree_path)?;
      let fitch = create_fitch_partition(&graph, 0, NUC_ALPHABET.clone(), &aln)?;
      infer_gtr_fitch(&fitch, &graph)?
    };

    let gtr_b = {
      let graph: GraphAncestral = nwk_read_file(&tree_path)?;
      let fitch = create_fitch_partition(&graph, 0, NUC_ALPHABET.clone(), &aln)?;
      infer_gtr_fitch(&fitch, &graph)?
    };

    pretty_assert_ulps_eq!(gtr_a.mu, gtr_b.mu, epsilon = 1e-15);
    pretty_assert_ulps_eq!(gtr_a.pi, gtr_b.pi, epsilon = 1e-15);
    pretty_assert_ulps_eq!(gtr_a.W, gtr_b.W, epsilon = 1e-15);
    Ok(())
  }
}

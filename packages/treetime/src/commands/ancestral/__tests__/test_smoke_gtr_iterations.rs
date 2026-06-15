#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::params::MethodAncestral;
  use crate::ancestral::pipeline::{AncestralInput, AncestralParams};
  use crate::ancestral::sample::SampleMode;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::progress::NoopProgress;
  use eyre::Report;
  use lazy_static::lazy_static;
  use std::path::PathBuf;
  use treetime_io::fasta::read_many_fasta;
  use treetime_io::nwk::nwk_read_file;

  lazy_static! {
    static ref PROJECT_ROOT: PathBuf = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .expect("Failed to find project root")
      .to_path_buf();
  }

  #[test]
  fn test_smoke_ancestral_gtr_iterations_sparse() -> Result<(), Report> {
    let alphabet = Alphabet::default();
    let graph = nwk_read_file(PROJECT_ROOT.join("data/flu/h3n2/20/tree.nwk"))?;
    let sequences = read_many_fasta(&[PROJECT_ROOT.join("data/flu/h3n2/20/aln.fasta.xz")], &alphabet)?;

    let params = AncestralParams {
      method: MethodAncestral::Marginal,
      model: GtrModelName::Infer,
      dense: None,
      reconstruct_tip_states: false,
      gtr_iterations: 3,
      site_specific_gtr: false,
      seed: None,
      sample_from_profile: SampleMode::Argmax,
      ignore_missing_alns: false,
    };

    let input = AncestralInput {
      graph,
      alphabet,
      sequences,
    };

    let result = crate::ancestral::pipeline::run(&params, input, |_, _| Ok(()), &NoopProgress)?;

    let gtr = result.output.gtr.expect("GTR should be fitted with --model=infer");
    assert!(
      gtr.mu > 0.0,
      "mu should be positive after GTR iterations, got {}",
      gtr.mu
    );

    Ok(())
  }

  #[test]
  fn test_smoke_ancestral_gtr_iterations_dense() -> Result<(), Report> {
    let alphabet = Alphabet::default();
    let graph = nwk_read_file(PROJECT_ROOT.join("data/flu/h3n2/20/tree.nwk"))?;
    let sequences = read_many_fasta(&[PROJECT_ROOT.join("data/flu/h3n2/20/aln.fasta.xz")], &alphabet)?;

    let params = AncestralParams {
      method: MethodAncestral::Marginal,
      model: GtrModelName::Infer,
      dense: Some(true),
      reconstruct_tip_states: false,
      gtr_iterations: 3,
      site_specific_gtr: false,
      seed: None,
      sample_from_profile: SampleMode::Argmax,
      ignore_missing_alns: false,
    };

    let input = AncestralInput {
      graph,
      alphabet,
      sequences,
    };

    let result = crate::ancestral::pipeline::run(&params, input, |_, _| Ok(()), &NoopProgress)?;

    let gtr = result.output.gtr.expect("GTR should be fitted with --model=infer");
    assert!(
      gtr.mu > 0.0,
      "mu should be positive after GTR iterations, got {}",
      gtr.mu
    );

    Ok(())
  }
}

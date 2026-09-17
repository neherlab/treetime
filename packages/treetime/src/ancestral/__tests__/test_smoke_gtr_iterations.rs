#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::attach::complete_alignment_for_leaves;
  use crate::ancestral::mask::create_mask;
  use crate::ancestral::params::MethodAncestral;
  use crate::ancestral::pipeline::AncestralParams;
  use crate::ancestral::sample::SampleMode;
  use crate::cancel::NoopCancel;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::progress::NoopProgress;
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use lazy_static::lazy_static;
  use std::path::PathBuf;
  use treetime_io::fasta::read_many_fasta_path;
  use treetime_io::nwk::{NwkFastaInput, nwk_read_file};

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
    let parse = nwk_read_file(PROJECT_ROOT.join("data/flu/h3n2/20/tree.nwk"))?;
    let sequences = read_many_fasta_path(&[PROJECT_ROOT.join("data/flu/h3n2/20/aln.fasta.xz")], &alphabet)?;

    let params = AncestralParams {
      method: MethodAncestral::Marginal,
      model: GtrModelName::Infer,
      dense: None,
      include_leaves: false,
      impute_missing_data: false,
      gtr_iterations: 3,
      site_specific_gtr: false,
      seed: None,
      sample_from_profile: SampleMode::Argmax,
      ignore_missing_alns: false,
    };
    let names = parse.names();
    let sequences = complete_alignment_for_leaves(&parse.graph, sequences, &alphabet, false, &names)?;
    let alignment_length = get_common_length(&sequences)?;
    let mask = create_mask(&sequences, alignment_length, &alphabet);
    let input = NwkFastaInput::from_parse_and_aln(parse, sequences);

    let result = crate::ancestral::pipeline::run(&params, &input, alphabet, mask, &NoopCancel, &NoopProgress)?;

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
    let parse = nwk_read_file(PROJECT_ROOT.join("data/flu/h3n2/20/tree.nwk"))?;
    let sequences = read_many_fasta_path(&[PROJECT_ROOT.join("data/flu/h3n2/20/aln.fasta.xz")], &alphabet)?;

    let params = AncestralParams {
      method: MethodAncestral::Marginal,
      model: GtrModelName::Infer,
      dense: Some(true),
      include_leaves: false,
      impute_missing_data: false,
      gtr_iterations: 3,
      site_specific_gtr: false,
      seed: None,
      sample_from_profile: SampleMode::Argmax,
      ignore_missing_alns: false,
    };
    let names = parse.names();
    let sequences = complete_alignment_for_leaves(&parse.graph, sequences, &alphabet, false, &names)?;
    let alignment_length = get_common_length(&sequences)?;
    let mask = create_mask(&sequences, alignment_length, &alphabet);
    let input = NwkFastaInput::from_parse_and_aln(parse, sequences);

    let result = crate::ancestral::pipeline::run(&params, &input, alphabet, mask, &NoopCancel, &NoopProgress)?;

    let gtr = result.output.gtr.expect("GTR should be fitted with --model=infer");
    assert!(
      gtr.mu > 0.0,
      "mu should be positive after GTR iterations, got {}",
      gtr.mu
    );

    Ok(())
  }
}

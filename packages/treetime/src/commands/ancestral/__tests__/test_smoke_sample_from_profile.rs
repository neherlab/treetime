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
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::path::PathBuf;
  use treetime_graph::node::Named;
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
  fn test_smoke_ancestral_sample_from_profile_root_reproducible() -> Result<(), Report> {
    let seqs_a = helpers::run_root_sampled(42)?;
    let seqs_b = helpers::run_root_sampled(42)?;
    assert_eq!(seqs_a, seqs_b);
    Ok(())
  }

  #[test]
  fn test_sample_from_profile_rejected_for_parsimony() {
    let alphabet = Alphabet::default();
    let graph = nwk_read_file(PROJECT_ROOT.join("data/flu/h3n2/20/tree.nwk")).unwrap();
    let sequences = read_many_fasta(&[PROJECT_ROOT.join("data/flu/h3n2/20/aln.fasta.xz")], &alphabet).unwrap();

    let params = AncestralParams {
      method: MethodAncestral::Parsimony,
      model: GtrModelName::Infer,
      dense: None,
      reconstruct_tip_states: false,
      gtr_iterations: 0,
      site_specific_gtr: false,
      seed: None,
      sample_from_profile: SampleMode::Root,
      ignore_missing_alns: false,
    };

    let input = AncestralInput {
      graph,
      alphabet,
      sequences,
    };

    let result = crate::ancestral::pipeline::run(&params, input, |_, _| Ok(()), &NoopProgress);
    assert!(result.is_err(), "parsimony with posterior sampling must be rejected");
    let err = result.err().unwrap().to_string();
    assert!(
      err.contains("requires --method-anc=marginal"),
      "expected sampling/method guard error, got: {err}"
    );
  }

  #[test]
  fn test_smoke_ancestral_sample_from_profile_all() -> Result<(), Report> {
    let seqs = helpers::run_sampled(SampleMode::All, 7)?;
    assert!(!seqs.is_empty(), "sampled sequences should not be empty");
    Ok(())
  }

  mod helpers {
    use super::*;

    pub fn run_root_sampled(seed: u64) -> Result<BTreeMap<String, String>, Report> {
      run_sampled(SampleMode::Root, seed)
    }

    pub fn run_sampled(mode: SampleMode, seed: u64) -> Result<BTreeMap<String, String>, Report> {
      let alphabet = Alphabet::default();
      let graph = nwk_read_file(PROJECT_ROOT.join("data/flu/h3n2/20/tree.nwk"))?;
      let sequences = read_many_fasta(&[PROJECT_ROOT.join("data/flu/h3n2/20/aln.fasta.xz")], &alphabet)?;

      let params = AncestralParams {
        method: MethodAncestral::Marginal,
        model: GtrModelName::Infer,
        dense: None,
        reconstruct_tip_states: false,
        gtr_iterations: 0,
        site_specific_gtr: false,
        seed: Some(seed),
        sample_from_profile: mode,
        ignore_missing_alns: false,
      };

      let input = AncestralInput {
        graph,
        alphabet,
        sequences,
      };

      let mut captured: BTreeMap<String, String> = BTreeMap::new();
      crate::ancestral::pipeline::run(
        &params,
        input,
        |node, seq| {
          let name = node.name().map_or_else(String::new, |n| n.as_ref().to_owned());
          captured.insert(name, seq.to_string());
          Ok(())
        },
        &NoopProgress,
      )?;

      Ok(captured)
    }
  }
}

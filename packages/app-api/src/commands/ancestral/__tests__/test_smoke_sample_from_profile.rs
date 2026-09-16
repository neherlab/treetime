#[cfg(test)]
mod tests {
  use treetime::alphabet::alphabet::Alphabet;
  use treetime::ancestral::attach::complete_alignment_for_leaves;
  use treetime::ancestral::mask::create_mask;
  use treetime::ancestral::params::MethodAncestral;
  use treetime::ancestral::pipeline::AncestralParams;
  use treetime::ancestral::sample::SampleMode;
  use treetime::gtr::get_gtr::GtrModelName;
  use treetime::progress::NoopProgress;
  use treetime::seq::alignment::get_common_length;
  use eyre::Report;
  use lazy_static::lazy_static;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
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
  fn test_smoke_ancestral_sample_from_profile_root_reproducible() -> Result<(), Report> {
    let seqs_a = helpers::run_root_sampled(42)?;
    let seqs_b = helpers::run_root_sampled(42)?;
    assert_eq!(seqs_a, seqs_b);
    Ok(())
  }

  #[test]
  fn test_sample_from_profile_rejected_for_parsimony() {
    let alphabet = Alphabet::default();
    let parse = nwk_read_file(PROJECT_ROOT.join("data/flu/h3n2/20/tree.nwk")).unwrap();
    let sequences = read_many_fasta_path(&[PROJECT_ROOT.join("data/flu/h3n2/20/aln.fasta.xz")], &alphabet).unwrap();

    let params = AncestralParams {
      method: MethodAncestral::Parsimony,
      model: GtrModelName::Infer,
      dense: None,
      include_leaves: false,
      impute_missing_data: false,
      gtr_iterations: 0,
      site_specific_gtr: false,
      seed: None,
      sample_from_profile: SampleMode::Root,
      ignore_missing_alns: false,
    };
    let names = parse.names();
    let sequences = complete_alignment_for_leaves(&parse.graph, sequences, &alphabet, false, &names).unwrap();
    let alignment_length = get_common_length(&sequences).unwrap();
    let mask = create_mask(&sequences, alignment_length, &alphabet);
    let input = NwkFastaInput::from_parse_and_aln(parse, sequences);

    let result = treetime::ancestral::pipeline::run(&params, &input, alphabet, mask, &NoopProgress);
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
      let parse = nwk_read_file(PROJECT_ROOT.join("data/flu/h3n2/20/tree.nwk"))?;
      let sequences = read_many_fasta_path(&[PROJECT_ROOT.join("data/flu/h3n2/20/aln.fasta.xz")], &alphabet)?;

      let params = AncestralParams {
        method: MethodAncestral::Marginal,
        model: GtrModelName::Infer,
        dense: None,
        include_leaves: false,
        impute_missing_data: false,
        gtr_iterations: 0,
        site_specific_gtr: false,
        seed: Some(seed),
        sample_from_profile: mode,
        ignore_missing_alns: false,
      };
      let names = parse.names();
      let sequences = complete_alignment_for_leaves(&parse.graph, sequences, &alphabet, false, &names)?;
      let alignment_length = get_common_length(&sequences)?;
      let mask = create_mask(&sequences, alignment_length, &alphabet);
      let input = NwkFastaInput::from_parse_and_aln(parse, sequences);

      // Read the reconstructed sequences back off the partition in the walk's emission order, as the
      // reconstructed-FASTA writer does. This exercises the same path the CLI streams to file.
      let result = treetime::ancestral::pipeline::run(&params, &input, alphabet, mask, &NoopProgress)?;
      let partition = result.partition.expect("marginal reconstruction produces a partition");
      let captured: BTreeMap<String, String> = result
        .output
        .emitted_nodes
        .iter()
        .map(|&key| {
          let name = input.nodes[&key].name.clone().unwrap_or_default();
          (name, partition.augur_node_sequence(key).to_string())
        })
        .collect();

      Ok(captured)
    }
  }
}

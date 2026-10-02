#[cfg(test)]
mod tests {
  use eyre::Report;
  use treetime_utils::io::json::{JsonPretty, json_write_str};

  #[test]
  fn test_smoke_ancestral_gtr_iterations_sparse() -> Result<(), Report> {
    let initial = helpers::run_flu_20(None, 0)?;
    let refined = helpers::run_flu_20(None, 3)?;

    assert!(
      refined.mu > 0.0,
      "mu should be positive after GTR iterations, got {}",
      refined.mu
    );
    assert_ne!(
      json_write_str(&initial, JsonPretty(false))?,
      json_write_str(&refined, JsonPretty(false))?,
      "GTR iterations must replace the parsimony-inferred GTR"
    );

    Ok(())
  }

  #[test]
  fn test_smoke_ancestral_gtr_iterations_dense() -> Result<(), Report> {
    let initial = helpers::run_flu_20(Some(true), 0)?;
    let refined = helpers::run_flu_20(Some(true), 3)?;

    assert!(
      refined.mu > 0.0,
      "mu should be positive after GTR iterations, got {}",
      refined.mu
    );
    assert_ne!(
      json_write_str(&initial, JsonPretty(false))?,
      json_write_str(&refined, JsonPretty(false))?,
      "GTR iterations must replace the parsimony-inferred GTR"
    );

    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::ancestral::attach::complete_alignment_for_leaves;
    use crate::ancestral::mask::create_mask;
    use crate::ancestral::params::{AncestralParams, MethodAncestral};
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::GtrModelName;
    use crate::gtr::gtr::GTR;
    use crate::partition::marginal::sample::SampleMode;
    use crate::progress::NoopProgress;
    use crate::seq::alignment::{AncestralInput, EdgeSeqInput, get_common_length, node_seq_inputs};
    use eyre::{OptionExt, Report};
    use std::path::PathBuf;
    use std::sync::LazyLock;
    use treetime_io::fasta::read_many_fasta_path;
    use treetime_io::nwk::nwk_read_file;
    use treetime_primitives::AlignmentRecord;

    static PROJECT_ROOT: LazyLock<PathBuf> = LazyLock::new(|| {
      PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("Failed to find project root")
        .to_path_buf()
    });

    pub(super) fn run_flu_20(dense: Option<bool>, gtr_iterations: usize) -> Result<GTR, Report> {
      let alphabet = Alphabet::default();
      let parse = nwk_read_file(PROJECT_ROOT.join("data/flu/h3n2/20/tree.nwk"))?;
      let sequences: Vec<AlignmentRecord> =
        read_many_fasta_path(&[PROJECT_ROOT.join("data/flu/h3n2/20/aln.fasta.xz")], &alphabet)?
          .into_iter()
          .map(AlignmentRecord::from)
          .collect();

      let params = AncestralParams {
        method: MethodAncestral::Marginal,
        model: GtrModelName::Infer,
        dense,
        include_leaves: false,
        impute_missing_data: false,
        gtr_iterations,
        site_specific_gtr: false,
        seed: None,
        sample_from_profile: SampleMode::Argmax,
      };
      let names = parse.names();
      let sequences = complete_alignment_for_leaves(&parse.graph, sequences, &alphabet, false, &names, &NoopProgress)?;
      let alignment_length = get_common_length(&sequences)?;
      let mask = create_mask(&sequences, alignment_length, &alphabet);
      let input = AncestralInput {
        nodes: node_seq_inputs(&parse.graph, &names, sequences),
        edges: parse
          .branch_lengths
          .into_iter()
          .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
          .collect(),
        graph: parse.graph,
      };

      let result = crate::ancestral::pipeline::run(
        &params,
        &input,
        alphabet,
        mask,
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?;
      result.output.gtr.ok_or_eyre("GTR should be fitted with --model=infer")
    }
  }
}

#[cfg(test)]
mod tests {
  use crate::ancestral::params::MethodAncestral;
  use crate::ancestral::sample::SampleMode;
  use crate::commands::ancestral::args::TreetimeAncestralArgs;
  use crate::commands::ancestral::run::run_ancestral_reconstruction;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::model::ModelArgs;
  use crate::commands::shared::output::OutputCoreArgs;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::progress::NoopProgress;
  use eyre::Report;
  use lazy_static::lazy_static;
  use pretty_assertions::assert_eq;
  use std::path::PathBuf;

  lazy_static! {
    static ref PROJECT_ROOT: PathBuf = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .expect("Failed to find project root")
      .to_path_buf();
  }

  /// Root sampling runs end-to-end through the ancestral command and is reproducible: two runs with
  /// the same seed produce byte-identical reconstructed sequences.
  #[test]
  fn test_smoke_ancestral_sample_from_profile_root_reproducible() -> Result<(), Report> {
    let seqs_a = helpers::run_root_sampled("tmp/test-sample-root-a", 42)?;
    let seqs_b = helpers::run_root_sampled("tmp/test-sample-root-b", 42)?;
    assert_eq!(seqs_a, seqs_b);
    Ok(())
  }

  mod helpers {
    use super::PROJECT_ROOT;
    use crate::ancestral::sample::SampleMode;
    use crate::commands::ancestral::args::TreetimeAncestralArgs;
    use crate::commands::ancestral::run::run_ancestral_reconstruction;
    use crate::commands::shared::alignment::AlignmentArgs;
    use crate::commands::shared::model::ModelArgs;
    use crate::commands::shared::output::OutputCoreArgs;
    use crate::gtr::get_gtr::GtrModelName;
    use crate::progress::NoopProgress;
    use eyre::Report;

    /// Run the ancestral command with root sampling at the given seed, returning the reconstructed
    /// `annotated_tree.reconstructed-nuc.fasta` contents for comparison across runs.
    pub fn run_root_sampled(out_subdir: &str, seed: u64) -> Result<String, Report> {
      let outdir = PROJECT_ROOT.join(out_subdir);
      std::fs::create_dir_all(&outdir)?;

      let args = TreetimeAncestralArgs {
        alignment: AlignmentArgs {
          alignment: vec![PROJECT_ROOT.join("data/flu/h3n2/20/aln.fasta.xz")],
        },
        tree: PROJECT_ROOT.join("data/flu/h3n2/20/tree.nwk"),
        model_args: ModelArgs {
          model: GtrModelName::Infer,
          ..ModelArgs::default()
        },
        sample_from_profile: SampleMode::Root,
        seed: Some(seed),
        output: OutputCoreArgs {
          output_all: Some(outdir.clone()),
          ..Default::default()
        },
        ..TreetimeAncestralArgs::default()
      };

      run_ancestral_reconstruction(&args, &NoopProgress)?;
      Ok(std::fs::read_to_string(
        outdir.join("annotated_tree.reconstructed-nuc.fasta"),
      )?)
    }
  }

  /// Posterior sampling with a non-marginal method is rejected by the pipeline.
  #[test]
  fn test_sample_from_profile_rejected_for_parsimony() {
    let args = TreetimeAncestralArgs {
      alignment: AlignmentArgs {
        alignment: vec![PROJECT_ROOT.join("data/flu/h3n2/20/aln.fasta.xz")],
      },
      tree: PROJECT_ROOT.join("data/flu/h3n2/20/tree.nwk"),
      method_anc: MethodAncestral::Parsimony,
      sample_from_profile: SampleMode::Root,
      output: OutputCoreArgs {
        output_all: Some(PROJECT_ROOT.join("tmp/test-sample-parsimony-reject")),
        ..Default::default()
      },
      ..TreetimeAncestralArgs::default()
    };

    let err = run_ancestral_reconstruction(&args, &NoopProgress)
      .expect_err("parsimony with posterior sampling must be rejected")
      .to_string();
    assert!(
      err.contains("requires --method-anc=marginal"),
      "expected sampling/method guard error, got: {err}"
    );
  }

  /// Sampling every node runs end-to-end without error.
  #[test]
  fn test_smoke_ancestral_sample_from_profile_all() -> Result<(), Report> {
    let outdir = PROJECT_ROOT.join("tmp/test-sample-all");
    std::fs::create_dir_all(&outdir)?;

    let args = TreetimeAncestralArgs {
      alignment: AlignmentArgs {
        alignment: vec![PROJECT_ROOT.join("data/flu/h3n2/20/aln.fasta.xz")],
      },
      tree: PROJECT_ROOT.join("data/flu/h3n2/20/tree.nwk"),
      model_args: ModelArgs {
        model: GtrModelName::Infer,
        ..ModelArgs::default()
      },
      sample_from_profile: SampleMode::All,
      seed: Some(7),
      output: OutputCoreArgs {
        output_all: Some(outdir.clone()),
        ..Default::default()
      },
      ..TreetimeAncestralArgs::default()
    };

    run_ancestral_reconstruction(&args, &NoopProgress)?;
    assert!(outdir.join("annotated_tree.reconstructed-nuc.fasta").exists());
    Ok(())
  }
}

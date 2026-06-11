#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::TreetimeAncestralArgs;
  use crate::commands::ancestral::run::run_ancestral_reconstruction;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::model::ModelArgs;
  use crate::commands::shared::output::OutputCoreArgs;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::progress::NoopProgress;
  use eyre::Report;
  use lazy_static::lazy_static;
  use std::path::PathBuf;

  lazy_static! {
    static ref PROJECT_ROOT: PathBuf = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .expect("Failed to find project root")
      .to_path_buf();
  }

  #[test]
  fn test_smoke_ancestral_gtr_iterations_sparse() -> Result<(), Report> {
    let outdir = PROJECT_ROOT.join("tmp/test-gtr-iter-sparse");
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
      gtr_iterations: 3,
      output: OutputCoreArgs {
        output_all: Some(outdir.clone()),
        ..Default::default()
      },
      ..TreetimeAncestralArgs::default()
    };

    let result = run_ancestral_reconstruction(&args, &NoopProgress)?;

    assert!(result.gtr.is_some());
    let gtr = result.gtr.unwrap();
    assert!(gtr.mu > 0.0, "mu should be positive after GTR iterations");

    assert!(outdir.join("annotated_tree.nwk").exists());
    assert!(outdir.join("annotated_tree.reconstructed-nuc.fasta").exists());

    Ok(())
  }

  #[test]
  fn test_smoke_ancestral_gtr_iterations_dense() -> Result<(), Report> {
    let outdir = PROJECT_ROOT.join("tmp/test-gtr-iter-dense");
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
      dense: Some(true),
      gtr_iterations: 3,
      output: OutputCoreArgs {
        output_all: Some(outdir.clone()),
        ..Default::default()
      },
      ..TreetimeAncestralArgs::default()
    };

    let result = run_ancestral_reconstruction(&args, &NoopProgress)?;

    assert!(result.gtr.is_some());
    let gtr = result.gtr.unwrap();
    assert!(gtr.mu > 0.0, "mu should be positive after GTR iterations");

    assert!(outdir.join("annotated_tree.nwk").exists());
    assert!(outdir.join("annotated_tree.reconstructed-nuc.fasta").exists());

    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::ancestral::run::run_ancestral_reconstruction;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::model::{GtrModelNameCli, ModelArgs};
  use eyre::{Report, WrapErr};
  use pretty_assertions::assert_eq;
  use std::fs;
  use std::path::{Path, PathBuf};
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_utils::{assert_error, vec_of_owned};

  const TREE: &str = "((A:0.1,B:0.1)AB:0.1,(C:0.1,D:0.1)CD:0.1)root;";
  const ALIGNMENT: &str = ">A\nATGACG\n>B\nACGACG\n>C\nACGACG\n>D\nACGACG\n";
  const TRANSLATION: &str = ">A\nMRL\n>B\nMKL\n>C\nMKL\n>D\nMKL\n";

  #[test]
  fn test_reconstructed_aa_fasta_written_to_the_expanded_file_name() -> Result<(), Report> {
    let dir = tempdir().wrap_err("When creating a temporary directory")?;
    let template = dir.path().join("aa.%GENE.fasta");
    let args = helpers::args(dir.path(), template.clone(), None)?;

    run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;

    assert_eq!(
      (true, false),
      (dir.path().join("aa.S.fasta").is_file(), template.exists())
    );
    Ok(())
  }

  #[test]
  fn test_reconstructed_aa_fasta_path_of_a_cds_shared_with_another_output_fails_before_the_run() -> Result<(), Report> {
    let dir = tempdir().wrap_err("When creating a temporary directory")?;
    let gtr = dir.path().join("S.json");
    let nuc_fasta = dir.path().join("nuc.fasta");
    let mut args = helpers::args(dir.path(), dir.path().join("%GENE.json"), Some(gtr.clone()))?;
    args.reconstructed_fasta_output.output_sequences = Some(nuc_fasta.clone());

    assert_error!(
      run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &NoopProgress),
      format!(
        "Output destination '{}' is selected more than once \
         (--output-model and --output-translations for CDS 'S')",
        gtr.display()
      )
    );
    assert_eq!((false, false), (gtr.exists(), nuc_fasta.exists()));
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::commands::shared::gtr_output::GtrOutputArgs;
    use crate::commands::shared::tree_input::TreeArgs;

    pub(super) fn args(dir: &Path, aa_fasta: PathBuf, gtr: Option<PathBuf>) -> Result<TreetimeAncestralArgs, Report> {
      let tree_path = dir.join("tree.nwk");
      let fasta_path = dir.join("aln.fasta");
      fs::write(&tree_path, TREE).wrap_err("When writing the tree fixture")?;
      fs::write(&fasta_path, ALIGNMENT).wrap_err("When writing the alignment fixture")?;
      fs::write(dir.join("S.fasta"), TRANSLATION).wrap_err("When writing the translation fixture")?;
      TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw {
        alignment: AlignmentArgs {
          alignment: vec![fasta_path],
        },
        tree: TreeArgs { tree: Some(tree_path) },
        model_args: ModelArgs {
          model: Some(GtrModelNameCli::JC69),
          ..ModelArgs::default()
        },
        translations: Some(format!("{}/%GENE.fasta", dir.display())),
        cds_selection: vec_of_owned!["S"],
        output_translations: Some(aa_fasta),
        gtr_output: GtrOutputArgs { output_model: gtr },
        ..TreetimeAncestralArgsRaw::default()
      })
    }
  }
}

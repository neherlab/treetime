#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::ancestral::run::run_ancestral_reconstruction;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::model::{GtrModelNameCli, ModelArgs};
  use clap::Parser;
  use eyre::{Report, WrapErr};
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
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
    let mut raw = helpers::raw(dir.path(), dir.path().join("%GENE.json"), Some(gtr.clone()))?;
    raw.reconstructed_fasta_output.output_sequences = Some(nuc_fasta.clone());
    let args = TreetimeAncestralArgs::try_from(raw)?;

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

  #[test]
  fn test_reconstructed_aa_fasta_of_an_alignment_split_into_partitions_equals_the_unsplit_run() -> Result<(), Report> {
    let dir = tempdir().wrap_err("When creating a temporary directory")?;
    helpers::write_fixtures(dir.path())?;
    let whole = dir.path().join("whole");
    let split = dir.path().join("split");

    helpers::run_with_output_all(dir.path(), &whole, &["--model", "jc69"])?;
    helpers::run_with_output_all(
      dir.path(),
      &split,
      &[
        "--partition",
        "{name: cp1, sites: ['1-6\\3'], model: jc69}",
        "--partition",
        "{name: cp2, sites: ['2-6\\3'], model: jc69}",
        "--partition",
        "{name: cp3, sites: ['3-6\\3'], model: jc69}",
      ],
    )?;

    let whole_files = helpers::output_files(&whole)?;
    assert_eq!(
      (true, whole_files.keys().collect::<Vec<_>>()),
      (
        whole_files.contains_key("translations.S.fasta"),
        helpers::output_files(&split)?.keys().collect::<Vec<_>>()
      )
    );
    let split_files = helpers::output_files(&split)?;
    for (name, content) in whole_files.iter().filter(|(name, _)| !name.ends_with(".model.json")) {
      assert_eq!(content, &split_files[name], "output file {name}");
    }
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::commands::shared::gtr_output::GtrOutputArgs;
    use crate::commands::shared::tree_input::TreeArgs;

    pub(super) fn write_fixtures(dir: &Path) -> Result<(), Report> {
      fs::write(dir.join("tree.nwk"), TREE).wrap_err("When writing the tree fixture")?;
      fs::write(dir.join("aln.fasta"), ALIGNMENT).wrap_err("When writing the alignment fixture")?;
      fs::write(dir.join("S.fasta"), TRANSLATION).wrap_err("When writing the translation fixture")
    }

    pub(super) fn run_with_output_all(dir: &Path, out: &Path, extra: &[&str]) -> Result<(), Report> {
      let path = |name: &str| dir.join(name).display().to_string();
      let output_all = out.display().to_string();
      let translations = format!("{}/%GENE.fasta", dir.display());
      let output_translations = out.join("translations.%GENE.fasta").display().to_string();
      let raw = TreetimeAncestralArgsRaw::try_parse_from(
        [
          &[
            "ancestral",
            "--tree",
            &path("tree.nwk"),
            "--alignment",
            &path("aln.fasta"),
            "--translations",
            &translations,
            "--cds-selection",
            "S",
            "--output-all",
            &output_all,
            "--output-translations",
            &output_translations,
          ][..],
          extra,
        ]
        .concat(),
      )?;
      run_ancestral_reconstruction(
        &TreetimeAncestralArgs::try_from(raw)?,
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )
    }

    pub(super) fn output_files(out: &Path) -> Result<BTreeMap<String, String>, Report> {
      fs::read_dir(out)
        .wrap_err_with(|| format!("When listing '{}'", out.display()))?
        .map(|entry| {
          let path = entry?.path();
          let name = path
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_default();
          let content = fs::read_to_string(&path).wrap_err_with(|| format!("When reading '{}'", path.display()))?;
          Ok((name, content))
        })
        .collect()
    }

    pub(super) fn args(dir: &Path, aa_fasta: PathBuf, gtr: Option<PathBuf>) -> Result<TreetimeAncestralArgs, Report> {
      TreetimeAncestralArgs::try_from(raw(dir, aa_fasta, gtr)?)
    }

    pub(super) fn raw(dir: &Path, aa_fasta: PathBuf, gtr: Option<PathBuf>) -> Result<TreetimeAncestralArgsRaw, Report> {
      let tree_path = dir.join("tree.nwk");
      let fasta_path = dir.join("aln.fasta");
      fs::write(&tree_path, TREE).wrap_err("When writing the tree fixture")?;
      fs::write(&fasta_path, ALIGNMENT).wrap_err("When writing the alignment fixture")?;
      fs::write(dir.join("S.fasta"), TRANSLATION).wrap_err("When writing the translation fixture")?;
      Ok(TreetimeAncestralArgsRaw {
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

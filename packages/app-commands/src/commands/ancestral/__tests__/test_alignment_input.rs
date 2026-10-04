#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::ancestral::run::run_ancestral_reconstruction;
  use eyre::Report;
  use std::fs;
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_utils::assert_error;

  #[test]
  fn test_alignment_input_is_required() -> Result<(), Report> {
    let dir = tempdir()?;
    fs::write(dir.path().join("tree.nwk"), "(A:0.1,B:0.2)root;\n")?;
    let args = TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw {
      tree: Some(dir.path().join("tree.nwk")),
      output_reconstructed_nuc_fasta: Some(dir.path().join("out.fasta")),
      ..TreetimeAncestralArgsRaw::default()
    })?;

    assert_error!(
      run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &NoopProgress),
      "--alignment is required: pass one or more FASTA files, or '-' to read standard input"
    );
    Ok(())
  }
}

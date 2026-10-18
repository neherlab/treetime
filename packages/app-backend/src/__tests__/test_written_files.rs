#[cfg(test)]
mod tests {
  use crate::command::{AppCommandExt, OutputFile, written_files};
  use app_commands::command::AppCommand;
  use app_output::output_plan::OutputSelection;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::fs;
  use tempfile::TempDir;

  #[test]
  fn test_written_files_lists_every_partition_file() -> Result<(), Report> {
    let dir = TempDir::new()?;
    let config = json!({
      "tree": "t.nwk",
      "alignment": ["a.fasta", "b.fasta"],
      "output_selection": ["reconstructed-nuc-fasta", "gtr"],
    });
    let prepared = AppCommand::Ancestral.prepare_run(&config, dir.path())?;
    let names = [
      "ancestral.1.auspice.json",
      "ancestral.1.reconstructed-nuc.fasta",
      "ancestral.2.auspice.json",
      "ancestral.2.reconstructed-nuc.fasta",
      "ancestral.gtr.json",
    ];
    for name in names {
      fs::write(dir.path().join(name), "")?;
    }

    let actual = written_files(&prepared.args.resolve_outputs()?)?;

    let kinds = [
      OutputSelection::Auspice,
      OutputSelection::ReconstructedNucFasta,
      OutputSelection::Auspice,
      OutputSelection::ReconstructedNucFasta,
      OutputSelection::Gtr,
    ];
    let expected = names
      .into_iter()
      .zip(kinds)
      .map(|(name, kind)| OutputFile {
        path: dir.path().join(name),
        kind,
      })
      .collect::<Vec<_>>();
    assert_eq!(expected, actual);
    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use crate::command::{AppCommandExt, OutputAlignment, OutputFile, written_files};
  use app_commands::command::AppCommand;
  use app_output::output_plan::OutputSelection;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::fs;
  use std::path::PathBuf;
  use tempfile::TempDir;

  #[test]
  fn test_written_files_lists_every_alignment_file() -> Result<(), Report> {
    let dir = TempDir::new()?;
    let config = json!({
      "tree": "t.nwk",
      "alignment": ["a.fasta", "b.fasta"],
      "output_selection": ["sequences", "model"],
    });
    let prepared = AppCommand::Ancestral.prepare_run(&config, dir.path())?;
    let names = [
      "ancestral.1.auspice.json",
      "ancestral.1.sequences.fasta",
      "ancestral.2.auspice.json",
      "ancestral.2.sequences.fasta",
      "ancestral.model.json",
    ];
    for name in names {
      fs::write(dir.path().join(name), "")?;
    }

    let actual = written_files(
      &prepared.args.resolve_outputs()?,
      &[PathBuf::from("a.fasta"), PathBuf::from("b.fasta")],
    )?;

    let alignment = |number: usize, path: &str| {
      Some(OutputAlignment {
        number,
        path: PathBuf::from(path),
      })
    };
    let kinds = [
      (OutputSelection::Auspice, alignment(1, "a.fasta")),
      (OutputSelection::Sequences, alignment(1, "a.fasta")),
      (OutputSelection::Auspice, alignment(2, "b.fasta")),
      (OutputSelection::Sequences, alignment(2, "b.fasta")),
      (OutputSelection::Model, None),
    ];
    let expected = names
      .into_iter()
      .zip(kinds)
      .map(|(name, (kind, alignment))| OutputFile {
        path: dir.path().join(name),
        kind,
        alignment,
      })
      .collect::<Vec<_>>();
    assert_eq!(expected, actual);
    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{NUC_ALPHABET, project_root, write_column_blocks};
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::ancestral::run::run_ancestral_reconstruction;
  use clap::Parser;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::fs;
  use std::path::{Path, PathBuf};
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_io::fasta::fasta_read_file;
  use treetime_utils::io::json::json_value_read_file;

  const ZIKA_COLUMNS: usize = 10807;
  const CUT: usize = 5000;

  #[test]
  fn test_ancestral_two_alignments_split_the_columns_of_one_alignment() -> Result<(), Report> {
    let dir = tempdir()?;
    let alignment = project_root().join("data/zika/20/aln.fasta.xz");
    let blocks = write_column_blocks(&alignment, CUT, dir.path())?;
    let whole = dir.path().join("whole.fasta");

    helpers::run_ancestral(
      &[alignment],
      &[format!("--output-reconstructed-nuc-fasta={}", whole.display())],
    )?;
    helpers::run_ancestral(
      &blocks,
      &[format!(
        "--output-reconstructed-nuc-fasta={}",
        dir.path().join("block.{partition}.fasta").display()
      )],
    )?;

    let whole = helpers::read_sequences(&whole)?;
    let first = helpers::read_sequences(&dir.path().join("block.1.fasta"))?;
    let second = helpers::read_sequences(&dir.path().join("block.2.fasta"))?;
    let joined: BTreeMap<String, String> = first
      .iter()
      .map(|(name, columns)| (name.clone(), format!("{columns}{}", second[name])))
      .collect();
    let lengths = [&first, &second].map(|sequences| sequences.values().map(String::len).unique().collect_vec());

    let expected = (
      [vec![CUT], vec![ZIKA_COLUMNS - CUT]],
      whole.keys().collect_vec(),
      whole.clone(),
    );
    assert_eq!(expected, (lengths, second.keys().collect_vec(), joined));
    Ok(())
  }

  #[test]
  fn test_ancestral_two_alignments_write_each_sequence_output_once_per_alignment() -> Result<(), Report> {
    let dir = tempdir()?;
    let blocks = write_column_blocks(&project_root().join("data/zika/20/aln.fasta.xz"), CUT, dir.path())?;
    let out = dir.path().join("out");

    helpers::run_ancestral(
      &blocks,
      &[
        format!("--output-all={}", out.display()),
        "--output-selection=all".to_owned(),
        "--output-nwk-style=plain".to_owned(),
        "--output-nwk-style=beast".to_owned(),
      ],
    )?;

    let files = fs::read_dir(&out)?
      .map(|entry| Ok(entry?.file_name().to_string_lossy().into_owned()))
      .collect::<Result<Vec<_>, Report>>()?
      .into_iter()
      .sorted()
      .collect_vec();
    let lengths = [1, 2]
      .map(|partition| helpers::sequence_output_lengths(&out, partition))
      .into_iter()
      .collect::<Result<Vec<_>, Report>>()?;
    let expected = (
      vec![
        "ancestral.1.annotated.nexus",
        "ancestral.1.annotated.nwk",
        "ancestral.1.augur-node-data.json",
        "ancestral.1.auspice.json",
        "ancestral.1.mat.json",
        "ancestral.1.mat.pb",
        "ancestral.1.reconstructed-nuc.fasta",
        "ancestral.2.annotated.nexus",
        "ancestral.2.annotated.nwk",
        "ancestral.2.augur-node-data.json",
        "ancestral.2.auspice.json",
        "ancestral.2.mat.json",
        "ancestral.2.mat.pb",
        "ancestral.2.reconstructed-nuc.fasta",
        "ancestral.dot",
        "ancestral.graph.json",
        "ancestral.gtr.json",
        "ancestral.nexus",
        "ancestral.nwk",
      ],
      vec![[CUT; 3], [ZIKA_COLUMNS - CUT; 3]],
    );
    assert_eq!(expected, (files.iter().map(String::as_str).collect_vec(), lengths));
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn run_ancestral(alignments: &[PathBuf], outputs: &[String]) -> Result<(), Report> {
      let tree = project_root().join("data/zika/20/tree.nwk");
      let args = ["ancestral".to_owned(), format!("--tree={}", tree.display())]
        .into_iter()
        .chain(alignments.iter().map(|path| format!("--alignment={}", path.display())))
        .chain(["--model=jc69".to_owned(), "--include-leaves".to_owned()])
        .chain(outputs.iter().cloned());
      let args = TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw::try_parse_from(args)?)?;
      run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &NoopProgress)
    }

    pub(super) fn read_sequences(path: &Path) -> Result<BTreeMap<String, String>, Report> {
      Ok(
        fasta_read_file(path, &*NUC_ALPHABET)?
          .into_iter()
          .map(|record| (record.seq_name, record.seq.to_string()))
          .collect(),
      )
    }

    pub(super) fn sequence_output_lengths(out: &Path, partition: usize) -> Result<[usize; 3], Report> {
      let node_data = json_value_read_file(out.join(format!("ancestral.{partition}.augur-node-data.json")))?;
      let auspice = json_value_read_file(out.join(format!("ancestral.{partition}.auspice.json")))?;
      let fasta = read_sequences(&out.join(format!("ancestral.{partition}.reconstructed-nuc.fasta")))?;
      let annotated_end = node_data["annotations"]["nuc"]["end"].as_u64().unwrap_or_default();
      let auspice_root = auspice["root_sequence"]["nuc"].as_str().map_or(0, str::len);
      Ok([
        usize::try_from(annotated_end)?,
        auspice_root,
        fasta.values().map(String::len).max().unwrap_or_default(),
      ])
    }
  }
}

#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::RecordingLog;
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::ancestral::run::run_ancestral_reconstruction;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::method_anc::{MethodAncArgs, MethodAncestralCli};
  use crate::commands::shared::reconstructed_fasta::ReconstructedFastaOutputArgs;
  use crate::commands::shared::tree_input::{SegmentArgs, TreeArgs};
  use eyre::{Report, WrapErr};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::fs;
  use std::path::{Path, PathBuf};
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_utils::assert_error;

  const SEGMENTS: [&str; 3] = ["seg1", "seg2", "seg3"];

  fn toy(file: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .join("../../data/toy/reassortment")
      .join(file)
  }

  fn segment_files(prefix: &str, extension: &str) -> Vec<PathBuf> {
    SEGMENTS
      .iter()
      .map(|segment| toy(&format!("{prefix}_{segment}.{extension}")))
      .collect()
  }

  #[expect(
    clippy::literal_string_with_formatting_args,
    reason = "the braces are a path template placeholder, not a format argument"
  )]
  fn raw(
    tree: Option<&str>,
    alignment: &[&str],
    segment_trees: &[&str],
    segment_alignments: &[&str],
  ) -> TreetimeAncestralArgsRaw {
    let paths = |paths: &[&str]| paths.iter().map(PathBuf::from).collect::<Vec<_>>();
    TreetimeAncestralArgsRaw {
      tree: TreeArgs {
        tree: tree.map(PathBuf::from),
      },
      alignment: AlignmentArgs {
        alignment: paths(alignment),
      },
      segments: SegmentArgs {
        segment_trees: paths(segment_trees),
        segment_alignments: paths(segment_alignments),
      },
      method_anc: MethodAncArgs {
        method_anc: MethodAncestralCli::Parsimony,
      },
      reconstructed_fasta_output: ReconstructedFastaOutputArgs {
        output_sequences: Some(PathBuf::from("out.{alignment}.fasta")),
      },
      ..TreetimeAncestralArgsRaw::default()
    }
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::count_mismatch(    None,           &[],        &["a.nwk", "b.nwk"], &["a.fasta"],            "--segment-trees and --segment-alignments pair one tree with one alignment, but they list 2 trees and 1 alignments")]
  #[case::with_tree(         Some("t.nwk"),  &[],        &["a.nwk"],          &["a.fasta"],            "--segment-trees replaces --tree: give one of them")]
  #[case::with_alignment(    None,           &["x.fa"],  &["a.nwk"],          &["a.fasta"],            "--segment-trees takes its alignments from --segment-alignments, which replaces --alignment: give one of them")]
  #[case::alignments_alone(  Some("t.nwk"),  &[],        &[],                 &["a.fasta"],            "--segment-alignments requires --segment-trees, one tree per alignment")]
  #[trace]
  fn test_segment_trees_arguments_are_checked(
    #[case] tree: Option<&str>,
    #[case] alignment: &[&str],
    #[case] segment_trees: &[&str],
    #[case] segment_alignments: &[&str],
    #[case] expected: &str,
  ) {
    assert_error!(
      TreetimeAncestralArgs::try_from(raw(tree, alignment, segment_trees, segment_alignments)),
      expected
    );
  }

  #[test]
  fn test_segment_trees_with_translations_is_error() {
    let raw = TreetimeAncestralArgsRaw {
      translations: Some("aa_{cds}.fasta".to_owned()),
      ..raw(None, &[], &["a.nwk"], &["a.fasta"])
    };
    assert_error!(
      TreetimeAncestralArgs::try_from(raw),
      "--translations needs one tree, so it does not combine with --segment-trees yet"
    );
  }

  #[test]
  fn test_segment_trees_parsimony_writes_the_sequences_of_each_separate_segment_run() -> Result<(), Report> {
    let dir = tempdir().wrap_err("When creating a temporary directory")?;
    let template = dir.path().join("merged.{alignment}.fasta");
    let args = TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw {
      segments: SegmentArgs {
        segment_trees: segment_files("tree", "nwk"),
        segment_alignments: segment_files("aln", "fasta"),
      },
      reconstructed_fasta_output: ReconstructedFastaOutputArgs {
        output_sequences: Some(template),
      },
      ..raw(None, &[], &[], &[])
    })?;
    run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &RecordingLog::default())?;

    for (index, segment) in SEGMENTS.iter().enumerate() {
      let separate = run_separate(dir.path(), segment)?;
      let merged = records(&fs::read_to_string(
        dir.path().join(format!("merged.{}.fasta", index + 1)),
      )?);
      assert_eq!(
        (false, records(&separate)),
        (merged.is_empty(), merged),
        "segment {segment}"
      );
    }
    Ok(())
  }

  fn records(fasta: &str) -> BTreeMap<String, String> {
    fasta
      .split('>')
      .filter_map(|record| record.split_once('\n'))
      .map(|(name, seq)| (name.trim().to_owned(), seq.replace('\n', "")))
      .collect()
  }

  fn run_separate(dir: &Path, segment: &str) -> Result<String, Report> {
    let path = dir.join(format!("separate_{segment}.fasta"));
    let args = TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw {
      tree: TreeArgs {
        tree: Some(toy(&format!("tree_{segment}.nwk"))),
      },
      alignment: AlignmentArgs {
        alignment: vec![toy(&format!("aln_{segment}.fasta"))],
      },
      reconstructed_fasta_output: ReconstructedFastaOutputArgs {
        output_sequences: Some(path.clone()),
      },
      ..raw(None, &[], &[], &[])
    })?;
    run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &RecordingLog::default())?;
    fs::read_to_string(&path).wrap_err("When reading the separate segment run")
  }
}

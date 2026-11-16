#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::RecordingLog;
  use crate::commands::shared::metadata::MetadataArgs;
  use crate::commands::shared::output_args::OutputCoreArgs;
  use crate::commands::shared::tree_input::{SegmentArgs, TreeArgs};
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use crate::commands::timetree::coalescent::CoalescentArg;
  use crate::commands::timetree::run::run_timetree_estimation;
  use eyre::{Report, WrapErr};
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::fs;
  use std::path::{Path, PathBuf};
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;

  const SEGMENTS: [&str; 3] = ["seg1", "seg2", "seg3"];

  const TOY_CLOCK_RATE: f64 = 0.005;

  fn toy(file: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .join("../../data/toy/reassortment")
      .join(file)
  }

  fn segment_raw(dates: &Path, output: OutputCoreArgs) -> TreetimeTimetreeArgsRaw {
    TreetimeTimetreeArgsRaw {
      segments: SegmentArgs {
        segment_trees: SEGMENTS
          .iter()
          .map(|segment| toy(&format!("tree_{segment}.nwk")))
          .collect(),
        segment_alignments: SEGMENTS
          .iter()
          .map(|segment| toy(&format!("aln_{segment}.fasta")))
          .collect(),
      },
      metadata: MetadataArgs {
        metadata: Some(toy("metadata.tsv")),
      },
      clock_rate: Some(TOY_CLOCK_RATE),
      confidence: true,
      output_dates: Some(dates.to_path_buf()),
      output,
      ..TreetimeTimetreeArgsRaw::default()
    }
  }

  fn read_table(path: &Path) -> Result<BTreeMap<String, Vec<String>>, Report> {
    Ok(
      fs::read_to_string(path)?
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
          let mut fields = line.split('\t').map(str::to_owned);
          let name = fields.next()?;
          Some((name, fields.collect()))
        })
        .collect(),
    )
  }

  #[test]
  fn test_segment_trees_timetree_dates_every_node_and_covers_the_true_dates() -> Result<(), Report> {
    let dir = tempdir().wrap_err("When creating a temporary directory")?;
    let dates = dir.path().join("dates.tsv");
    let args = TreetimeTimetreeArgs::try_from(segment_raw(
      &dates,
      OutputCoreArgs {
        output_tree_nexus: Some(dir.path().join("timetree.{alignment}.nexus")),
        ..OutputCoreArgs::default()
      },
    ))?;
    run_timetree_estimation(&args, &NoopCancel, &NoopProgress, &RecordingLog::default())?;

    let inferred = read_table(&dates)?;
    let mut truth = read_table(&toy("nodes.tsv"))?;
    truth.remove("name");
    let internal = truth
      .iter()
      .filter(|(name, _)| name.len() > 1 || name.as_str() >= "K")
      .map(|(name, fields)| {
        let date: f64 = fields[1].parse()?;
        let row = &inferred[name];
        let lower: f64 = row[2].parse()?;
        let upper: f64 = row[3].parse()?;
        Ok((lower..=upper).contains(&date))
      })
      .collect::<Result<Vec<bool>, Report>>()?;
    let covered = internal.iter().filter(|&&inside| inside).count();
    let nexus = (1..=SEGMENTS.len())
      .filter(|number| dir.path().join(format!("timetree.{number}.nexus")).exists())
      .count();
    assert_eq!(
      (truth.len(), 12, true, SEGMENTS.len()),
      (inferred.len(), internal.len(), covered >= 10, nexus)
    );
    Ok(())
  }

  #[test]
  fn test_segment_trees_timetree_requires_a_clock_rate() {
    let raw = TreetimeTimetreeArgsRaw {
      clock_rate: None,
      ..segment_raw(Path::new("dates.tsv"), OutputCoreArgs::default())
    };
    let error = TreetimeTimetreeArgs::try_from(raw).err().map(|report| report.to_string());
    assert!(
      error.as_deref().is_some_and(|message| message.contains("--clock-rate")),
      "{error:?}"
    );
  }

  #[test]
  fn test_segment_trees_timetree_rejects_a_coalescent_prior() {
    let raw = TreetimeTimetreeArgsRaw {
      coalescent: Some(CoalescentArg::Skyline),
      ..segment_raw(Path::new("dates.tsv"), OutputCoreArgs::default())
    };
    let error = TreetimeTimetreeArgs::try_from(raw).err().map(|report| report.to_string());
    assert!(
      error.as_deref().is_some_and(|message| message.contains("--coalescent")),
      "{error:?}"
    );
  }

  #[test]
  fn test_tree_timetree_rejects_the_cross_parent_switch() {
    let raw = TreetimeTimetreeArgsRaw {
      cross_parent_messages: true,
      segments: SegmentArgs::default(),
      tree: TreeArgs {
        tree: Some(toy("tree_seg1.nwk")),
      },
      ..segment_raw(Path::new("dates.tsv"), OutputCoreArgs::default())
    };
    let error = TreetimeTimetreeArgs::try_from(raw).err().map(|report| report.to_string());
    assert!(
      error.as_deref().is_some_and(|message| message.contains("segment")),
      "{error:?}"
    );
  }
}

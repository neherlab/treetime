#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::RecordingLog;
  use crate::commands::optimize::args::{TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw};
  use crate::commands::optimize::run::run_optimize;
  use crate::commands::shared::output_args::OutputCoreArgs;
  use crate::commands::shared::tree_input::SegmentArgs;
  use eyre::{Report, WrapErr};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::collections::{BTreeMap, BTreeSet};
  use std::fs;
  use std::path::PathBuf;
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_io::nwk::TREE_DIALECT_DEFAULT;
  use treetime_io::tree::tree_read;
  use treetime_utils::make_report;

  const SEGMENTS: [&str; 3] = ["seg1", "seg2", "seg3"];

  fn toy(file: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .join("../../data/toy/reassortment")
      .join(file)
  }

  fn edge_lengths(newick: &str) -> Result<BTreeMap<(String, String), f64>, Report> {
    let tree = tree_read(newick.as_bytes(), TREE_DIALECT_DEFAULT)?;
    tree
      .graph
      .get_edges()
      .map(|edge| {
        let length = tree.branch_lengths[&edge.key()].ok_or_else(|| make_report!("edge without length"))?;
        Ok((
          (
            tree.names[edge.source()].to_owned(),
            tree.names[edge.target()].to_owned(),
          ),
          length,
        ))
      })
      .collect()
  }

  #[test]
  fn test_segment_trees_optimize_writes_one_tree_per_segment_with_one_length_per_shared_edge() -> Result<(), Report> {
    let dir = tempdir().wrap_err("When creating a temporary directory")?;
    let args = TreetimeOptimizeArgs::try_from(TreetimeOptimizeArgsRaw {
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
      output: OutputCoreArgs {
        output_tree_nwk: Some(dir.path().join("optimized.{alignment}.nwk")),
        ..OutputCoreArgs::default()
      },
      ..TreetimeOptimizeArgsRaw::default()
    })?;
    run_optimize(&args, &NoopCancel, &NoopProgress, &RecordingLog::default())?;

    let trees = (1..=SEGMENTS.len())
      .map(|number| edge_lengths(&fs::read_to_string(dir.path().join(format!("optimized.{number}.nwk")))?))
      .collect::<Result<Vec<_>, Report>>()?;
    let mut by_edge: BTreeMap<&(String, String), Vec<f64>> = BTreeMap::new();
    for tree in &trees {
      for (edge, &length) in tree {
        by_edge.entry(edge).or_default().push(length);
      }
    }
    let shared_with_one_length = by_edge
      .iter()
      .filter(|(_, lengths)| lengths.len() > 1)
      .map(|(edge, lengths)| {
        (
          (*edge).clone(),
          lengths.iter().all(|length| length.to_bits() == lengths[0].to_bits()),
        )
      })
      .collect_vec();
    let leaves = trees
      .iter()
      .map(|tree| {
        tree
          .keys()
          .map(|(_, child)| child.as_str())
          .filter(|name| name.len() == 1 && *name < "K")
          .collect::<BTreeSet<_>>()
      })
      .collect_vec();
    assert_eq!(
      (
        18,
        true,
        vec![BTreeSet::from(["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]); 3]
      ),
      (
        shared_with_one_length.len(),
        shared_with_one_length.iter().all(|(_, same)| *same),
        leaves
      )
    );
    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::ancestral::run::run_ancestral_reconstruction;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::model::{GtrModelNameCli, ModelArgs};
  use crate::commands::shared::output_args::OutputCoreArgs;
  use crate::job::{JobEvent, JobProgress};
  use eyre::{Report, WrapErr};
  use maplit::btreemap;
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::fs;
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::{LogLevel, NoopProgress};
  use treetime_io::usher_mat::UsherTree;
  use treetime_utils::io::json::json_read_file;
  use treetime_utils::{o, vec_of_owned};

  const TREE: &str = "((A:0.1,B:0.1)AB:0.1,(C:0.1,D:0.1)CD:0.1)root;";
  const ALIGNMENT: &str = ">A\nA-GTTA\n>B\nACG-TC\n>C\nACG-TC\n>D\nACG-TC\n";

  #[rustfmt::skip]
  #[rstest]
  #[case::dense( Some(true))]
  #[case::sparse(Some(false))]
  #[trace]
  fn test_mat_gaps_written_as_missing_data_with_one_warning(#[case] dense: Option<bool>) -> Result<(), Report> {
    let dir = tempdir().wrap_err("When creating a temporary directory")?;
    let tree_path = dir.path().join("tree.nwk");
    let fasta_path = dir.path().join("aln.fasta");
    let mat_json_path = dir.path().join("out.mat.json");
    let mat_pb_path = dir.path().join("out.mat.pb");
    fs::write(&tree_path, TREE).wrap_err("When writing the tree fixture")?;
    fs::write(&fasta_path, ALIGNMENT).wrap_err("When writing the alignment fixture")?;
    let args = TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw {
      alignment: AlignmentArgs {
        alignment: vec![fasta_path],
      },
      tree: Some(tree_path),
      model_args: ModelArgs {
        model: GtrModelNameCli::JC69,
        ..ModelArgs::default()
      },
      dense,
      output: OutputCoreArgs {
        output_tree_mat_json: Some(mat_json_path.clone()),
        output_tree_mat_pb: Some(mat_pb_path.clone()),
        ..OutputCoreArgs::default()
      },
      ..TreetimeAncestralArgsRaw::default()
    })?;
    let warnings = Mutex::new(vec![]);
    let log = JobProgress::new(|event: JobEvent| {
      if let JobEvent::Log { data: event } = event
        && event.level == LogLevel::Warn
      {
        warnings.lock().push(event.message);
      }
    });

    run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &log)?;

    let mat: UsherTree = json_read_file(&mat_json_path)?;
    let actual: BTreeMap<String, Vec<(i32, i32, i32, Vec<i32>)>> = mat
      .condensed_nodes
      .into_iter()
      .zip(mat.node_mutations)
      .filter(|(_, mutations)| !mutations.mutation.is_empty())
      .map(|(node, mutations)| {
        let mutations = mutations
          .mutation
          .into_iter()
          .map(|mutation| (mutation.position, mutation.ref_nuc, mutation.par_nuc, mutation.mut_nuc))
          .collect();
        (node.node_name, mutations)
      })
      .collect();
    let expected = btreemap! { o!("A") => vec![(6, 1, 1, vec![0])] };
    assert_eq!(expected, actual);
    let expected_warnings = vec_of_owned![
      "UShER MAT has no gap state: wrote 1 deletion(s) as missing data (N); left out 1 insertion(s) and 0 substitution(s) in alignment columns where the root sequence, which is the MAT reference, has a gap"
    ];
    assert_eq!(expected_warnings, warnings.into_inner());
    assert!(mat_pb_path.is_file());
    Ok(())
  }
}

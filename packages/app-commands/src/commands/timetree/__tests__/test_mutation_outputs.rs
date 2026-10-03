#[cfg(test)]
mod tests {
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::branch_length_mode::BranchLengthModeCli;
  use crate::commands::shared::output_args::{DivergenceUnits, OutputCoreArgs};
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use crate::commands::timetree::run::run_timetree_estimation;
  use eyre::{Report, WrapErr};
  use helpers::{node_mutations, run_timetree};
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_io::auspice_types::{AuspiceTree, AuspiceTreeNode};
  use treetime_utils::io::json::json_read_file;
  use treetime_utils::{assert_error, o, vec_of_owned};
  use util_augur_node_data_json::AugurNodeDataJsonRefine;

  const TREE: &str = "((A:0.06,B:0.09)AB:0.03,(C:0.07,(D:0.03,E:0.06)DE:0.02)CDE:0.03,F:0.04)root;";

  const ALIGNMENT: &str = indoc! {"
    >A
    ACGTACGTGTATACGTACGTGTATACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
    >B
    ACGTACGTGTATACGTACGTACGTACGTACACGTACACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
    >C
    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGTATACGTACACGTGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
    >D
    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGTATACGTACGTACGTACGTGTGTACGTACACGCGTAKGNACGTACGTACGTACGTACGT
    >E
    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGTATACGTACGTACGTACGTGTGTACGTACGTACGTACGTACGTGTACGTGTACGTACGT
    >F
    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACACGTGT
  "};

  const METADATA: &str = "strain\tdate\nA\t2006\nB\t2009\nC\t2007\nD\t2008\nE\t2011\nF\t2004\n";

  #[test]
  fn test_timetree_mutation_units_count_only_state_changes() -> Result<(), Report> {
    let dir = tempfile::tempdir().wrap_err("When creating a temporary directory")?;
    let augur_path = dir.path().join("node_data.json");
    run_timetree(dir.path(), false, |raw| {
      raw.divergence_units = DivergenceUnits::Mutations;
      raw.output_augur_node_data = Some(augur_path.clone());
    })?;

    let node_data: AugurNodeDataJsonRefine = json_read_file(&augur_path)?;
    let actual: BTreeMap<String, Option<f64>> = node_data
      .nodes
      .into_iter()
      .map(|(name, node)| (name, node.mutation_length))
      .collect();
    let expected = btreemap! {
      o!("root") => Some(0.0),
      o!("AB") => Some(3.0),
      o!("A") => Some(3.0),
      o!("B") => Some(6.0),
      o!("CDE") => Some(3.0),
      o!("C") => Some(4.0),
      o!("DE") => Some(2.0),
      o!("D") => Some(4.0),
      o!("E") => Some(6.0),
      o!("F") => Some(4.0),
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_timetree_mutation_units_reject_input_branch_lengths_before_inference() {
    let dir = tempfile::tempdir().unwrap();
    let result = run_timetree(dir.path(), false, |raw| {
      raw.divergence_units = DivergenceUnits::Mutations;
      raw.branch_length_mode = BranchLengthModeCli::Input;
      raw.output_augur_node_data = Some(dir.path().join("node_data.json"));
    });
    assert_error!(
      result,
      "--divergence-units=mutations requires ancestral reconstruction; incompatible with --branch-length-mode=input"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::default_hides_unknown( false, vec_of_owned!["G71A", "T72C", "A73G", "C78K"])]
  #[case::report_ambiguous(      true,  vec_of_owned!["G71A", "T72C", "A73G", "C78K", "T80N"])]
  #[trace]
  fn test_timetree_report_ambiguous_reaches_auspice_mutations(
    #[case] report_ambiguous: bool,
    #[case] expected_d: Vec<String>,
  ) -> Result<(), Report> {
    let dir = tempfile::tempdir().wrap_err("When creating a temporary directory")?;
    let auspice_path = dir.path().join("timetree.auspice.json");
    run_timetree(dir.path(), report_ambiguous, |raw| {
      raw.output = OutputCoreArgs {
        output_tree_auspice: Some(auspice_path.clone()),
        ..OutputCoreArgs::default()
      };
    })?;

    let auspice: AuspiceTree = json_read_file(&auspice_path)?;
    let expected = btreemap! {
      o!("AB") => vec_of_owned!["A9G", "C10T", "G11A"],
      o!("A") => vec_of_owned!["A21G", "C22T", "G23A"],
      o!("B") => vec_of_owned!["G31A", "T32C", "A33G", "C34T", "G35A", "T36C"],
      o!("CDE") => vec_of_owned!["A41G", "C42T", "G43A"],
      o!("C") => vec_of_owned!["G51A", "T52C", "A53G", "C54T"],
      o!("DE") => vec_of_owned!["A61G", "C62T"],
      o!("D") => expected_d,
      o!("E") => vec_of_owned!["A85G", "C86T", "G87A", "T88C", "A89G", "C90T"],
      o!("F") => vec_of_owned!["G95A", "T96C", "A97G", "C98T"],
    };
    assert_eq!(expected, node_mutations(&auspice.tree));
    Ok(())
  }

  mod helpers {
    use super::*;
    use std::path::Path;

    pub(super) fn run_timetree(
      dir: &Path,
      report_ambiguous: bool,
      configure: impl FnOnce(&mut TreetimeTimetreeArgsRaw),
    ) -> Result<(), Report> {
      let tree_path = dir.join("tree.nwk");
      let fasta_path = dir.join("aln.fasta");
      let metadata_path = dir.join("metadata.tsv");
      std::fs::write(&tree_path, TREE).wrap_err("When writing the tree fixture")?;
      std::fs::write(&fasta_path, ALIGNMENT).wrap_err("When writing the alignment fixture")?;
      std::fs::write(&metadata_path, METADATA).wrap_err("When writing the metadata fixture")?;

      let mut raw = TreetimeTimetreeArgsRaw {
        alignment: AlignmentArgs {
          alignment: vec![fasta_path],
        },
        tree: Some(tree_path),
        metadata: Some(metadata_path),
        keep_root: true,
        clock_rate: Some(0.01),
        max_iter: 2,
        report_ambiguous,
        ..TreetimeTimetreeArgsRaw::default()
      };
      configure(&mut raw);
      let args = TreetimeTimetreeArgs::try_from(raw)?;
      run_timetree_estimation(&args, &NoopCancel, &NoopProgress, &NoopProgress)
    }

    pub(super) fn node_mutations(node: &AuspiceTreeNode) -> BTreeMap<String, Vec<String>> {
      let own = node
        .branch_attrs
        .mutations
        .get("nuc")
        .filter(|mutations| !mutations.is_empty())
        .map(|mutations| (node.name.clone(), mutations.clone()));
      own
        .into_iter()
        .chain(node.children.iter().flat_map(node_mutations))
        .collect()
    }
  }
}

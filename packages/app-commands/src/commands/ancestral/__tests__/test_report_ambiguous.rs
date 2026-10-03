#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::ancestral::run::run_ancestral_reconstruction;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::model::{GtrModelNameCli, ModelArgs};
  use crate::commands::shared::output_args::{NwkStyleArg, OutputCoreArgs};
  use eyre::{Report, WrapErr};
  use helpers::written_mutations;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_io::usher_mat::UsherTree;
  use treetime_utils::io::fs::read_file_to_string;
  use treetime_utils::io::json::json_read_file;
  use treetime_utils::{o, vec_of_owned};
  use util_augur_node_data_json::AugurNodeDataJsonAncestral;

  const STAR_TREE: &str = "(A:0.1,B:0.1,C:0.1)root;";
  const STAR_ALIGNMENT: &str = ">A\nACGT\n>B\nACGT\n>C\nANGK\n";

  const UNKNOWN_CLADE_TREE: &str = "(((C1:0.02,C2:0.03)P:0.02,C3:0.04)P2:0.01,((D:0.02,E:0.03)Q:0.02,F:0.02)Q2:0.02)R;";
  const UNKNOWN_CLADE_ALIGNMENT: &str = ">C1\nANGT\n>C2\nANGT\n>C3\nANGT\n>D\nAAGT\n>E\nAAGT\n>F\nAAGT\n";
  const UNKNOWN_CLADE_WITH_SUB_ALIGNMENT: &str = ">C1\nANTT\n>C2\nANGT\n>C3\nANGT\n>D\nAAGT\n>E\nAAGT\n>F\nAAGT\n";

  #[rustfmt::skip]
  #[rstest]
  #[case::default_drops_unknown(  false, vec!["T4K"],        "(A:0.1,B:0.1,C[&mutations=T4K]:0.1)root;")]
  #[case::report_ambiguous(       true,  vec!["C2N", "T4K"], r#"(A:0.1,B:0.1,C[&mutations="C2N,T4K"]:0.1)root;"#)]
  #[trace]
  fn test_report_ambiguous_reaches_mutation_writers(
    #[case] report_ambiguous: bool,
    #[case] expected: Vec<&str>,
    #[case] expected_nwk: &str,
  ) -> Result<(), Report> {
    let written = written_mutations(STAR_TREE, STAR_ALIGNMENT, report_ambiguous, false, None)?;
    assert_eq!(expected, written.augur_muts.get("C").cloned().unwrap_or_default());
    assert_eq!(expected_nwk, written.nwk.trim());
    Ok(())
  }

  #[rstest]
  #[case::dense_default_hides_unknown(Some(true), false, false, btreemap! {})]
  #[case::sparse_default_hides_unknown(Some(false), false, false, btreemap! {})]
  #[case::dense_default_hides_imputed_from_unknown(Some(true), false, true, btreemap! {})]
  #[case::dense_reports_unknown_above_clade(Some(true), true, false, btreemap! {
    o!("P2") => vec_of_owned!["A2N"],
  })]
  #[case::sparse_reports_unknown_above_clade(Some(false), true, false, btreemap! {
    o!("P2") => vec_of_owned!["A2N"],
  })]
  #[case::dense_reports_imputed_tips_below_unknown(Some(true), true, true, btreemap! {
    o!("C1") => vec_of_owned!["N2A"],
    o!("C2") => vec_of_owned!["N2A"],
    o!("C3") => vec_of_owned!["N2A"],
    o!("P2") => vec_of_owned!["A2N"],
  })]
  #[trace]
  fn test_report_ambiguous_shows_transitions_of_internal_unknown_state(
    #[case] dense: Option<bool>,
    #[case] report_ambiguous: bool,
    #[case] impute: bool,
    #[case] expected: BTreeMap<String, Vec<String>>,
  ) -> Result<(), Report> {
    let written = written_mutations(
      UNKNOWN_CLADE_TREE,
      UNKNOWN_CLADE_ALIGNMENT,
      report_ambiguous,
      impute,
      dense,
    )?;
    let actual: BTreeMap<String, Vec<String>> = written
      .augur_muts
      .into_iter()
      .filter(|(_, muts)| !muts.is_empty())
      .collect();
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  #[ignore = "sparse imputation skips tips below an unknown node: kb/issues/M-ancestral-sparse-imputation-skips-tips-below-unknown-node.md"]
  fn test_report_ambiguous_sparse_reports_imputed_tips_below_unknown() -> Result<(), Report> {
    let written = written_mutations(UNKNOWN_CLADE_TREE, UNKNOWN_CLADE_ALIGNMENT, true, true, Some(false))?;
    let actual: BTreeMap<String, Vec<String>> = written
      .augur_muts
      .into_iter()
      .filter(|(_, muts)| !muts.is_empty())
      .collect();
    let expected = btreemap! {
      o!("C1") => vec_of_owned!["N2A"],
      o!("C2") => vec_of_owned!["N2A"],
      o!("C3") => vec_of_owned!["N2A"],
      o!("P2") => vec_of_owned!["A2N"],
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::leaf_unknown(          STAR_TREE,          STAR_ALIGNMENT,                    false, btreemap! { o!("C") => vec![(4, 3)] })]
  #[case::internal_unknown(      UNKNOWN_CLADE_TREE, UNKNOWN_CLADE_ALIGNMENT,           false, btreemap! {})]
  #[case::determined_in_clade(   UNKNOWN_CLADE_TREE, UNKNOWN_CLADE_WITH_SUB_ALIGNMENT,  false, btreemap! { o!("C1") => vec![(3, 2)] })]
  #[case::imputed_below_unknown( UNKNOWN_CLADE_TREE, UNKNOWN_CLADE_ALIGNMENT,           true,  btreemap! {})]
  #[trace]
  fn test_report_ambiguous_mat_stores_unknown_state_as_missing_data(
    #[case] tree: &str,
    #[case] alignment: &str,
    #[case] impute: bool,
    #[case] expected: BTreeMap<String, Vec<(i32, i32)>>,
  ) -> Result<(), Report> {
    let written = written_mutations(tree, alignment, true, impute, Some(true))?;
    assert_eq!(expected, written.mat_mutations);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) struct WrittenMutations {
      pub(super) augur_muts: BTreeMap<String, Vec<String>>,
      pub(super) nwk: String,
      pub(super) mat_mutations: BTreeMap<String, Vec<(i32, i32)>>,
    }

    pub(super) fn written_mutations(
      tree: &str,
      alignment: &str,
      report_ambiguous: bool,
      impute_missing_data: bool,
      dense: Option<bool>,
    ) -> Result<WrittenMutations, Report> {
      let dir = tempdir().wrap_err("When creating a temporary directory")?;
      let tree_path = dir.path().join("tree.nwk");
      let fasta_path = dir.path().join("aln.fasta");
      let augur_path = dir.path().join("node_data.json");
      let nwk_path = dir.path().join("out.nwk");
      let mat_path = dir.path().join("out.mat.json");
      std::fs::write(&tree_path, tree).wrap_err("When writing the tree fixture")?;
      std::fs::write(&fasta_path, alignment).wrap_err("When writing the alignment fixture")?;

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
        impute_missing_data,
        report_ambiguous,
        output_augur_node_data: Some(augur_path.clone()),
        output: OutputCoreArgs {
          output_nwk_style: vec![NwkStyleArg::Beast],
          output_tree_nwk: Some(nwk_path.clone()),
          output_tree_mat_json: Some(mat_path.clone()),
          ..OutputCoreArgs::default()
        },
        ..TreetimeAncestralArgsRaw::default()
      })?;

      run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;

      let node_data: AugurNodeDataJsonAncestral = json_read_file(&augur_path)?;
      let augur_muts = node_data
        .nodes
        .into_iter()
        .map(|(name, node)| (name, node.muts))
        .collect();
      let mat: UsherTree = json_read_file(&mat_path)?;
      let mat_mutations = mat
        .condensed_nodes
        .into_iter()
        .zip(mat.node_mutations)
        .filter(|(_, mutations)| !mutations.mutation.is_empty())
        .map(|(node, mutations)| {
          let positions = mutations
            .mutation
            .iter()
            .map(|mutation| (mutation.position, mutation.par_nuc))
            .collect();
          (node.node_name, positions)
        })
        .collect();
      Ok(WrittenMutations {
        augur_muts,
        nwk: read_file_to_string(&nwk_path)?,
        mat_mutations,
      })
    }
  }
}

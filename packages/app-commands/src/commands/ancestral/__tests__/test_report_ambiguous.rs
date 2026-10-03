#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::ancestral::run::run_ancestral_reconstruction;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::model::{GtrModelNameCli, ModelArgs};
  use crate::commands::shared::output_args::{NwkStyleArg, OutputCoreArgs};
  use eyre::{Report, WrapErr};
  use helpers::written_leaf_mutations;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_utils::io::fs::read_file_to_string;
  use treetime_utils::io::json::json_read_file;
  use util_augur_node_data_json::AugurNodeDataJsonAncestral;

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
    let (augur_muts, nwk) = written_leaf_mutations(report_ambiguous)?;
    assert_eq!(expected, augur_muts);
    assert_eq!(expected_nwk, nwk.trim());
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn written_leaf_mutations(report_ambiguous: bool) -> Result<(Vec<String>, String), Report> {
      let dir = tempdir().wrap_err("When creating a temporary directory")?;
      let tree_path = dir.path().join("tree.nwk");
      let fasta_path = dir.path().join("aln.fasta");
      let augur_path = dir.path().join("node_data.json");
      let nwk_path = dir.path().join("out.nwk");
      std::fs::write(&tree_path, "(A:0.1,B:0.1,C:0.1)root;").wrap_err("When writing the tree fixture")?;
      std::fs::write(&fasta_path, ">A\nACGT\n>B\nACGT\n>C\nANGK\n").wrap_err("When writing the alignment fixture")?;

      let args = TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw {
        alignment: AlignmentArgs {
          alignment: vec![fasta_path],
        },
        tree: Some(tree_path),
        model_args: ModelArgs {
          model: GtrModelNameCli::JC69,
          ..ModelArgs::default()
        },
        report_ambiguous,
        output_augur_node_data: Some(augur_path.clone()),
        output: OutputCoreArgs {
          output_nwk_style: vec![NwkStyleArg::Beast],
          output_tree_nwk: Some(nwk_path.clone()),
          ..OutputCoreArgs::default()
        },
        ..TreetimeAncestralArgsRaw::default()
      })?;

      run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;

      let mut node_data: AugurNodeDataJsonAncestral = json_read_file(&augur_path)?;
      let augur_muts = node_data.nodes.remove("C").map(|node| node.muts).unwrap_or_default();
      Ok((augur_muts, read_file_to_string(&nwk_path)?))
    }
  }
}

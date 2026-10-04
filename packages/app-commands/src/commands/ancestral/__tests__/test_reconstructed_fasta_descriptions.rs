#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::ancestral::run::run_ancestral_reconstruction;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::method_anc::MethodAncestralCli;
  use crate::commands::shared::model::GtrModelNameCli;
  use crate::commands::shared::model::ModelArgs;
  use eyre::{Report, WrapErr};
  use helpers::reconstructed_descriptions;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use tempfile::tempdir;
  use treetime::alphabet::alphabet::Alphabet;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_io::fasta::fasta_read_file;
  use treetime_utils::o;

  #[test]
  fn test_reconstructed_fasta_descriptions_parsimony() -> Result<(), Report> {
    let descriptions = reconstructed_descriptions(MethodAncestralCli::Parsimony, None, true)?;
    assert_eq!(Some(&Some("sample description".to_owned())), descriptions.get("A"));
    assert_eq!(Some(&None), descriptions.get("B"));
    assert_eq!(Some(&None), descriptions.get("root"));
    Ok(())
  }

  #[test]
  fn test_reconstructed_fasta_descriptions_marginal() -> Result<(), Report> {
    let descriptions = reconstructed_descriptions(MethodAncestralCli::Marginal, Some(false), true)?;
    assert_eq!(Some(&Some("sample description".to_owned())), descriptions.get("A"));
    assert_eq!(Some(&None), descriptions.get("B"));
    assert_eq!(Some(&None), descriptions.get("root"));
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::parsimony(      MethodAncestralCli::Parsimony, None)]
  #[case::marginal_sparse(MethodAncestralCli::Marginal,  Some(false))]
  #[case::marginal_dense( MethodAncestralCli::Marginal,  Some(true))]
  #[trace]
  fn test_reconstructed_fasta_without_leaves_holds_only_internal_nodes(
    #[case] method: MethodAncestralCli,
    #[case] dense: Option<bool>,
  ) -> Result<(), Report> {
    let descriptions = reconstructed_descriptions(method, dense, false)?;
    let expected = btreemap! { o!("root") => None };
    assert_eq!(expected, descriptions);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn reconstructed_descriptions(
      method: MethodAncestralCli,
      dense: Option<bool>,
      include_leaves: bool,
    ) -> Result<BTreeMap<String, Option<String>>, Report> {
      let dir = tempdir().wrap_err("When creating a temporary directory")?;
      let tree_path = dir.path().join("tree.nwk");
      let fasta_path = dir.path().join("aln.fasta");
      let out_path = dir.path().join("reconstructed-nuc.fasta");
      std::fs::write(&tree_path, "(A:0.1,B:0.1)root;").wrap_err("When writing the tree fixture")?;
      std::fs::write(&fasta_path, ">A sample description\nACGT\n>B\nACGT\n")
        .wrap_err("When writing the alignment fixture")?;

      let args = TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw {
        alignment: AlignmentArgs {
          alignment: vec![fasta_path],
        },
        tree: Some(tree_path),
        method_anc: method,
        dense,
        model_args: ModelArgs {
          model: GtrModelNameCli::JC69,
          ..ModelArgs::default()
        },
        include_leaves,
        output_reconstructed_nuc_fasta: Some(out_path.clone()),
        ..TreetimeAncestralArgsRaw::default()
      })?;

      run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;

      Ok(
        fasta_read_file(out_path, &Alphabet::default())?
          .into_iter()
          .map(|record| (record.seq_name, record.desc))
          .collect(),
      )
    }
  }
}

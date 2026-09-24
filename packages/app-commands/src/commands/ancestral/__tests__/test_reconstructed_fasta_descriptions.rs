#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::ancestral::run::run_ancestral_reconstruction;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::method_anc::MethodAncestralCli;
  use crate::commands::shared::model::GtrModelNameCli;
  use crate::commands::shared::model::ModelArgs;
  use eyre::{Report, WrapErr};
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use tempfile::tempdir;
  use treetime::alphabet::alphabet::Alphabet;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_io::fasta::read_many_fasta_path;

  fn reconstructed_descriptions(
    method: MethodAncestralCli,
    dense: Option<bool>,
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
      include_leaves: true,
      output_reconstructed_nuc_fasta: Some(out_path.clone()),
      ..TreetimeAncestralArgsRaw::default()
    })?;

    run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress)?;

    Ok(
      read_many_fasta_path(&[out_path], &Alphabet::default())?
        .into_iter()
        .map(|record| (record.seq_name, record.desc))
        .collect(),
    )
  }

  #[test]
  fn test_reconstructed_fasta_descriptions_parsimony() -> Result<(), Report> {
    let descriptions = reconstructed_descriptions(MethodAncestralCli::Parsimony, None)?;
    assert_eq!(Some(&Some("sample description".to_owned())), descriptions.get("A"));
    assert_eq!(Some(&None), descriptions.get("B"));
    assert_eq!(Some(&None), descriptions.get("root"));
    Ok(())
  }

  #[test]
  fn test_reconstructed_fasta_descriptions_marginal() -> Result<(), Report> {
    let descriptions = reconstructed_descriptions(MethodAncestralCli::Marginal, Some(false))?;
    assert_eq!(Some(&Some("sample description".to_owned())), descriptions.get("A"));
    assert_eq!(Some(&None), descriptions.get("B"));
    assert_eq!(Some(&None), descriptions.get("root"));
    Ok(())
  }
}

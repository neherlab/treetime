#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::ancestral::run::run_ancestral_reconstruction;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::method_anc::MethodAncestralCli;
  use crate::commands::shared::model::GtrModelNameCli;
  use crate::commands::shared::model::ModelArgs;
  use crate::commands::shared::representation::RepresentationCli;
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
    let descriptions = reconstructed_descriptions(MethodAncestralCli::Parsimony, None, None, true)?;
    assert_eq!(Some(&Some("sample description".to_owned())), descriptions.get("A"));
    assert_eq!(Some(&None), descriptions.get("B"));
    assert_eq!(Some(&None), descriptions.get("root"));
    Ok(())
  }

  #[test]
  fn test_reconstructed_fasta_descriptions_marginal() -> Result<(), Report> {
    let descriptions = reconstructed_descriptions(
      MethodAncestralCli::Marginal,
      Some(RepresentationCli::Sparse),
      Some(GtrModelNameCli::JC69),
      true,
    )?;
    assert_eq!(Some(&Some("sample description".to_owned())), descriptions.get("A"));
    assert_eq!(Some(&None), descriptions.get("B"));
    assert_eq!(Some(&None), descriptions.get("root"));
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::parsimony(      (MethodAncestralCli::Parsimony, None,                            None))]
  #[case::marginal_sparse((MethodAncestralCli::Marginal,  Some(RepresentationCli::Sparse), Some(GtrModelNameCli::JC69)))]
  #[case::marginal_dense( (MethodAncestralCli::Marginal,  Some(RepresentationCli::Dense),  Some(GtrModelNameCli::JC69)))]
  #[trace]
  fn test_reconstructed_fasta_without_leaves_holds_only_internal_nodes(
    #[case] (method, representation, model): (MethodAncestralCli, Option<RepresentationCli>, Option<GtrModelNameCli>),
  ) -> Result<(), Report> {
    let descriptions = reconstructed_descriptions(method, representation, model, false)?;
    let expected = btreemap! { o!("root") => None };
    assert_eq!(expected, descriptions);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::commands::shared::method_anc::MethodAncArgs;
    use crate::commands::shared::reconstructed_fasta::ReconstructedFastaOutputArgs;
    use crate::commands::shared::representation::RepresentationArgs;
    use crate::commands::shared::tip_sequences::TipSequenceArgs;
    use crate::commands::shared::tree_input::TreeArgs;

    pub(super) fn reconstructed_descriptions(
      method: MethodAncestralCli,
      representation: Option<RepresentationCli>,
      model: Option<GtrModelNameCli>,
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
        tree: TreeArgs { tree: Some(tree_path) },
        method_anc: MethodAncArgs { method_anc: method },
        representation: RepresentationArgs { representation },
        model_args: ModelArgs {
          model,
          ..ModelArgs::default()
        },
        tip_sequences: TipSequenceArgs {
          include_leaves,
          ..TipSequenceArgs::default()
        },
        reconstructed_fasta_output: ReconstructedFastaOutputArgs {
          output_reconstructed_nuc_fasta: Some(out_path.clone()),
        },
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

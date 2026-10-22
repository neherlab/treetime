#[cfg(test)]
mod tests {
  use crate::commands::shared::ambiguity::AncestralAmbiguityOutput;
  use eyre::Report;
  use helpers::{AmbiguityRun, owned, run};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;

  const REPORTED: &[(&str, &str)] = &[("AB", "ARYM"), ("CD", "ARYK"), ("root", "ARYN")];
  const MOST_LIKELY: &[(&str, &str)] = &[("AB", "AACA"), ("CD", "AACG"), ("root", "AACA")];

  #[rustfmt::skip]
  #[rstest]
  #[case::every_output(   vec![],                                                (REPORTED,    "ARYN"))]
  #[case::fasta_only(     vec![AncestralAmbiguityOutput::ReconstructedNucFasta], (REPORTED,    "AACA"))]
  #[case::node_data_only( vec![AncestralAmbiguityOutput::AugurNodeData],         (MOST_LIKELY, "ARYN"))]
  #[trace]
  fn test_ambiguity_report_reaches_only_the_chosen_sequence_outputs(
    #[case] outputs: Vec<AncestralAmbiguityOutput>,
    #[case] (expected_fasta, expected_node_data_root): (&[(&str, &str)], &str),
  ) -> Result<(), Report> {
    let AmbiguityRun { fasta, node_data } = run(Some(0.0), outputs)?;

    assert_eq!(owned(&expected_fasta.iter().copied().collect()), fasta);
    assert_eq!(expected_node_data_root, node_data["root"]);
    Ok(())
  }

  #[test]
  fn test_ambiguity_report_is_off_without_a_threshold() -> Result<(), Report> {
    let AmbiguityRun { fasta, node_data } = run(None, vec![])?;

    let expected: BTreeMap<&str, &str> = btreemap! { "AB" => "AACA", "CD" => "AACG", "root" => "AACA" };
    assert_eq!(owned(&expected), fasta);
    assert_eq!("AACA", node_data["root"]);
    Ok(())
  }

  mod helpers {
    use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
    use crate::commands::ancestral::run::run_ancestral_reconstruction;
    use crate::commands::shared::alignment::AlignmentArgs;
    use crate::commands::shared::ambiguity::{AncestralAmbiguityArgs, AncestralAmbiguityOutput};
    use crate::commands::shared::model::{GtrModelNameCli, ModelArgs};
    use crate::commands::shared::reconstructed_fasta::ReconstructedFastaOutputArgs;
    use crate::commands::shared::tree_input::TreeArgs;
    use eyre::{Report, WrapErr};
    use std::collections::BTreeMap;
    use std::fs;
    use tempfile::tempdir;
    use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
    use treetime::cancel::NoopCancel;
    use treetime::partition::marginal::posterior_states::AmbiguityThreshold;
    use treetime::progress::NoopProgress;
    use treetime_io::fasta::fasta_read_file;
    use treetime_utils::io::json::json_read_file;
    use util_augur_node_data_json::AugurNodeDataJsonAncestral;

    const TREE: &str = "((A:0.1,B:0.1)AB:0.1,(C:0.1,D:0.1)CD:0.1)root;";
    const ALIGNMENT: &str = ">A\nAACA\n>B\nAGTC\n>C\nAACG\n>D\nAGTT\n";

    pub(crate) struct AmbiguityRun {
      pub(crate) fasta: BTreeMap<String, String>,
      pub(crate) node_data: BTreeMap<String, String>,
    }

    pub(crate) fn run(threshold: Option<f64>, outputs: Vec<AncestralAmbiguityOutput>) -> Result<AmbiguityRun, Report> {
      let dir = tempdir().wrap_err("When creating a temporary directory")?;
      let tree_path = dir.path().join("tree.nwk");
      let fasta_path = dir.path().join("aln.fasta");
      let output_fasta = dir.path().join("reconstructed.fasta");
      let output_node_data = dir.path().join("node-data.json");
      fs::write(&tree_path, TREE).wrap_err("When writing the tree fixture")?;
      fs::write(&fasta_path, ALIGNMENT).wrap_err("When writing the alignment fixture")?;
      let args = TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw {
        alignment: AlignmentArgs {
          alignment: vec![fasta_path],
        },
        tree: TreeArgs { tree: Some(tree_path) },
        model_args: ModelArgs {
          model: GtrModelNameCli::JC69,
          ..ModelArgs::default()
        },
        reconstructed_fasta_output: ReconstructedFastaOutputArgs {
          output_reconstructed_nuc_fasta: Some(output_fasta.clone()),
        },
        output_augur_node_data: Some(output_node_data.clone()),
        ambiguity: AncestralAmbiguityArgs {
          ambiguity_threshold: threshold.map(AmbiguityThreshold::new).transpose()?,
          ambiguity_outputs: outputs,
        },
        ..TreetimeAncestralArgsRaw::default()
      })?;

      run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;

      let fasta = fasta_read_file(&output_fasta, &Alphabet::new(AlphabetName::Nuc)?)?
        .into_iter()
        .map(|record| (record.seq_name, record.seq.to_string()))
        .collect();
      let node_data: AugurNodeDataJsonAncestral = json_read_file(&output_node_data)?;
      let node_data = node_data
        .nodes
        .into_iter()
        .filter_map(|(name, node)| node.sequence.map(|sequence| (name, sequence)))
        .collect();
      Ok(AmbiguityRun { fasta, node_data })
    }

    pub(crate) fn owned(expected: &BTreeMap<&str, &str>) -> BTreeMap<String, String> {
      expected
        .iter()
        .map(|(&name, &seq)| (name.to_owned(), seq.to_owned()))
        .collect()
    }
  }
}

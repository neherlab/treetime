#[cfg(test)]
mod tests {
  use crate::commands::shared::ambiguity::MugrationAmbiguityOutput;
  use eyre::Report;
  use helpers::{MugrationRun, run};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::every_output( vec![],                                     ("node,country,country_ambiguous\nA,usa,\nB,germany,\nroot,germany,germany|usa\n", Some("germany|usa")))]
  #[case::csv_only(     vec![MugrationAmbiguityOutput::Traits],     ("node,country,country_ambiguous\nA,usa,\nB,germany,\nroot,germany,germany|usa\n", None))]
  #[case::node_data_only(vec![MugrationAmbiguityOutput::NodeData], ("node,country\nA,usa\nB,germany\nroot,germany\n",                                Some("germany|usa")))]
  #[trace]
  fn test_ambiguity_report_names_plausible_trait_states_in_the_chosen_outputs(
    #[case] outputs: Vec<MugrationAmbiguityOutput>,
    #[case] (expected_csv, expected_root_field): (&str, Option<&str>),
  ) -> Result<(), Report> {
    let MugrationRun { traits_csv, root_ambiguous } = run(outputs)?;

    assert_eq!(expected_csv, traits_csv);
    assert_eq!(expected_root_field.map(str::to_owned), root_ambiguous);
    Ok(())
  }

  mod helpers {
    use crate::commands::mugration::args::{TreetimeMugrationArgs, TreetimeMugrationArgsRaw};
    use crate::commands::mugration::run::run_mugration;
    use crate::commands::shared::ambiguity::{MugrationAmbiguityArgs, MugrationAmbiguityOutput};
    use crate::commands::shared::metadata::MetadataArgs;
    use crate::commands::shared::tree_input::TreeArgs;
    use eyre::{Report, WrapErr};
    use std::fs;
    use tempfile::tempdir;
    use treetime::cancel::NoopCancel;
    use treetime::partition::marginal::posterior_states::AmbiguityThreshold;
    use treetime::progress::NoopProgress;
    use treetime_utils::io::json::json_read_file;
    use treetime_utils::o;
    use util_augur_node_data_json::AugurNodeDataJsonTraits;

    const TREE: &str = "(A:0.1,B:0.1)root;";
    const METADATA: &str = "strain\tcountry\nA\tusa\nB\tgermany\n";

    pub(super) struct MugrationRun {
      pub(super) traits_csv: String,
      pub(super) root_ambiguous: Option<String>,
    }

    pub(super) fn run(outputs: Vec<MugrationAmbiguityOutput>) -> Result<MugrationRun, Report> {
      let dir = tempdir().wrap_err("When creating a temporary directory")?;
      let tree_path = dir.path().join("tree.nwk");
      let metadata_path = dir.path().join("metadata.tsv");
      let traits_csv = dir.path().join("traits.csv");
      let node_data = dir.path().join("node-data.json");
      fs::write(&tree_path, TREE).wrap_err("When writing the tree fixture")?;
      fs::write(&metadata_path, METADATA).wrap_err("When writing the metadata fixture")?;
      let args = TreetimeMugrationArgs::try_from(TreetimeMugrationArgsRaw {
        tree: TreeArgs { tree: Some(tree_path) },
        metadata: MetadataArgs {
          metadata: Some(metadata_path),
        },
        attribute: Some(o!("country")),
        output_traits: Some(traits_csv.clone()),
        output_node_data: Some(node_data.clone()),
        ambiguity: MugrationAmbiguityArgs {
          ambiguity_threshold: Some(AmbiguityThreshold::new(0.5)?),
          ambiguity_outputs: outputs,
        },
        ..TreetimeMugrationArgsRaw::default()
      })?;

      run_mugration(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;

      let node_data: AugurNodeDataJsonTraits = json_read_file(&node_data)?;
      let root_ambiguous = node_data.nodes["root"]
        .fields
        .get("country_ambiguous")
        .and_then(|value| value.as_str())
        .map(str::to_owned);
      Ok(MugrationRun {
        traits_csv: fs::read_to_string(&traits_csv).wrap_err("When reading the traits CSV")?,
        root_ambiguous,
      })
    }
  }
}

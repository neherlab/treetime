#[cfg(test)]
mod tests {
  use crate::commands::mugration::run::{MugrationTraits, pair_traits};
  use app_output::annotated_graph::{AnnotatedGraph, AnnotatedTreeView, Divergence, TreeTraits};
  use app_output::output_plan::{CommandKind, TreeWriteKind};
  use app_output::tree_output::write_tree_outputs;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::path::Path;
  use tempfile::TempDir;
  use treetime::cancel::NoopCancel;
  use treetime::mugration::pipeline::{self, MugrationInput, MugrationParams};
  use treetime::progress::NoopProgress;
  use treetime_io::nwk::{NwkStyle, nwk_read};
  use treetime_utils::io::fs::read_file_to_string;
  use treetime_utils::o;

  #[test]
  fn test_mugration_annotated_tree_has_trait_comments() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2)root;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };
    let params = MugrationParams {
      missing_data: o!("?"),
      pc: None,
      missing_weights_threshold: 0.5,
      iterations: 5,
      sampling_bias_correction: None,
      smooth_initial_pi: false,
      filter_uninformative_root: false,
    };
    let MugrationTraits {
      traits,
      observed_values,
    } = pair_traits(
      traits.into_iter().collect(),
      &graph,
      &names,
      Path::new("metadata.tsv"),
      &NoopProgress,
    );
    let input = MugrationInput {
      graph,
      traits,
      observed_values,
      weights: None,
      branch_lengths: branch_lengths.clone(),
    };
    let output = pipeline::run(&params, input, &names, &NoopCancel, &NoopProgress).map_err(|err| err.into_report())?;
    let annotated = AnnotatedGraph {
      graph: &output.graph,
      names: &names,
      divergence_branch_lengths: &branch_lengths,
      time_branch_lengths: None,
      divergence: Divergence::CumulativeBranchLength,
      sequences: None,
      dates: None,
      traits: Some(TreeTraits {
        attribute: "country",
        states: &output.states,
        values: &output.reconstructed_traits,
        profiles: &output.confidences,
      }),
    };
    let dir = TempDir::new()?;
    let path = dir.path().join("mugration.nexus");

    write_tree_outputs(
      &AnnotatedTreeView::new(&annotated)?,
      &btreemap! { TreeWriteKind::Nexus(NwkStyle::Beast) => path.clone() },
      CommandKind::Mugration,
      &NoopProgress,
    )?;

    let actual = read_file_to_string(&path)?;
    let expected = indoc! {r#"
      #NEXUS

      Begin Taxa;
        Dimensions NTax=2;
        TaxLabels
          A
          B
        ;
      End;

      Begin Trees;
        Tree tree1 = (A[&country="usa"]:0.1,B[&country="germany"]:0.2)root[&country="usa"];
      End;
    "#};
    assert_eq!(expected, actual);

    Ok(())
  }
}

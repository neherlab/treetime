#[cfg(test)]
mod tests {
  use app_output::discrete_trait_comment::DiscreteTraitCommentProvider;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime::cancel::NoopCancel;
  use treetime::mugration::pipeline::{self, MugrationInput, MugrationParams};
  use treetime_io::nex::NexWriteOptions;
  use treetime_io::nwk::{CommentProviders, NwkStyle, nwk_read_str};
  use treetime_utils::o;

  #[test]
  fn test_mugration_annotated_tree_has_trait_comments() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2)root;")?;
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
    let input = MugrationInput {
      graph,
      traits,
      weights: None,
      branch_lengths: branch_lengths.clone(),
    };
    let output = pipeline::run(&params, input, &names, &NoopCancel).map_err(|err| err.into_report())?;
    let provider = DiscreteTraitCommentProvider::new(&output.reconstructed_traits, "country");
    let providers = CommentProviders::new().with(&provider);

    let options = NexWriteOptions {
      style: NwkStyle::Beast,
      ..NexWriteOptions::default()
    };
    let actual = treetime_io::nex::nex_write_str_with(&output.graph, &names, &branch_lengths, &options, &providers)?;
    let expected = indoc! {r#"
      #NEXUS
      Begin Taxa;
        Dimensions NTax=2;
        TaxLabels A B;
      End;
      Begin Trees;
        Tree tree1=(A[&country=usa]:0.1,B[&country=germany]:0.2)root[&country=usa];
      End;

    "#};
    assert_eq!(expected, actual);

    Ok(())
  }
}

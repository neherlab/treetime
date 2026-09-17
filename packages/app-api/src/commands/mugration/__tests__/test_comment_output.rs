#[cfg(test)]
mod tests {
  use app_output::discrete_trait_comment::DiscreteTraitCommentProvider;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime::cancel::NoopCancel;
  use treetime::mugration::mugration::execute_mugration;
  use treetime_io::nex::NexWriteOptions;
  use treetime_io::nwk::{CommentProviders, NwkStyle, nwk_read_str};
  use treetime_utils::o;

  #[test]
  fn test_mugration_annotated_tree_has_trait_comments() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };
    let names_tt_1 = names.clone();
    let (result, maps) = execute_mugration(
      graph,
      &confidences,
      &names_tt_1,
      &branch_lengths,
      &traits,
      "country",
      None,
      "?",
      None,
      0.5,
      5,
      None,
      false,
      false,
      &NoopCancel,
    )?;
    let provider = DiscreteTraitCommentProvider::new(&maps.reconstructed_traits, &result.traits.attribute);
    let providers = CommentProviders::new().with(&provider);

    let options = NexWriteOptions {
      style: NwkStyle::Beast,
      ..NexWriteOptions::default()
    };
    let actual = treetime_io::nex::nex_write_str_with(&result.graph, &names, &branch_lengths, &options, &providers)?;
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

#[cfg(test)]
mod tests {
  use crate::mugration::input::MugrationInput;
  use crate::mugration::mugration::execute_mugration;
  use crate::partition::marginal_discrete::DiscreteCommentProvider;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_io::nex::{NexWriteOptions, nex_write_str_with};
  use treetime_io::nwk::{CommentProviders, nwk_read_str};
  use treetime_utils::o;

  #[test]
  fn test_mugration_annotated_tree_has_trait_comments() -> Result<(), Report> {
    let graph = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };
    let input = MugrationInput {
      graph,
      traits,
      attribute: o!("country"),
      weights: None,
      missing_data: o!("?"),
      pc: None,
      missing_weights_threshold: 0.5,
      iterations: 5,
      sampling_bias_correction: None,
    };
    let result = execute_mugration(input)?;
    let provider = DiscreteCommentProvider::new(&result.partition, &result.traits.attribute);
    let providers = CommentProviders::new().with(&provider);

    let actual = nex_write_str_with(&result.graph, &NexWriteOptions::default(), &providers)?;
    let expected = indoc! {r#"
      #NEXUS
      Begin Taxa;
        Dimensions NTax=2;
        TaxLabels A B;
      End;
      Begin Trees;
        Tree tree1=(A:0.1[&country="usa"],B:0.2[&country="germany"])root[&country="usa"];;
      End;

    "#};
    assert_eq!(expected, actual);

    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use crate::commands::mugration::comment_provider::PartitionCommentProvider;
  use crate::commands::mugration::input::MugrationInput;
  use crate::commands::mugration::run::execute_mugration;
  use eyre::Report;
  use maplit::btreemap;
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
    };
    let result = execute_mugration(input)?;
    let provider = PartitionCommentProvider::new(&result.partition, &result.traits.attribute);
    let providers = CommentProviders::new().with(&provider);

    let actual = nex_write_str_with(&result.graph, &NexWriteOptions::default(), &providers)?;
    assert!(actual.starts_with("#NEXUS\nBegin Taxa;\n"));
    assert!(actual.contains("  Dimensions NTax=2;\n"));
    assert!(actual.contains("  TaxLabels A B;\n"));
    assert!(actual.contains(
      "  Tree tree1=(A:0.1[&country=\"usa\"],B:0.2[&country=\"germany\"])root[&country=\"usa\"];;\n"
    ));
    assert!(actual.ends_with("End;\n\n"));

    Ok(())
  }
}

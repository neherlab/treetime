#[cfg(test)]
mod tests {
  use crate::mugration::mugration::execute_mugration;
  use crate::partition::marginal_discrete::DiscreteCommentProvider;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_io::nex::NexWriteOptions;
  use treetime_io::nwk::{CommentProviders, NwkStyle, nwk_read_str};
  use treetime_utils::o;

  #[test]
  fn test_mugration_annotated_tree_has_trait_comments() -> Result<(), Report> {
    let graph = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };
    let result = execute_mugration(graph, &traits, "country", None, "?", None, 0.5, 5, None, false, false)?;
    let provider = DiscreteCommentProvider::new(&result.partition, &result.traits.attribute);
    let providers = CommentProviders::new().with(&provider);

    let options = NexWriteOptions {
      style: NwkStyle::Beast,
      ..NexWriteOptions::default()
    };
    let actual = treetime_io::nex::nex_write_str_with(&result.graph, &options, &providers)?;
    let expected = indoc! {r#"
      #NEXUS
      Begin Taxa;
        Dimensions NTax=2;
        TaxLabels A B;
      End;
      Begin Trees;
        Tree tree1=(A[&country=usa]:0.1,B[&country=germany]:0.2)root[&country=usa];;
      End;

    "#};
    assert_eq!(expected, actual);

    Ok(())
  }
}

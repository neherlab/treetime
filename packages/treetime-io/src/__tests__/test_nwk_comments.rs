#[cfg(test)]
mod tests {
  use crate::nwk::{NwkNodeComments, NwkParse, NwkStyle, NwkWriteOptions, nwk_read, nwk_write_str};
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::tree_view::TreeView;
  use treetime_utils::o;

  #[rustfmt::skip]
  #[rstest]
  #[case::beast_single(     (NwkStyle::Beast, vec![("country", "usa")]),                  "((A[&country=usa]:0.1,B:0.2)inner:0.3,C:0.4)root;")]
  #[case::beast_two_keys(   (NwkStyle::Beast, vec![("country", "usa"), ("region", "na")]), "((A[&country=usa,region=na]:0.1,B:0.2)inner:0.3,C:0.4)root;")]
  #[case::nhx(              (NwkStyle::Nhx,   vec![("S", "human"), ("D", "Y")]),           "((A[&&NHX:D=Y:S=human]:0.1,B:0.2)inner:0.3,C:0.4)root;")]
  #[case::plain_suppresses( (NwkStyle::Plain, vec![("country", "usa")]),                  "((A:0.1,B:0.2)inner:0.3,C:0.4)root;")]
  #[case::beast_number(     (NwkStyle::Beast, vec![("date", "2020.50")]),                 "((A[&date=2020.5]:0.1,B:0.2)inner:0.3,C:0.4)root;")]
  #[case::beast_special(    (NwkStyle::Beast, vec![("label", "New York, USA")]),          r#"((A[&label="New York, USA"]:0.1,B:0.2)inner:0.3,C:0.4)root;"#)]
  #[case::beast_no_comments((NwkStyle::Beast, vec![]),                                     "((A:0.1,B:0.2)inner:0.3,C:0.4)root;")]
  #[trace]
  fn test_nwk_comments_written_on_their_node(
    #[case] (style, comments): (NwkStyle, Vec<(&str, &str)>),
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let parse = nwk_read(b"((A:0.1,B:0.2)inner:0.3,C:0.4)root;".as_slice())?;
    let names = parse.names();
    let node_a = helpers::find_node_key_by_name(&names, "A");
    let comments: NwkNodeComments = btreemap! {
      node_a => comments.into_iter().map(|(key, value)| (o!(key), o!(value))).collect(),
    };
    let NwkParse {
      graph, branch_lengths, ..
    } = parse;
    let options = NwkWriteOptions {
      style,
      ..NwkWriteOptions::default()
    };

    let actual = nwk_write_str(&TreeView::new(&graph)?, &names, &branch_lengths, &options, &comments)?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_nwk_comments_name_quoting_special_chars() -> Result<(), Report> {
    let parse = nwk_read(b"('node (1)':0.1,B:0.2)root;".as_slice())?;
    let names = parse.names();
    let NwkParse {
      graph, branch_lengths, ..
    } = parse;

    let actual = nwk_write_str(
      &TreeView::new(&graph)?,
      &names,
      &branch_lengths,
      &NwkWriteOptions::default(),
      &btreemap! {},
    )?;

    assert_eq!("('node (1)':0.1,B:0.2)root;", actual);
    Ok(())
  }

  mod helpers {
    use std::collections::BTreeMap;
    use treetime_graph::node::GraphNodeKey;

    pub(super) fn find_node_key_by_name(names: &BTreeMap<GraphNodeKey, Option<String>>, name: &str) -> GraphNodeKey {
      names
        .iter()
        .find_map(|(key, node_name)| (node_name.as_deref() == Some(name)).then_some(*key))
        .unwrap_or_else(|| panic!("Missing test node '{name}'"))
    }
  }
}

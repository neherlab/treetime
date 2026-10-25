#[cfg(test)]
mod tests {
  use crate::nwk::{NewickValue, NwkParse, NwkStyle, NwkWriteOptions, nwk_read, nwk_write_str};
  use eyre::Report;
  use treetime_graph::keymap;

  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::schedule::GraphSchedule;
  use treetime_utils::o;

  #[rustfmt::skip]
  #[rstest]
  #[case::beast_single(     (NwkStyle::Beast, vec![("country", "usa")]),                  r#"((A[&country="usa"]:0.1,B:0.2)inner:0.3,C:0.4)root;"#)]
  #[case::beast_caller_order((NwkStyle::Beast, vec![("region", "na"), ("country", "usa")]), r#"((A[&region="na",country="usa"]:0.1,B:0.2)inner:0.3,C:0.4)root;"#)]
  #[case::nhx(              (NwkStyle::Nhx,   vec![("S", "human"), ("D", "Y")]),           "((A[&&NHX:S=human:D=Y]:0.1,B:0.2)inner:0.3,C:0.4)root;")]
  #[case::plain_suppresses( (NwkStyle::Plain, vec![("country", "usa")]),                  "((A:0.1,B:0.2)inner:0.3,C:0.4)root;")]
  #[case::beast_special(    (NwkStyle::Beast, vec![("label", "New York, USA")]),          r#"((A[&label="New York, USA"]:0.1,B:0.2)inner:0.3,C:0.4)root;"#)]
  #[case::beast_no_comments((NwkStyle::Beast, vec![]),                                     "((A:0.1,B:0.2)inner:0.3,C:0.4)root;")]
  #[trace]
  fn test_nwk_comments_written_on_their_node(
    #[case] (style, comments): (NwkStyle, Vec<(&str, &str)>),
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let comments = comments
      .into_iter()
      .map(|(key, value)| (o!(key), NewickValue::String(o!(value))))
      .collect();

    let actual = helpers::write_with_comments_on_a(style, comments)?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_nwk_beast_numeric_value_written_bare() -> Result<(), Report> {
    let comments = vec![(o!("date"), NewickValue::NumberText(o!("2020.50")))];

    let actual = helpers::write_with_comments_on_a(NwkStyle::Beast, comments)?;

    assert_eq!("((A[&date=2020.50]:0.1,B:0.2)inner:0.3,C:0.4)root;", actual);
    Ok(())
  }

  #[test]
  fn test_nwk_comments_name_quoting_special_chars() -> Result<(), Report> {
    let parse = nwk_read(b"('node (1)':0.1,B:0.2)root;".as_slice())?;
    let names = parse.names.clone();
    let NwkParse {
      graph, branch_lengths, ..
    } = parse;

    let actual = nwk_write_str(
      &GraphSchedule::new(&graph)?.into_tree()?,
      &names,
      &branch_lengths,
      &NwkWriteOptions::default(),
      &keymap! {},
    )?;

    assert_eq!("('node (1)':0.1,B:0.2)root;", actual);
    Ok(())
  }

  mod helpers {
    use crate::nwk::{NewickValue, NwkNodeComments, NwkParse, NwkStyle, NwkWriteOptions, nwk_read, nwk_write_str};
    use eyre::Report;
    use treetime_graph::keymap;

    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::node_names::NodeNames;
    use treetime_graph::schedule::GraphSchedule;

    pub(super) fn write_with_comments_on_a(
      style: NwkStyle,
      comments: Vec<(String, NewickValue)>,
    ) -> Result<String, Report> {
      let parse = nwk_read(b"((A:0.1,B:0.2)inner:0.3,C:0.4)root;".as_slice())?;
      let NwkParse {
        graph,
        names,
        branch_lengths,
        ..
      } = parse;
      let comments: NwkNodeComments = keymap! { find_node_key_by_name(&graph, &names, "A") => comments };
      let options = NwkWriteOptions {
        style,
        ..NwkWriteOptions::default()
      };
      nwk_write_str(
        &GraphSchedule::new(&graph)?.into_tree()?,
        &names,
        &branch_lengths,
        &options,
        &comments,
      )
    }

    fn find_node_key_by_name(graph: &Graph, names: &NodeNames, name: &str) -> GraphNodeKey {
      graph
        .get_nodes()
        .map(|node| node.key())
        .find(|key| &names[*key] == name)
        .unwrap_or_else(|| panic!("Missing test node '{name}'"))
    }
  }
}

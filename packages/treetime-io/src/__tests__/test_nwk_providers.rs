#[cfg(test)]
mod tests {
  use crate::nwk::{
    CommentProviders, NodeCommentProvider, NwkParse, NwkStyle, NwkWriteOptions, nwk_read_str, nwk_write_str,
    nwk_write_str_with,
  };
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;

  fn beast_options() -> NwkWriteOptions {
    NwkWriteOptions {
      style: NwkStyle::Beast,
      ..NwkWriteOptions::default()
    }
  }

  #[test]
  fn test_nwk_provider_empty_produces_same_output() -> Result<(), Report> {
    let (graph, names, branch_lengths) = helpers::make_graph()?;
    let providers = CommentProviders::new();

    let expected = nwk_write_str(&graph, &names, &branch_lengths, &beast_options())?;
    let actual = nwk_write_str_with(&graph, &names, &branch_lengths, &beast_options(), &providers)?;

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_provider_single_beast_annotation() -> Result<(), Report> {
    let (graph, names, branch_lengths) = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&names, "A");
    let provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "country".to_owned() => "usa".to_owned(),
      },
    });
    let providers = CommentProviders::new().with(&provider);

    let actual = nwk_write_str_with(&graph, &names, &branch_lengths, &beast_options(), &providers)?;
    let expected = "((A[&country=usa]:0.1,B:0.2)inner:0.3,C:0.4)root;";

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_provider_multiple_merged_into_single_block() -> Result<(), Report> {
    let (graph, names, branch_lengths) = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&names, "A");
    let country_provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "country".to_owned() => "usa".to_owned(),
      },
    });
    let region_provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "region".to_owned() => "na".to_owned(),
      },
    });
    let providers = CommentProviders::new().with(&country_provider).with(&region_provider);

    let actual = nwk_write_str_with(&graph, &names, &branch_lengths, &beast_options(), &providers)?;
    let expected = "((A[&country=usa,region=na]:0.1,B:0.2)inner:0.3,C:0.4)root;";

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_provider_later_overrides_earlier_on_key_collision() -> Result<(), Report> {
    let (graph, names, branch_lengths) = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&names, "A");
    let first_provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "country".to_owned() => "usa".to_owned(),
      },
    });
    let second_provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "country".to_owned() => "canada".to_owned(),
      },
    });
    let providers = CommentProviders::new().with(&first_provider).with(&second_provider);

    let actual = nwk_write_str_with(&graph, &names, &branch_lengths, &beast_options(), &providers)?;
    let expected = "((A[&country=canada]:0.1,B:0.2)inner:0.3,C:0.4)root;";

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_provider_multikey_partial_override() -> Result<(), Report> {
    let (graph, names, branch_lengths) = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&names, "A");
    let base_provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "country".to_owned() => "base".to_owned(),
        "note".to_owned() => "base".to_owned(),
      },
    });
    let override_provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "country".to_owned() => "usa".to_owned(),
      },
    });
    let providers = CommentProviders::new().with(&base_provider).with(&override_provider);

    let actual = nwk_write_str_with(&graph, &names, &branch_lengths, &beast_options(), &providers)?;
    let expected = "((A[&country=usa,note=base]:0.1,B:0.2)inner:0.3,C:0.4)root;";

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_plain_style_suppresses_all_annotations() -> Result<(), Report> {
    let (graph, names, branch_lengths) = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&names, "A");
    let provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "country".to_owned() => "usa".to_owned(),
      },
    });
    let providers = CommentProviders::new().with(&provider);

    let actual = nwk_write_str_with(&graph, &names, &branch_lengths, &NwkWriteOptions::default(), &providers)?;
    let expected = "((A:0.1,B:0.2)inner:0.3,C:0.4)root;";

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_nhx_style_annotation() -> Result<(), Report> {
    let (graph, names, branch_lengths) = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&names, "A");
    let provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "S".to_owned() => "human".to_owned(),
        "D".to_owned() => "Y".to_owned(),
      },
    });
    let providers = CommentProviders::new().with(&provider);

    let options = NwkWriteOptions {
      style: NwkStyle::Nhx,
      ..NwkWriteOptions::default()
    };
    let actual = nwk_write_str_with(&graph, &names, &branch_lengths, &options, &providers)?;
    let expected = "((A[&&NHX:D=Y:S=human]:0.1,B:0.2)inner:0.3,C:0.4)root;";

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_beast_numeric_value_written_bare() -> Result<(), Report> {
    let (graph, names, branch_lengths) = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&names, "A");
    let provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "date".to_owned() => "2020.50".to_owned(),
      },
    });
    let providers = CommentProviders::new().with(&provider);

    let actual = nwk_write_str_with(&graph, &names, &branch_lengths, &beast_options(), &providers)?;
    let expected = "((A[&date=2020.5]:0.1,B:0.2)inner:0.3,C:0.4)root;";

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_beast_quoted_value_with_special_chars() -> Result<(), Report> {
    let (graph, names, branch_lengths) = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&names, "A");
    let provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "label".to_owned() => "New York, USA".to_owned(),
      },
    });
    let providers = CommentProviders::new().with(&provider);

    let actual = nwk_write_str_with(&graph, &names, &branch_lengths, &beast_options(), &providers)?;
    assert!(actual.contains(r#""New York, USA""#));

    Ok(())
  }

  #[test]
  fn test_nwk_name_quoting_special_chars() -> Result<(), Report> {
    let parse = nwk_read_str("('node (1)':0.1,B:0.2)root;")?;
    let names = parse.names();
    let NwkParse {
      graph, branch_lengths, ..
    } = parse;

    let actual = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!("('node (1)':0.1,B:0.2)root;", actual);

    Ok(())
  }

  mod helpers {
    use super::*;
    use treetime_graph::edge::GraphEdgeKey;

    pub(super) struct MockCommentProvider {
      comments: BTreeMap<GraphNodeKey, BTreeMap<String, String>>,
    }

    impl MockCommentProvider {
      pub(super) fn new(comments: BTreeMap<GraphNodeKey, BTreeMap<String, String>>) -> Self {
        Self { comments }
      }
    }

    impl NodeCommentProvider for MockCommentProvider {
      fn node_comments(&self, key: GraphNodeKey) -> Result<BTreeMap<String, String>, Report> {
        Ok(self.comments.get(&key).cloned().unwrap_or_default())
      }
    }

    #[allow(clippy::type_complexity)]
    pub(super) fn make_graph() -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let parse = nwk_read_str("((A:0.1,B:0.2)inner:0.3,C:0.4)root;")?;
      let names = parse.names();
      let NwkParse {
        graph, branch_lengths, ..
      } = parse;
      Ok((graph, names, branch_lengths))
    }

    pub(super) fn find_node_key_by_name(names: &BTreeMap<GraphNodeKey, Option<String>>, name: &str) -> GraphNodeKey {
      names
        .iter()
        .find_map(|(key, node_name)| (node_name.as_deref() == Some(name)).then_some(*key))
        .unwrap_or_else(|| panic!("Missing test node '{name}'"))
    }
  }
}

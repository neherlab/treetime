#[cfg(test)]
mod tests {
  use crate::nwk::{
    CommentProviders, EdgeFromNwk, EdgeToNwk, NodeCommentProvider, NodeFromNwk, NodeToNwk, NwkStyle, NwkWriteOptions,
    nwk_read_str, nwk_write_str, nwk_write_str_with,
  };
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::{GraphEdge, HasBranchLength};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::{GraphNode, GraphNodeKey, Named};

  fn beast_options() -> NwkWriteOptions {
    NwkWriteOptions {
      style: NwkStyle::Beast,
      ..NwkWriteOptions::default()
    }
  }

  #[test]
  fn test_nwk_provider_empty_produces_same_output() -> Result<(), Report> {
    let graph = helpers::make_graph()?;
    let providers = CommentProviders::new();

    let expected = nwk_write_str(&graph, &beast_options())?;
    let actual = nwk_write_str_with(&graph, &beast_options(), &providers)?;

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_provider_single_beast_annotation() -> Result<(), Report> {
    let graph = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&graph, "A");
    let provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "country".to_owned() => "usa".to_owned(),
      },
    });
    let providers = CommentProviders::new().with(&provider);

    let actual = nwk_write_str_with(&graph, &beast_options(), &providers)?;
    let expected = "((A[&country=usa]:0.1,B:0.2)inner:0.3,C:0.4)root;";

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_provider_multiple_merged_into_single_block() -> Result<(), Report> {
    let graph = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&graph, "A");
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

    let actual = nwk_write_str_with(&graph, &beast_options(), &providers)?;
    let expected = "((A[&country=usa,region=na]:0.1,B:0.2)inner:0.3,C:0.4)root;";

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_provider_later_overrides_earlier_on_key_collision() -> Result<(), Report> {
    let graph = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&graph, "A");
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

    let actual = nwk_write_str_with(&graph, &beast_options(), &providers)?;
    let expected = "((A[&country=canada]:0.1,B:0.2)inner:0.3,C:0.4)root;";

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_provider_overrides_payload_comments() -> Result<(), Report> {
    let graph = helpers::make_graph()?;
    helpers::set_payload_comments(
      &graph,
      "A",
      btreemap! {
        "country".to_owned() => "payload".to_owned(),
        "note".to_owned() => "payload".to_owned(),
      },
    );
    let node_a = helpers::find_node_key_by_name(&graph, "A");
    let provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "country".to_owned() => "usa".to_owned(),
      },
    });
    let providers = CommentProviders::new().with(&provider);

    let actual = nwk_write_str_with(&graph, &beast_options(), &providers)?;
    let expected = "((A[&country=usa,note=payload]:0.1,B:0.2)inner:0.3,C:0.4)root;";

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_plain_style_suppresses_all_annotations() -> Result<(), Report> {
    let graph = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&graph, "A");
    let provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "country".to_owned() => "usa".to_owned(),
      },
    });
    let providers = CommentProviders::new().with(&provider);

    let actual = nwk_write_str_with(&graph, &NwkWriteOptions::default(), &providers)?;
    let expected = "((A:0.1,B:0.2)inner:0.3,C:0.4)root;";

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_nhx_style_annotation() -> Result<(), Report> {
    let graph = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&graph, "A");
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
    let actual = nwk_write_str_with(&graph, &options, &providers)?;
    let expected = "((A[&&NHX:D=Y:S=human]:0.1,B:0.2)inner:0.3,C:0.4)root;";

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_beast_numeric_value_written_bare() -> Result<(), Report> {
    let graph = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&graph, "A");
    let provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "date".to_owned() => "2020.50".to_owned(),
      },
    });
    let providers = CommentProviders::new().with(&provider);

    let actual = nwk_write_str_with(&graph, &beast_options(), &providers)?;
    let expected = "((A[&date=2020.5]:0.1,B:0.2)inner:0.3,C:0.4)root;";

    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_nwk_beast_quoted_value_with_special_chars() -> Result<(), Report> {
    let graph = helpers::make_graph()?;
    let node_a = helpers::find_node_key_by_name(&graph, "A");
    let provider = helpers::MockCommentProvider::new(btreemap! {
      node_a => btreemap! {
        "label".to_owned() => "New York, USA".to_owned(),
      },
    });
    let providers = CommentProviders::new().with(&provider);

    let actual = nwk_write_str_with(&graph, &beast_options(), &providers)?;
    assert!(actual.contains(r#""New York, USA""#));

    Ok(())
  }

  #[test]
  fn test_nwk_name_quoting_special_chars() -> Result<(), Report> {
    let graph: Graph<helpers::TestNode, helpers::TestEdge, ()> = nwk_read_str("('node (1)':0.1,B:0.2)root;")?;

    let actual = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!("('node (1)':0.1,B:0.2)root;", actual);

    Ok(())
  }

  mod helpers {
    use super::*;

    #[derive(Clone, Debug, PartialEq, Eq)]
    pub(super) struct TestNode {
      name: Option<String>,
      comments: BTreeMap<String, String>,
    }

    impl TestNode {
      fn new(name: Option<String>) -> Self {
        Self {
          name,
          comments: BTreeMap::new(),
        }
      }
    }

    impl GraphNode for TestNode {}

    impl Named for TestNode {
      fn name(&self) -> Option<impl AsRef<str>> {
        self.name.as_deref()
      }

      fn set_name(&mut self, name: Option<impl AsRef<str>>) {
        self.name = name.map(|name| name.as_ref().to_owned());
      }
    }

    impl NodeFromNwk for TestNode {
      fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
        Ok(Self::new(name.map(|name| name.as_ref().to_owned())))
      }
    }

    impl NodeToNwk for TestNode {
      fn nwk_name(&self) -> Option<impl AsRef<str>> {
        self.name.as_deref()
      }

      fn nwk_comments(&self) -> BTreeMap<String, String> {
        self.comments.clone()
      }
    }

    #[derive(Clone, Debug, PartialEq)]
    pub(super) struct TestEdge(Option<f64>);

    impl GraphEdge for TestEdge {}

    impl HasBranchLength for TestEdge {
      fn branch_length(&self) -> Option<f64> {
        self.0
      }

      fn set_branch_length(&mut self, weight: Option<f64>) {
        self.0 = weight;
      }
    }

    impl EdgeFromNwk for TestEdge {
      fn from_nwk(weight: Option<f64>) -> Result<Self, Report> {
        Ok(Self(weight))
      }
    }

    impl EdgeToNwk for TestEdge {
      fn nwk_weight(&self) -> Option<f64> {
        self.0
      }
    }

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

    pub(super) fn make_graph() -> Result<Graph<TestNode, TestEdge, ()>, Report> {
      nwk_read_str("((A:0.1,B:0.2)inner:0.3,C:0.4)root;")
    }

    pub(super) fn find_node_key_by_name(graph: &Graph<TestNode, TestEdge, ()>, name: &str) -> GraphNodeKey {
      graph
        .get_nodes()
        .iter()
        .find_map(|node| {
          let node = node.read_arc();
          let payload = node.payload().read_arc();
          (payload.name().map(|node_name| node_name.as_ref().to_owned()) == Some(name.to_owned())).then_some(node.key())
        })
        .unwrap_or_else(|| panic!("Missing test node '{name}'"))
    }

    pub(super) fn set_payload_comments(
      graph: &Graph<TestNode, TestEdge, ()>,
      name: &str,
      comments: BTreeMap<String, String>,
    ) {
      let node_key = find_node_key_by_name(graph, name);
      let node = graph
        .get_node(node_key)
        .unwrap_or_else(|| panic!("Missing test node '{name}'"));
      let node = node.read_arc();
      let mut payload = node.payload().write_arc();
      payload.comments = comments;
    }
  }
}

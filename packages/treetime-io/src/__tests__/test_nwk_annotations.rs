#[cfg(test)]
mod tests {
  use crate::nwk::{
    EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkStyle, NwkWriteOptions, nwk_read_str, nwk_write_str,
  };
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::{GraphEdge, HasBranchLength};
  use treetime_graph::node::{GraphNode, Named};

  #[derive(Clone, Debug, PartialEq, Eq)]
  struct AnnotNode {
    name: Option<String>,
    comments: BTreeMap<String, String>,
  }

  impl GraphNode for AnnotNode {}

  impl Named for AnnotNode {
    fn name(&self) -> Option<impl AsRef<str>> {
      self.name.as_deref()
    }

    fn set_name(&mut self, name: Option<impl AsRef<str>>) {
      self.name = name.map(|n| n.as_ref().to_owned());
    }
  }

  impl NodeFromNwk for AnnotNode {
    fn from_nwk(
      name: Option<impl AsRef<str>>,
      _confidence: Option<f64>,
      comments: &BTreeMap<String, String>,
    ) -> Result<Self, Report> {
      Ok(Self {
        name: name.map(|n| n.as_ref().to_owned()),
        comments: comments.clone(),
      })
    }
  }

  impl NodeToNwk for AnnotNode {
    fn nwk_name(&self) -> Option<impl AsRef<str>> {
      self.name.as_deref()
    }

    fn nwk_comments(&self) -> BTreeMap<String, String> {
      self.comments.clone()
    }
  }

  #[derive(Clone, Debug, PartialEq)]
  struct AnnotEdge(Option<f64>);

  impl GraphEdge for AnnotEdge {}

  impl HasBranchLength for AnnotEdge {
    fn branch_length(&self) -> Option<f64> {
      self.0
    }

    fn set_branch_length(&mut self, weight: Option<f64>) {
      self.0 = weight;
    }
  }

  impl EdgeFromNwk for AnnotEdge {
    fn from_nwk(weight: Option<f64>) -> Result<Self, Report> {
      Ok(Self(weight))
    }
  }

  impl EdgeToNwk for AnnotEdge {
    fn nwk_weight(&self) -> Option<f64> {
      self.0
    }
  }

  fn find_node_comments(
    graph: &treetime_graph::graph::Graph<AnnotNode, AnnotEdge, ()>,
    name: &str,
  ) -> BTreeMap<String, String> {
    graph
      .get_nodes()
      .iter()
      .find_map(|node| {
        let node = node.read_arc();
        let payload = node.payload().read_arc();
        payload
          .name()
          .is_some_and(|n| n.as_ref() == name)
          .then(|| payload.comments.clone())
      })
      .unwrap_or_else(|| panic!("Node '{name}' not found"))
  }

  #[test]
  fn test_nwk_annotations_beast_single_attribute() -> Result<(), Report> {
    let graph = nwk_read_str::<AnnotNode, AnnotEdge, ()>("(A[&prob=0.95]:0.1,B:0.2);")?;
    let comments = find_node_comments(&graph, "A");
    let expected = btreemap! { "prob".to_owned() => "0.95".to_owned() };
    assert_eq!(expected, comments);
    Ok(())
  }

  #[test]
  fn test_nwk_annotations_beast_multiple_attributes() -> Result<(), Report> {
    let graph = nwk_read_str::<AnnotNode, AnnotEdge, ()>(r#"(A[&prob=0.95,country="USA"]:0.1,B:0.2);"#)?;
    let comments = find_node_comments(&graph, "A");
    assert_eq!(comments.get("prob").map(String::as_str), Some("0.95"));
    assert_eq!(comments.get("country").map(String::as_str), Some("USA"));
    Ok(())
  }

  #[test]
  fn test_nwk_annotations_nhx_attributes() -> Result<(), Report> {
    let graph = nwk_read_str::<AnnotNode, AnnotEdge, ()>("(A[&&NHX:S=human:B=90]:0.1,B:0.2);")?;
    let comments = find_node_comments(&graph, "A");
    assert_eq!(comments.get("S").map(String::as_str), Some("human"));
    assert_eq!(comments.get("B").map(String::as_str), Some("90"));
    Ok(())
  }

  #[test]
  fn test_nwk_annotations_plain_comment_not_wired() -> Result<(), Report> {
    let graph = nwk_read_str::<AnnotNode, AnnotEdge, ()>("(A[a note]:0.1,B:0.2);")?;
    let comments = find_node_comments(&graph, "A");
    assert!(comments.is_empty());
    Ok(())
  }

  #[test]
  fn test_nwk_annotations_unannotated_tree_unchanged() -> Result<(), Report> {
    let graph = nwk_read_str::<AnnotNode, AnnotEdge, ()>("(A:0.1,B:0.2)root;")?;
    let comments_a = find_node_comments(&graph, "A");
    let comments_b = find_node_comments(&graph, "B");
    let comments_root = find_node_comments(&graph, "root");
    assert!(comments_a.is_empty());
    assert!(comments_b.is_empty());
    assert!(comments_root.is_empty());
    Ok(())
  }

  #[test]
  fn test_nwk_annotations_mixed_dialects_per_node() -> Result<(), Report> {
    let graph = nwk_read_str::<AnnotNode, AnnotEdge, ()>("(A[&k=v]:0.1,B[&&NHX:S=human]:0.2);")?;
    let comments_a = find_node_comments(&graph, "A");
    let comments_b = find_node_comments(&graph, "B");
    assert_eq!(comments_a.get("k").map(String::as_str), Some("v"));
    assert_eq!(comments_b.get("S").map(String::as_str), Some("human"));
    Ok(())
  }

  #[test]
  fn test_nwk_annotations_beast_boolean_value() -> Result<(), Report> {
    let graph = nwk_read_str::<AnnotNode, AnnotEdge, ()>("(A[&fixed=TRUE]:0.1,B:0.2);")?;
    let comments = find_node_comments(&graph, "A");
    assert_eq!(comments.get("fixed").map(String::as_str), Some("true"));
    Ok(())
  }

  #[test]
  fn test_nwk_annotations_beast_array_value() -> Result<(), Report> {
    let graph = nwk_read_str::<AnnotNode, AnnotEdge, ()>("(A[&hpd={1.0,2.0}]:0.1,B:0.2);")?;
    let comments = find_node_comments(&graph, "A");
    assert_eq!(comments.get("hpd").map(String::as_str), Some("{1,2}"));
    Ok(())
  }

  #[test]
  fn test_nwk_annotations_roundtrip_beast() -> Result<(), Report> {
    let input = "(A[&prob=0.95]:0.1,B:0.2)root;";
    let graph = nwk_read_str::<AnnotNode, AnnotEdge, ()>(input)?;
    let options = NwkWriteOptions {
      style: NwkStyle::Beast,
      ..NwkWriteOptions::default()
    };
    let output = nwk_write_str(&graph, &options)?;
    assert_eq!(input, output);
    Ok(())
  }

  #[test]
  fn test_nwk_annotations_root_annotated() -> Result<(), Report> {
    let graph = nwk_read_str::<AnnotNode, AnnotEdge, ()>("(A:0.1,B:0.2)root[&date=2020.5];")?;
    let comments = find_node_comments(&graph, "root");
    assert_eq!(comments.get("date").map(String::as_str), Some("2020.5"));
    Ok(())
  }
}

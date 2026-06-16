#[cfg(test)]
mod tests {
  use crate::parse::newick_from_string;
  use crate::types::{NewickEdgeData, NewickGraph, NewickNodeData, NewickValue, NewickWriteOptions, NwkStyle};
  use crate::write::newick_to_string;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;

  fn opts(style: NwkStyle) -> NewickWriteOptions {
    NewickWriteOptions {
      style,
      significant_digits: None,
      decimal_digits: None,
    }
  }

  // Plain style
  #[test]
  fn test_write_plain_simple() {
    let g = newick_from_string("(A:0.1,B:0.2)root;").unwrap();
    let s = newick_to_string(&g, &opts(NwkStyle::Plain)).unwrap();
    assert_eq!("(A:0.1,B:0.2)root;", s);
  }

  #[test]
  fn test_write_plain_strips_annotations() {
    let g = newick_from_string("(A[&prob=0.95]:0.1,B:0.2);").unwrap();
    let s = newick_to_string(&g, &opts(NwkStyle::Plain)).unwrap();
    assert_eq!("(A:0.1,B:0.2);", s);
  }

  #[test]
  fn test_write_plain_strips_raw_comments() {
    let g = newick_from_string("(A[some comment]:0.1,B:0.2);").unwrap();
    let s = newick_to_string(&g, &opts(NwkStyle::Plain)).unwrap();
    assert_eq!("(A:0.1,B:0.2);", s);
  }

  // Beast style
  #[test]
  fn test_write_beast_node_attrs() {
    let mut g = NewickGraph::new();
    let mut attrs = BTreeMap::new();
    attrs.insert("prob".to_owned(), NewickValue::Number(0.95));
    let root = g.add_node(NewickNodeData {
      name: None,
      confidence: None,
      node_attrs: BTreeMap::new(),
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    });
    let a = g.add_node(NewickNodeData {
      name: Some("A".to_owned()),
      confidence: None,
      node_attrs: attrs,
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    });
    let b = g.add_node(NewickNodeData::new().with_name("B"));
    g.add_edge(root, a, NewickEdgeData::new().with_length(0.1));
    g.add_edge(root, b, NewickEdgeData::new().with_length(0.2));
    g.root = root;

    let s = newick_to_string(&g, &opts(NwkStyle::Beast)).unwrap();
    assert_eq!("(A[&prob=0.95]:0.1,B:0.2);", s);
  }

  #[test]
  fn test_write_beast_branch_attrs() {
    let mut g = NewickGraph::new();
    let root = g.add_node(NewickNodeData::new());
    let a = g.add_node(NewickNodeData::new().with_name("A"));
    let mut branch_attrs = BTreeMap::new();
    branch_attrs.insert("rate".to_owned(), NewickValue::Number(1.5));
    g.add_edge(
      root,
      a,
      NewickEdgeData {
        branch_length: Some(0.1),
        branch_attrs,
        raw_comments: Vec::new(),
        is_acceptor: false,
      },
    );
    g.root = root;

    let s = newick_to_string(&g, &opts(NwkStyle::Beast)).unwrap();
    assert_eq!("(A:[&rate=1.5]0.1);", s);
  }

  #[test]
  fn test_write_beast_boolean_values() {
    let mut g = NewickGraph::new();
    let mut attrs = BTreeMap::new();
    attrs.insert("active".to_owned(), NewickValue::Boolean(true));
    attrs.insert("fixed".to_owned(), NewickValue::Boolean(false));
    let node = g.add_node(NewickNodeData {
      name: Some("A".to_owned()),
      confidence: None,
      node_attrs: attrs,
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    });
    g.root = node;

    let s = newick_to_string(&g, &opts(NwkStyle::Beast)).unwrap();
    assert_eq!("A[&active=TRUE,fixed=FALSE];", s);
  }

  #[test]
  fn test_write_beast_array_values() {
    let mut g = NewickGraph::new();
    let mut attrs = BTreeMap::new();
    attrs.insert(
      "hpd".to_owned(),
      NewickValue::Array(vec![NewickValue::Number(1.0), NewickValue::Number(2.0)]),
    );
    let node = g.add_node(NewickNodeData {
      name: Some("A".to_owned()),
      confidence: None,
      node_attrs: attrs,
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    });
    g.root = node;

    let s = newick_to_string(&g, &opts(NwkStyle::Beast)).unwrap();
    assert_eq!("A[&hpd={1,2}];", s);
  }

  #[test]
  fn test_write_beast_quoted_string() {
    let mut g = NewickGraph::new();
    let mut attrs = BTreeMap::new();
    attrs.insert("label".to_owned(), NewickValue::String("hello, world".to_owned()));
    let node = g.add_node(NewickNodeData {
      name: Some("A".to_owned()),
      confidence: None,
      node_attrs: attrs,
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    });
    g.root = node;

    let s = newick_to_string(&g, &opts(NwkStyle::Beast)).unwrap();
    assert_eq!("A[&label=\"hello, world\"];", s);
  }

  #[test]
  fn test_write_beast_raw_comments() {
    let mut g = NewickGraph::new();
    let node = g.add_node(NewickNodeData {
      name: Some("A".to_owned()),
      confidence: None,
      node_attrs: BTreeMap::new(),
      raw_comments: vec!["[some note]".to_owned()],
      hybrid: None,
      children: Vec::new(),
    });
    g.root = node;

    let s = newick_to_string(&g, &opts(NwkStyle::Beast)).unwrap();
    assert_eq!("A[some note];", s);
  }

  // NHX style
  #[test]
  fn test_write_nhx_node_attrs() {
    let mut g = NewickGraph::new();
    let mut attrs = BTreeMap::new();
    attrs.insert("S".to_owned(), NewickValue::String("human".to_owned()));
    attrs.insert("T".to_owned(), NewickValue::String("9606".to_owned()));
    let node = g.add_node(NewickNodeData {
      name: Some("A".to_owned()),
      confidence: None,
      node_attrs: attrs,
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    });
    g.root = node;

    let s = newick_to_string(&g, &opts(NwkStyle::Nhx)).unwrap();
    assert_eq!("A[&&NHX:S=human:T=9606];", s);
  }

  // Name quoting
  #[test]
  fn test_write_name_quoting_special_chars() {
    let mut g = NewickGraph::new();
    let node = g.add_node(NewickNodeData::new().with_name("node (1)"));
    g.root = node;
    let s = newick_to_string(&g, &opts(NwkStyle::Plain)).unwrap();
    assert_eq!("'node (1)';", s);
  }

  #[test]
  fn test_write_name_quoting_single_quote() {
    let mut g = NewickGraph::new();
    let node = g.add_node(NewickNodeData::new().with_name("it's"));
    g.root = node;
    let s = newick_to_string(&g, &opts(NwkStyle::Plain)).unwrap();
    assert_eq!("'it''s';", s);
  }

  #[test]
  fn test_write_name_no_quoting_needed() {
    let mut g = NewickGraph::new();
    let node = g.add_node(NewickNodeData::new().with_name("simple_name"));
    g.root = node;
    let s = newick_to_string(&g, &opts(NwkStyle::Plain)).unwrap();
    assert_eq!("simple_name;", s);
  }

  // Rooting prefix
  #[test]
  fn test_write_rooted_prefix() {
    let mut g = newick_from_string("(A,B);").unwrap();
    g.rooted = Some(true);
    let s = newick_to_string(&g, &opts(NwkStyle::Plain)).unwrap();
    assert_eq!("[&R](A,B);", s);
  }

  #[test]
  fn test_write_unrooted_prefix() {
    let mut g = newick_from_string("(A,B);").unwrap();
    g.rooted = Some(false);
    let s = newick_to_string(&g, &opts(NwkStyle::Plain)).unwrap();
    assert_eq!("[&U](A,B);", s);
  }

  // Float formatting
  #[test]
  fn test_write_significant_digits() {
    let g = newick_from_string("(A:0.123456789,B:0.2);").unwrap();

    // Full precision by default
    let full = newick_to_string(&g, &NewickWriteOptions::default()).unwrap();
    assert_eq!("(A:0.123456789,B:0.2);", full);

    // Explicit 3 significant digits
    let o = NewickWriteOptions {
      style: NwkStyle::Plain,
      significant_digits: Some(3),
      decimal_digits: None,
    };
    let s = newick_to_string(&g, &o).unwrap();
    assert_eq!("(A:0.123,B:0.2);", s);
  }

  // Empty branches
  #[test]
  fn test_write_empty_branches() {
    let g = newick_from_string("(,);").unwrap();
    let s = newick_to_string(&g, &opts(NwkStyle::Plain)).unwrap();
    assert_eq!("(,);", s);
  }

  // No branch length
  #[test]
  fn test_write_no_branch_length() {
    let g = newick_from_string("(A,B);").unwrap();
    let s = newick_to_string(&g, &opts(NwkStyle::Plain)).unwrap();
    assert_eq!("(A,B);", s);
  }

  // eNewick writer
  #[test]
  fn test_write_enewick_hybrid() {
    let g = newick_from_string("(A,B,((C,(Y)x#H1)c,(x#H1,D)d)e)f;").unwrap();
    let s = newick_to_string(&g, &opts(NwkStyle::Plain)).unwrap();
    // Should contain x#H1 marker twice but subtree only once
    let count = s.matches("x#H1").count();
    assert_eq!(count, 2, "hybrid marker should appear twice: {s}");
  }

  // BEAST string that looks like boolean/number must round-trip via quoting
  #[test]
  fn test_write_beast_string_true_roundtrip() {
    let mut g = NewickGraph::new();
    let mut attrs = BTreeMap::new();
    attrs.insert("val".to_owned(), NewickValue::String("TRUE".to_owned()));
    let node = g.add_node(NewickNodeData {
      name: Some("A".to_owned()),
      confidence: None,
      node_attrs: attrs,
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    });
    g.root = node;

    let written = newick_to_string(&g, &opts(NwkStyle::Beast)).unwrap();
    let parsed = newick_from_string(&written).unwrap();
    let a = parsed.nodes.iter().find(|n| n.name.as_deref() == Some("A")).unwrap();
    assert_eq!(Some(&NewickValue::String("TRUE".to_owned())), a.node_attrs.get("val"));
  }

  #[test]
  fn test_write_beast_string_numeric_roundtrip() {
    let mut g = NewickGraph::new();
    let mut attrs = BTreeMap::new();
    attrs.insert("id".to_owned(), NewickValue::String("123".to_owned()));
    let node = g.add_node(NewickNodeData {
      name: Some("A".to_owned()),
      confidence: None,
      node_attrs: attrs,
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    });
    g.root = node;

    let written = newick_to_string(&g, &opts(NwkStyle::Beast)).unwrap();
    let parsed = newick_from_string(&written).unwrap();
    let a = parsed.nodes.iter().find(|n| n.name.as_deref() == Some("A")).unwrap();
    assert_eq!(Some(&NewickValue::String("123".to_owned())), a.node_attrs.get("id"));
  }

  // Embedded double quotes in BEAST strings are escaped by doubling
  #[test]
  fn test_write_beast_string_embedded_quote_roundtrip() {
    let mut g = NewickGraph::new();
    let mut attrs = BTreeMap::new();
    attrs.insert("note".to_owned(), NewickValue::String("say \"hello\"".to_owned()));
    let node = g.add_node(NewickNodeData {
      name: Some("A".to_owned()),
      confidence: None,
      node_attrs: attrs,
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    });
    g.root = node;

    let written = newick_to_string(&g, &opts(NwkStyle::Beast)).unwrap();
    let parsed = newick_from_string(&written).unwrap();
    let a = parsed.nodes.iter().find(|n| n.name.as_deref() == Some("A")).unwrap();
    assert_eq!(
      Some(&NewickValue::String("say \"hello\"".to_owned())),
      a.node_attrs.get("note")
    );
  }

  // NHX rejects values containing reserved characters
  #[test]
  fn test_write_nhx_rejects_colon_in_value() {
    let mut g = NewickGraph::new();
    let mut attrs = BTreeMap::new();
    attrs.insert("key".to_owned(), NewickValue::String("a:b".to_owned()));
    let node = g.add_node(NewickNodeData {
      name: Some("A".to_owned()),
      confidence: None,
      node_attrs: attrs,
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    });
    g.root = node;

    let mut buf = Vec::new();
    let result = crate::write::newick_to_writer(
      &mut buf,
      &g,
      &NewickWriteOptions {
        style: NwkStyle::Nhx,
        ..Default::default()
      },
    );
    assert!(result.is_err(), "NHX writer should reject values containing ':'");
  }

  // All three styles produce different output for same annotated graph
  #[test]
  fn test_write_all_styles_differ() {
    let g = newick_from_string("(A[&prob=0.9]:0.1,B:0.2);").unwrap();
    let plain = newick_to_string(&g, &opts(NwkStyle::Plain)).unwrap();
    let beast = newick_to_string(&g, &opts(NwkStyle::Beast)).unwrap();
    let nhx = newick_to_string(&g, &opts(NwkStyle::Nhx)).unwrap();
    assert_ne!(plain, beast);
    assert_ne!(beast, nhx);
    assert_ne!(plain, nhx);
  }

  #[test]
  fn test_write_beast_key_escaping() {
    let mut g = NewickGraph::new();
    let mut attrs = BTreeMap::new();
    attrs.insert("posterior prob".to_owned(), NewickValue::Number(0.95));
    let node = g.add_node(NewickNodeData {
      name: Some("A".to_owned()),
      confidence: None,
      node_attrs: attrs,
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    });
    g.root = node;
    let s = newick_to_string(&g, &opts(NwkStyle::Beast)).unwrap();
    let parsed = newick_from_string(&s).unwrap();
    assert_eq!(
      Some(&NewickValue::Number(0.95)),
      parsed.nodes[0].node_attrs.get("posterior prob")
    );
  }

  #[test]
  fn test_write_newick_to_string_nhx_error() {
    let mut g = NewickGraph::new();
    let mut attrs = BTreeMap::new();
    attrs.insert("k".to_owned(), NewickValue::String("a:b".to_owned()));
    let node = g.add_node(NewickNodeData {
      name: Some("A".to_owned()),
      confidence: None,
      node_attrs: attrs,
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    });
    g.root = node;
    let result = newick_to_string(&g, &opts(NwkStyle::Nhx));
    assert!(
      result.is_err(),
      "newick_to_string should return Err for NHX reserved chars"
    );
  }
}

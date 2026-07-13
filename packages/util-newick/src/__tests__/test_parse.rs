#[cfg(test)]
mod tests {
  use crate::parse::newick_from_string;
  use crate::types::NewickValue;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[test]
  fn test_parse_single_leaf() {
    let g = newick_from_string("A;").unwrap();
    assert_eq!(g.nodes.len(), 1);
    assert_eq!(g.nodes[0].name.as_deref(), Some("A"));
    assert!(g.nodes[0].children.is_empty());
  }

  #[test]
  fn test_parse_two_leaves() {
    let g = newick_from_string("(A,B);").unwrap();
    assert_eq!(g.nodes.len(), 3);
    let root = &g.nodes[g.root];
    assert_eq!(root.children.len(), 2);
  }

  #[test]
  fn test_parse_named_internal() {
    let g = newick_from_string("(A,B)root;").unwrap();
    assert_eq!(g.nodes[g.root].name.as_deref(), Some("root"));
  }

  #[test]
  fn test_parse_branch_lengths() {
    let g = newick_from_string("(A:0.1,B:0.2):0.0;").unwrap();
    assert_eq!(g.edges.len(), 2);
    let lengths: Vec<f64> = g.edges.iter().filter_map(|e| e.data.branch_length).collect();
    assert!(lengths.contains(&0.1));
    assert!(lengths.contains(&0.2));
  }

  #[test]
  fn test_parse_scientific_notation() {
    let g = newick_from_string("(A:1.5e-3,B:2E4);").unwrap();
    let lengths: Vec<f64> = g.edges.iter().filter_map(|e| e.data.branch_length).collect();
    assert!(lengths.contains(&1.5e-3));
    assert!(lengths.contains(&2e4));
  }

  #[test]
  fn test_parse_negative_branch_length() {
    let g = newick_from_string("(A:-0.01,B:0.2);").unwrap();
    assert!(
      g.edges
        .iter()
        .filter_map(|e| e.data.branch_length)
        .any(|x| (x - (-0.01)).abs() < 1e-15)
    );
  }

  #[test]
  fn test_parse_quoted_label() {
    let g = newick_from_string("'node with spaces';").unwrap();
    assert_eq!(g.nodes[0].name.as_deref(), Some("node with spaces"));
  }

  #[test]
  fn test_parse_quoted_label_escaped_quote() {
    let g = newick_from_string("'it''s a name';").unwrap();
    assert_eq!(g.nodes[0].name.as_deref(), Some("it's a name"));
  }

  #[test]
  fn test_parse_empty_branches() {
    let g = newick_from_string("(,);").unwrap();
    assert_eq!(g.nodes.len(), 3);
    let root = &g.nodes[g.root];
    assert_eq!(root.children.len(), 2);
    // Anonymous children
    for &ei in &root.children {
      let child_idx = g.edges[ei].child;
      assert_eq!(g.nodes[child_idx].name, None);
    }
  }

  #[test]
  fn test_parse_empty_branches_three() {
    let g = newick_from_string("(A,,B);").unwrap();
    assert_eq!(g.nodes.len(), 4); // root + A + anonymous + B
    let root = &g.nodes[g.root];
    assert_eq!(root.children.len(), 3);
  }

  #[test]
  fn test_parse_rooted_marker() {
    let g = newick_from_string("[&R](A,B);").unwrap();
    assert_eq!(g.rooted, Some(true));
  }

  #[test]
  fn test_parse_unrooted_marker() {
    let g = newick_from_string("[&U](A,B);").unwrap();
    assert_eq!(g.rooted, Some(false));
  }

  #[test]
  fn test_parse_rooted_case_insensitive() {
    let g = newick_from_string("[&r](A,B);").unwrap();
    assert_eq!(g.rooted, Some(true));
  }

  #[test]
  fn test_parse_no_rooting_marker() {
    let g = newick_from_string("(A,B);").unwrap();
    assert_eq!(g.rooted, None);
  }

  // BEAST dialect tests
  #[test]
  fn test_parse_beast_node_attrs() {
    let g = newick_from_string("(A[&prob=0.95,rate=1.2],B);").unwrap();
    let a_idx = g.nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
    let attrs = &g.nodes[a_idx].node_attrs;
    assert_eq!(attrs.get("prob"), Some(&NewickValue::Number(0.95)));
    assert_eq!(attrs.get("rate"), Some(&NewickValue::Number(1.2)));
  }

  #[test]
  fn test_parse_beast_boolean_values() {
    let g = newick_from_string("(A[&fixed=TRUE,active=FALSE],B);").unwrap();
    let a_idx = g.nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
    let attrs = &g.nodes[a_idx].node_attrs;
    assert_eq!(attrs.get("fixed"), Some(&NewickValue::Boolean(true)));
    assert_eq!(attrs.get("active"), Some(&NewickValue::Boolean(false)));
  }

  #[test]
  fn test_parse_beast_boolean_case_insensitive() {
    let g = newick_from_string("(A[&x=true,y=False],B);").unwrap();
    let a_idx = g.nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
    let attrs = &g.nodes[a_idx].node_attrs;
    assert_eq!(attrs.get("x"), Some(&NewickValue::Boolean(true)));
    assert_eq!(attrs.get("y"), Some(&NewickValue::Boolean(false)));
  }

  #[test]
  fn test_parse_beast_array_values() {
    let g = newick_from_string("(A[&hpd={1.0,2.0,3.0}],B);").unwrap();
    let a_idx = g.nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
    let hpd = &g.nodes[a_idx].node_attrs["hpd"];
    match hpd {
      NewickValue::Array(arr) => {
        assert_eq!(arr.len(), 3);
        assert_eq!(arr[0], NewickValue::Number(1.0));
        assert_eq!(arr[1], NewickValue::Number(2.0));
        assert_eq!(arr[2], NewickValue::Number(3.0));
      },
      _ => panic!("Expected Array"),
    }
  }

  #[test]
  fn test_parse_beast_string_value() {
    let g = newick_from_string("(A[&country=USA],B);").unwrap();
    let a_idx = g.nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
    assert_eq!(
      g.nodes[a_idx].node_attrs.get("country"),
      Some(&NewickValue::String("USA".to_owned()))
    );
  }

  #[test]
  fn test_parse_beast_quoted_string_value() {
    let g = newick_from_string("(A[&label=\"hello, world\"],B);").unwrap();
    let a_idx = g.nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
    assert_eq!(
      g.nodes[a_idx].node_attrs.get("label"),
      Some(&NewickValue::String("hello, world".to_owned()))
    );
  }

  #[test]
  fn test_parse_beast_bare_key_boolean() {
    let g = newick_from_string("(A[&flagged],B);").unwrap();
    let a_idx = g.nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
    assert_eq!(
      g.nodes[a_idx].node_attrs.get("flagged"),
      Some(&NewickValue::Boolean(true))
    );
  }

  // BEAST2 canonical branch attr placement
  #[test]
  fn test_parse_beast2_branch_attrs() {
    let g = newick_from_string("(A:[&rate=1.5]0.1,B:0.2);").unwrap();
    let a_edge = g
      .edges
      .iter()
      .find(|e| g.nodes[e.child].name.as_deref() == Some("A"))
      .unwrap();
    assert_eq!(a_edge.data.branch_attrs.get("rate"), Some(&NewickValue::Number(1.5)));
    assert_eq!(a_edge.data.branch_length, Some(0.1));
  }

  // MrBayes-style: branch attrs after length
  #[test]
  fn test_parse_mrbayes_branch_attrs() {
    let g = newick_from_string("(A:0.1[&prob=0.99],B:0.2);").unwrap();
    let a_edge = g
      .edges
      .iter()
      .find(|e| g.nodes[e.child].name.as_deref() == Some("A"))
      .unwrap();
    assert_eq!(a_edge.data.branch_attrs.get("prob"), Some(&NewickValue::Number(0.99)));
  }

  // NHX dialect tests
  #[test]
  fn test_parse_nhx_attrs() {
    let g = newick_from_string("(A[&&NHX:S=human:B=90],B);").unwrap();
    let a_idx = g.nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
    let attrs = &g.nodes[a_idx].node_attrs;
    assert_eq!(attrs.get("S"), Some(&NewickValue::String("human".to_owned())));
    assert_eq!(attrs.get("B"), Some(&NewickValue::String("90".to_owned())));
  }

  // Raw comments
  #[test]
  fn test_parse_raw_comment() {
    let g = newick_from_string("(A[some comment],B);").unwrap();
    let a_idx = g.nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
    assert_eq!(g.nodes[a_idx].raw_comments, vec!["[some comment]"]);
    assert!(g.nodes[a_idx].node_attrs.is_empty());
  }

  #[test]
  fn test_parse_nested_bracket_comment() {
    let g = newick_from_string("(A[outer[inner]],B);").unwrap();
    let a_idx = g.nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
    assert_eq!(g.nodes[a_idx].raw_comments, vec!["[outer[inner]]"]);
  }

  // Malformed annotation falls back to raw
  #[test]
  fn test_parse_malformed_beast_fallback() {
    let g = newick_from_string("(A[&broken=],B);").unwrap();
    let a_idx = g.nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
    // Even malformed, the BEAST parser attempts to parse it
    // [&broken=] has key "broken" with empty value -> String("")
    assert!(g.nodes[a_idx].node_attrs.contains_key("broken"));
  }

  // eNewick hybrid tests
  #[test]
  fn test_parse_enewick_hybrid() {
    let g = newick_from_string("(A,B,((C,(Y)x#H1)c,(x#H1,D)d)e)f;").unwrap();
    // x#H1 should appear as a single node with two parents
    let hybrid_nodes: Vec<_> = g.nodes.iter().filter(|n| n.hybrid.is_some()).collect();
    assert_eq!(hybrid_nodes.len(), 1);
    let h = hybrid_nodes[0].hybrid.as_ref().unwrap();
    assert_eq!(h.kind.as_deref(), Some("H"));
    assert_eq!(h.index, 1);
    // Node x should have name "x"
    assert_eq!(hybrid_nodes[0].name.as_deref(), Some("x"));
  }

  #[test]
  fn test_parse_enewick_acceptor() {
    let g = newick_from_string("((A)x##LGT1,(x#LGT1,B));").unwrap();
    let hybrid_nodes: Vec<_> = g.nodes.iter().filter(|n| n.hybrid.is_some()).collect();
    assert_eq!(1, hybrid_nodes.len());
    let h = hybrid_nodes[0].hybrid.as_ref().unwrap();
    assert_eq!(Some("LGT"), h.kind.as_deref());
    assert_eq!(1, h.index);

    // Acceptor is edge-specific: exactly one of the two parent edges should be acceptor
    let hybrid_idx = g.nodes.iter().position(|n| n.hybrid.is_some()).unwrap();
    let parent_edges: Vec<_> = g.edges.iter().filter(|e| e.child == hybrid_idx).collect();
    assert_eq!(2, parent_edges.len());
    let acceptor_count = parent_edges.iter().filter(|e| e.data.is_acceptor).count();
    assert_eq!(1, acceptor_count, "exactly one parent edge should be the acceptor");
  }

  #[test]
  fn test_parse_enewick_no_name() {
    let g = newick_from_string("((A)#H1,(#H1,B));").unwrap();
    let hybrid_nodes: Vec<_> = g.nodes.iter().filter(|n| n.hybrid.is_some()).collect();
    assert_eq!(hybrid_nodes.len(), 1);
    assert_eq!(hybrid_nodes[0].name, None);
  }

  #[test]
  fn test_parse_hash_no_match_is_plain_name() {
    let g = newick_from_string("node#;").unwrap();
    assert_eq!(g.nodes[0].name.as_deref(), Some("node#"));
    assert!(g.nodes[0].hybrid.is_none());
  }

  #[test]
  fn test_parse_spec_no_names() {
    let g = newick_from_string("(,,(,));").unwrap();
    assert_eq!(6, g.nodes.len());
    assert_eq!(3, g.nodes[g.root].children.len());
  }

  #[test]
  fn test_parse_spec_leaf_names() {
    let g = newick_from_string("(A,B,(C,D));").unwrap();
    assert_eq!(6, g.nodes.len());
    let leaf_names: Vec<_> = g
      .nodes
      .iter()
      .filter(|n| n.children.is_empty())
      .filter_map(|n| n.name.as_deref())
      .collect();
    assert!(leaf_names.contains(&"A"));
    assert!(leaf_names.contains(&"B"));
    assert!(leaf_names.contains(&"C"));
    assert!(leaf_names.contains(&"D"));
  }

  #[test]
  fn test_parse_spec_all_names() {
    let g = newick_from_string("(A,B,(C,D)E)F;").unwrap();
    assert_eq!(6, g.nodes.len());
    assert_eq!(Some("F"), g.nodes[g.root].name.as_deref());
    assert_eq!(3, g.nodes[g.root].children.len());
  }

  #[test]
  fn test_parse_spec_distances_and_leaf_names() {
    let g = newick_from_string("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);").unwrap();
    assert_eq!(6, g.nodes.len());
    let a_edge = g
      .edges
      .iter()
      .find(|e| g.nodes[e.child].name.as_deref() == Some("A"))
      .unwrap();
    assert_eq!(Some(0.1), a_edge.data.branch_length);
  }

  #[test]
  fn test_parse_spec_distances_and_all_names() {
    let g = newick_from_string("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    assert_eq!(6, g.nodes.len());
    assert_eq!(Some("F"), g.nodes[g.root].name.as_deref());
    let e_edge = g
      .edges
      .iter()
      .find(|e| g.nodes[e.child].name.as_deref() == Some("E"))
      .unwrap();
    assert_eq!(Some(0.5), e_edge.data.branch_length);
  }

  #[test]
  fn test_parse_missing_semicolon() {
    let err = format!("{}", newick_from_string("(A,B)").unwrap_err());
    assert!(err.contains("Failed to parse Newick string"), "unexpected error: {err}");
  }

  #[test]
  fn test_parse_unmatched_paren() {
    let err = format!("{}", newick_from_string("(A,B;").unwrap_err());
    assert!(err.contains("Failed to parse Newick string"), "unexpected error: {err}");
  }

  #[test]
  fn test_parse_empty_string() {
    let err = format!("{}", newick_from_string("").unwrap_err());
    assert!(err.contains("Failed to parse Newick string"), "unexpected error: {err}");
  }

  // PartialEq tests
  #[test]
  fn test_eq_order_insensitive() {
    let g1 = newick_from_string("(A,B,C);").unwrap();
    let g2 = newick_from_string("(C,A,B);").unwrap();
    assert_eq!(g1, g2);
  }

  #[test]
  fn test_eq_ordered_sensitive() {
    let g1 = newick_from_string("(A,B,C);").unwrap();
    let g2 = newick_from_string("(C,A,B);").unwrap();
    assert!(!g1.eq_ordered(&g2));
  }

  #[test]
  fn test_eq_ordered_same() {
    let g1 = newick_from_string("(A,B,C);").unwrap();
    let g2 = newick_from_string("(A,B,C);").unwrap();
    assert!(g1.eq_ordered(&g2));
  }

  #[test]
  fn test_eq_with_branch_lengths() {
    let g1 = newick_from_string("(A:0.1,B:0.2);").unwrap();
    let g2 = newick_from_string("(B:0.2,A:0.1);").unwrap();
    assert_eq!(g1, g2);
  }

  #[test]
  fn test_eq_unnamed_internal_reorder() {
    let g1 = newick_from_string("((A,B),(C,D));").unwrap();
    let g2 = newick_from_string("((C,D),(A,B));").unwrap();
    assert_eq!(g1, g2);
  }

  #[test]
  fn test_neq_different_topology() {
    let g1 = newick_from_string("((A,B),C);").unwrap();
    let g2 = newick_from_string("(A,(B,C));").unwrap();
    assert_ne!(g1, g2);
  }

  #[test]
  fn test_eq_rooting_matters() {
    let g1 = newick_from_string("[&R](A,B);").unwrap();
    let g2 = newick_from_string("[&U](A,B);").unwrap();
    assert_ne!(g1, g2);
  }

  // Whitespace handling
  #[test]
  fn test_parse_whitespace() {
    let g = newick_from_string("  ( A : 0.1 , B : 0.2 ) ; ").unwrap();
    assert_eq!(g.nodes.len(), 3);
  }

  #[test]
  fn test_parse_newlines() {
    let g = newick_from_string(indoc! {"
      (
        A:0.1,
        B:0.2
      );
    "})
    .unwrap();
    assert_eq!(g.nodes.len(), 3);
  }

  // Single-child internal node
  #[test]
  fn test_parse_single_child_internal() {
    let g = newick_from_string("((A));").unwrap();
    assert_eq!(g.nodes.len(), 3);
  }

  // Deep nesting
  #[test]
  fn test_parse_deep_nesting() {
    let g = newick_from_string("((((A,B),C),D),E);").unwrap();
    assert_eq!(g.nodes.len(), 9);
  }

  // Comment on edge with no length
  #[test]
  fn test_parse_comment_no_length() {
    let g = newick_from_string("(A[&note=yes],B);").unwrap();
    let a_idx = g.nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
    assert_eq!(
      g.nodes[a_idx].node_attrs.get("note"),
      Some(&NewickValue::String("yes".to_owned()))
    );
  }

  // Confidence parsing: Biopython heuristic on internal node labels

  #[rustfmt::skip]
  #[rstest]
  #[case::float(        "(A:0.1,B:0.2)0.999:0.3;", (None,          Some(0.999)))]
  #[case::integer(      "(A,B)100;",                (None,          Some(100.0)))]
  #[case::zero(         "(A,B)0;",                  (None,          Some(0.0)  ))]
  #[case::scientific(   "(A,B)1e-5;",               (None,          Some(1e-5) ))]
  #[case::text_label(   "(A,B)root;",               (Some("root"),  None       ))]
  #[case::no_label(     "(A,B);",                   (None,          None       ))]
  #[trace]
  fn test_parse_confidence_on_root(
    #[case] nwk: &str,
    #[case] (expected_name, expected_confidence): (Option<&str>, Option<f64>),
  ) {
    let g = newick_from_string(nwk).unwrap();
    let root = &g.nodes[g.root];
    assert_eq!(root.name.as_deref(), expected_name);
    assert_eq!(root.confidence, expected_confidence);
  }

  #[test]
  fn test_parse_confidence_not_applied_to_leaves() {
    let g = newick_from_string("(0.999:0.1,B:0.2);").unwrap();
    let leaf = g.nodes.iter().find(|n| n.name.as_deref() == Some("0.999")).unwrap();
    assert_eq!(leaf.confidence, None);
  }

  #[test]
  fn test_parse_confidence_nested_internals() {
    let g = newick_from_string("((A,B)0.95:0.1,(C,D)0.80:0.2)0.50;").unwrap();
    let internal_count = g
      .nodes
      .iter()
      .filter(|n| !n.children.is_empty() && n.confidence.is_some())
      .count();
    assert_eq!(internal_count, 3);
    assert_eq!(g.nodes[g.root].confidence, Some(0.50));
  }

  #[test]
  fn test_parse_confidence_mixed_named_and_float_internals() {
    let g = newick_from_string("((A,B)clade1:0.1,(C,D)0.90:0.2)root;").unwrap();
    let clade1 = g.nodes.iter().find(|n| n.name.as_deref() == Some("clade1")).unwrap();
    assert_eq!(clade1.confidence, None);

    let float_node = g.nodes.iter().find(|n| n.confidence == Some(0.90)).unwrap();
    assert_eq!(float_node.name, None);

    let root = &g.nodes[g.root];
    assert_eq!(root.name.as_deref(), Some("root"));
    assert_eq!(root.confidence, None);
  }

  #[test]
  fn test_parse_confidence_roundtrip() {
    let g = newick_from_string("(A:0.1,B:0.2)0.999;").unwrap();
    let written = crate::write::newick_to_string(&g, &crate::types::NewickWriteOptions::default()).unwrap();
    assert_eq!(written, "(A:0.1,B:0.2)0.999;");
  }

  #[test]
  fn test_parse_beast_infinity_as_string() {
    let g = newick_from_string("(A[&score=-infinity],B);").unwrap();
    let a_idx = g.nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
    assert!(
      matches!(g.nodes[a_idx].node_attrs.get("score"), Some(NewickValue::String(_))),
      "Non-finite value should parse as String, not Number"
    );
  }
}

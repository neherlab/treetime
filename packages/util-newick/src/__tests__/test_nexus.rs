#[cfg(test)]
mod tests {
  use crate::nexus::{nexus_from_string, nexus_to_string};
  use crate::types::{NewickWriteOptions, NwkStyle};
  use indoc::indoc;
  use pretty_assertions::assert_eq;

  fn opts(style: NwkStyle) -> NewickWriteOptions {
    NewickWriteOptions {
      style,
      significant_digits: None,
      decimal_digits: None,
    }
  }

  #[test]
  fn test_nexus_simple() {
    let input = indoc! {"
      #NEXUS

      Begin Trees;
        Tree tree1 = (A:0.1,B:0.2);
      End;
    "};

    let trees = nexus_from_string(input).unwrap();
    assert_eq!(1, trees.len());
    assert_eq!("tree1", trees[0].name);
    assert_eq!(3, trees[0].graph.nodes.len());
  }

  #[test]
  fn test_nexus_multiple_trees() {
    let input = indoc! {"
      #NEXUS

      Begin Trees;
        Tree t1 = (A,B);
        Tree t2 = (C,D);
      End;
    "};

    let trees = nexus_from_string(input).unwrap();
    assert_eq!(2, trees.len());
    assert_eq!("t1", trees[0].name);
    assert_eq!("t2", trees[1].name);
  }

  #[test]
  fn test_nexus_translate() {
    let input = indoc! {"
      #NEXUS

      Begin Trees;
        Translate
          1 Homo_sapiens,
          2 Pan_troglodytes,
          3 Gorilla_gorilla
        ;
        Tree tree1 = (1:0.1,2:0.2,3:0.3);
      End;
    "};

    let trees = nexus_from_string(input).unwrap();
    assert_eq!(1, trees.len());
    let names: Vec<_> = trees[0].graph.nodes.iter().filter_map(|n| n.name.as_deref()).collect();
    assert!(names.contains(&"Homo_sapiens"));
    assert!(names.contains(&"Pan_troglodytes"));
    assert!(names.contains(&"Gorilla_gorilla"));
    assert!(!names.contains(&"1"));
  }

  #[test]
  fn test_nexus_rooted_marker() {
    let input = indoc! {"
      #NEXUS
      Begin Trees;
        Tree tree1 = [&R] (A,B);
      End;
    "};

    let trees = nexus_from_string(input).unwrap();
    assert_eq!(Some(true), trees[0].graph.rooted);
  }

  #[test]
  fn test_nexus_with_annotations() {
    let input = indoc! {"
      #NEXUS
      Begin Trees;
        Tree tree1 = (A[&prob=0.9]:0.1,B:0.2);
      End;
    "};

    let trees = nexus_from_string(input).unwrap();
    let a_idx = trees[0]
      .graph
      .nodes
      .iter()
      .position(|n| n.name.as_deref() == Some("A"))
      .unwrap();
    assert!(trees[0].graph.nodes[a_idx].node_attrs.contains_key("prob"));
  }

  #[test]
  fn test_nexus_unknown_blocks_ignored() {
    let input = indoc! {"
      #NEXUS

      Begin Data;
        something irrelevant;
      End;

      Begin Trees;
        Tree tree1 = (A,B);
      End;

      Begin figtree;
        display settings;
      End;
    "};

    let trees = nexus_from_string(input).unwrap();
    assert_eq!(1, trees.len());
    assert_eq!("tree1", trees[0].name);
  }

  #[test]
  fn test_nexus_case_insensitive_header() {
    let input = indoc! {"
      #nexus
      begin trees;
        tree t1 = (A,B);
      end;
    "};

    let trees = nexus_from_string(input).unwrap();
    assert_eq!(1, trees.len());
  }

  #[test]
  fn test_nexus_missing_header_error() {
    let input = indoc! {"
      Begin Trees;
        Tree t1 = (A,B);
      End;
    "};
    let err = format!("{}", nexus_from_string(input).unwrap_err());
    assert!(err.contains("#NEXUS"), "unexpected error: {err}");
  }

  #[test]
  fn test_nexus_roundtrip() {
    let input = indoc! {"
      #NEXUS
      Begin Trees;
        Tree tree1 = (A:0.1,B:0.2);
      End;
    "};

    let trees = nexus_from_string(input).unwrap();
    let output = nexus_to_string(&trees, &opts(NwkStyle::Plain)).unwrap();
    let trees2 = nexus_from_string(&output).unwrap();

    assert_eq!(trees.len(), trees2.len());
    assert_eq!(trees[0].name, trees2[0].name);
    assert_eq!(trees[0].graph, trees2[0].graph);
  }

  #[test]
  fn test_nexus_translate_quoted_names() {
    let input = indoc! {"
      #NEXUS
      Begin Trees;
        Translate
          1 'Homo sapiens',
          2 'Pan troglodytes'
        ;
        Tree tree1 = (1:0.1,2:0.2);
      End;
    "};

    let trees = nexus_from_string(input).unwrap();
    let names: Vec<_> = trees[0].graph.nodes.iter().filter_map(|n| n.name.as_deref()).collect();
    assert!(names.contains(&"Homo sapiens"));
    assert!(names.contains(&"Pan troglodytes"));
  }

  #[test]
  fn test_nexus_translate_comma_in_quoted_name() {
    let input = indoc! {"
      #NEXUS
      Begin Trees;
        Translate
          1 'A,B',
          2 C
        ;
        Tree tree1 = (1:0.1,2:0.2);
      End;
    "};

    let trees = nexus_from_string(input).unwrap();
    let names: Vec<_> = trees[0].graph.nodes.iter().filter_map(|n| n.name.as_deref()).collect();
    assert!(names.contains(&"A,B"), "Expected 'A,B' in {names:?}");
    assert!(names.contains(&"C"));
  }

  #[test]
  fn test_nexus_writer_taxa_block() {
    let input = indoc! {"
      #NEXUS
      Begin Trees;
        Tree tree1 = (A:0.1,B:0.2);
      End;
    "};

    let trees = nexus_from_string(input).unwrap();
    let output = nexus_to_string(&trees, &opts(NwkStyle::Plain)).unwrap();

    assert!(output.contains("Begin Taxa;"));
    assert!(output.contains("ntax=2"));
    assert!(output.contains('A'));
    assert!(output.contains('B'));
  }

  #[test]
  fn test_nexus_mixed_case_header() {
    let input = indoc! {"
      #Nexus
      Begin Trees;
        Tree t1 = (A,B);
      End;
    "};
    let trees = nexus_from_string(input).unwrap();
    assert_eq!(1, trees.len());
  }

  #[test]
  fn test_nexus_semicolon_inside_comment() {
    let input = indoc! {"
      #NEXUS
      Begin Trees;
        Tree t1 = (A[comment;still]:0.1,B:0.2);
      End;
    "};
    let trees = nexus_from_string(input).unwrap();
    assert_eq!(1, trees.len());
    assert_eq!(3, trees[0].graph.nodes.len());
  }
}

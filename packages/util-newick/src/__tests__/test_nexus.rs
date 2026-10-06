#[cfg(test)]
mod tests {
  use crate::nexus::{nexus_from_string, nexus_to_string};
  use crate::types::{NewickReadOptions, NewickWriteOptions, NwkStyle};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

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

    let trees = nexus_from_string(input, &NewickReadOptions::default()).unwrap();
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

    let trees = nexus_from_string(input, &NewickReadOptions::default()).unwrap();
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

    let trees = nexus_from_string(input, &NewickReadOptions::default()).unwrap();
    assert_eq!(1, trees.len());
    let names: Vec<_> = trees[0].graph.nodes.iter().filter_map(|n| n.name()).collect();
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

    let trees = nexus_from_string(input, &NewickReadOptions::default()).unwrap();
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

    let trees = nexus_from_string(input, &NewickReadOptions::default()).unwrap();
    let a_idx = trees[0].graph.nodes.iter().position(|n| n.name() == Some("A")).unwrap();
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

    let trees = nexus_from_string(input, &NewickReadOptions::default()).unwrap();
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

    let trees = nexus_from_string(input, &NewickReadOptions::default()).unwrap();
    assert_eq!(1, trees.len());
  }

  #[test]
  fn test_nexus_missing_header_error() {
    let input = indoc! {"
      Begin Trees;
        Tree t1 = (A,B);
      End;
    "};
    let err = format!(
      "{}",
      nexus_from_string(input, &NewickReadOptions::default()).unwrap_err()
    );
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

    let trees = nexus_from_string(input, &NewickReadOptions::default()).unwrap();
    let output = nexus_to_string(&trees, &opts(NwkStyle::Plain)).unwrap();
    let trees2 = nexus_from_string(&output, &NewickReadOptions::default()).unwrap();

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

    let trees = nexus_from_string(input, &NewickReadOptions::default()).unwrap();
    let names: Vec<_> = trees[0].graph.nodes.iter().filter_map(|n| n.name()).collect();
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

    let trees = nexus_from_string(input, &NewickReadOptions::default()).unwrap();
    let names: Vec<_> = trees[0].graph.nodes.iter().filter_map(|n| n.name()).collect();
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

    let trees = nexus_from_string(input, &NewickReadOptions::default()).unwrap();
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
    let trees = nexus_from_string(input, &NewickReadOptions::default()).unwrap();
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
    let trees = nexus_from_string(input, &NewickReadOptions::default()).unwrap();
    assert_eq!(1, trees.len());
    assert_eq!(3, trees[0].graph.nodes.len());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::translate_name_ending_in_end(  "Begin Trees;\n Translate 1 A, 2 Ostend;\n Tree t = (1,2);\nEnd;",                 vec![("t", vec!["A", "Ostend"])])]
  #[case::comment_with_begin(            "[ begin here; ]\nBegin Trees;\n Tree t = (A,B);\nEnd;",                         vec![("t", vec!["A", "B"])])]
  #[case::tree_name_with_translate(      "Begin Trees;\n Tree untranslated = (d:1,B:1);\nEnd;",                            vec![("untranslated", vec!["d", "B"])])]
  #[case::beast1_tree_comment(           "Begin trees;\n tree STATE_0 [&lnP=-3195.24] = [&R] (A:[&rate=1.0]1.0,B:1.0);\nEnd;", vec![("STATE_0", vec!["A", "B"])])]
  #[case::translate_per_block(           "Begin Trees;\n Translate 1 A, 2 B;\n Tree t1 = (1,2);\nEnd;\nBegin Trees;\n Translate 1 X, 2 Y;\n Tree t2 = (1,2);\nEnd;", vec![("t1", vec!["A", "B"]), ("t2", vec!["X", "Y"])])]
  #[case::endblock(                      "Begin Trees;\n Tree t = (A,B);\nEndBlock;",                                     vec![("t", vec!["A", "B"])])]
  #[case::tab_after_begin(               "Begin\tTrees;\n Tree t = (A,B);\nEnd;",                                         vec![("t", vec!["A", "B"])])]
  #[case::mrbayes_note(                  "Begin Trees;\n [Note: This tree contains information]\n Tree con_50 = (A,B);\nEnd;", vec![("con_50", vec!["A", "B"])])]
  #[case::paup_star(                     "Begin Trees;\n Tree * PAUP_1 = [&U] (A,B);\nEnd;",                              vec![("PAUP_1", vec!["A", "B"])])]
  #[case::utree(                         "Begin Trees;\n UTree t = (A,B);\nEnd;",                                          vec![("t", vec!["A", "B"])])]
  #[case::quoted_name_with_equals(       "Begin Trees;\n Tree 'a=b' = (A,B);\nEnd;",                                       vec![("a=b", vec!["A", "B"])])]
  #[case::data_block_skipped(            "Begin Data;\n Matrix\n A ACGT\n B 'AC;GT'\n ;\nEnd;\nBegin Trees;\n Tree t = (A,B);\nEnd;", vec![("t", vec!["A", "B"])])]
  #[trace]
  fn test_nexus_reads_trees(#[case] body: &str, #[case] expected: Vec<(&str, Vec<&str>)>) {
    let trees = nexus_from_string(&format!("#NEXUS\n{body}\n"), &NewickReadOptions::default()).unwrap();

    let actual: Vec<(&str, Vec<&str>)> = trees
      .iter()
      .map(|tree| (tree.name.as_str(), tree.graph.nodes.iter().filter_map(|node| node.name()).collect()))
      .collect();
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_end(       "Begin Trees;\n Tree t = (A,B);\n")]
  #[case::no_semicolon( "Begin Trees;\n Tree t = (A,B)\nEnd")]
  #[trace]
  fn test_nexus_unterminated_block_is_error(#[case] body: &str) {
    let actual = nexus_from_string(&format!("#NEXUS\n{body}"), &NewickReadOptions::default());

    assert!(format!("{:#}", actual.unwrap_err()).starts_with("Failed to parse Nexus input: "));
  }

  #[test]
  fn test_nexus_roundtrip_quotes_names_with_punctuation() {
    let graph = crate::parse::newick_from_string("('hCoV-19/USA/1',B);", &NewickReadOptions::default()).unwrap();
    let trees = vec![crate::types::NexusTree {
      name: "a=b".to_owned(),
      graph,
    }];

    let written = nexus_to_string(&trees, &opts(NwkStyle::Plain)).unwrap();
    let parsed = nexus_from_string(&written, &NewickReadOptions::default()).unwrap();

    let expected = indoc! {"
      #NEXUS

      Begin Taxa;
        Dimensions ntax=2;
        TaxLabels
          B
          'hCoV-19/USA/1'
        ;
      End;

      Begin Trees;
        Tree 'a=b' = (hCoV-19/USA/1,B);
      End;
    "};
    assert_eq!(
      (expected, "a=b", true),
      (
        written.as_str(),
        parsed[0].name.as_str(),
        parsed[0].graph == trees[0].graph
      )
    );
  }
}

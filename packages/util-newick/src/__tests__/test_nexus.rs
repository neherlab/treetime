#[cfg(test)]
mod tests {
  use crate::__tests__::test_read_basic::tests::helpers::summary;
  use crate::dialect::NewickDialect;
  use crate::model::comment::{MrBayesComment, MrBayesKind, NewickComment};
  use crate::model::value::NewickValue;
  use crate::nexus::read::{is_nexus, nexus_from_str, nexus_trees};
  use crate::nexus::types::{NexusCommand, NexusTreeRef, NexusWriteOptions};
  use crate::nexus::write::nexus_to_string;
  use crate::read::options::NewickReadOptions;
  use crate::write::options::NewickWriteOptions;
  use helpers::{read_error, tolerant, trees};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::simple(                        "Begin Trees;\n Tree tree1 = (A:0.1,B:0.2);\nEnd;",                        vec![("tree1", vec!["- ", "A :0.1", "B :0.2"])])]
  #[case::two_trees(                     "Begin Trees;\n Tree t1 = (A,B);\n Tree t2 = (C,D);\nEnd;",                 vec![("t1", vec!["- ", "A ", "B "]), ("t2", vec!["- ", "C ", "D "])])]
  #[case::translate(                     "Begin Trees;\n Translate\n  1 Homo_sapiens,\n  2 'Pan troglodytes'\n ;\n Tree t = (1:0.1,2:0.2);\nEnd;", vec![("t", vec!["- ", "Homo_sapiens :0.1", "Pan troglodytes :0.2"])])]
  #[case::translate_comma_in_name(       "Begin Trees;\n Translate 1 'A,B', 2 C;\n Tree t = (1,2);\nEnd;",               vec![("t", vec!["- ", "A,B ", "C "])])]
  #[case::translate_internal_unchanged(  "Begin Trees;\n Translate 1 A, 2 B;\n Tree t = ((1,2)1);\nEnd;",              vec![("t", vec!["- ", "- support=1(label)", "A ", "B "])])]
  #[case::translate_name_ending_in_end(  "Begin Trees;\n Translate 1 A, 2 Ostend;\n Tree t = (1,2);\nEnd;",              vec![("t", vec!["- ", "A ", "Ostend "])])]
  #[case::comment_with_begin(            "[ begin here; ]\nBegin Trees;\n Tree t = (A,B);\nEnd;",                    vec![("t", vec!["- ", "A ", "B "])])]
  #[case::semicolon_in_comment(          "Begin Trees;\n Tree t = (A[x;y]:0.1,B);\nEnd;",                           vec![("t", vec!["- ", "A :0.1", "B "])])]
  #[case::translate_per_block(           "Begin Trees;\n Translate 1 A, 2 B;\n Tree t1 = (1,2);\nEnd;\nBegin Trees;\n Translate 1 X, 2 Y;\n Tree t2 = (1,2);\nEnd;", vec![("t1", vec!["- ", "A ", "B "]), ("t2", vec!["- ", "X ", "Y "])])]
  #[case::endblock(                      "Begin Trees;\n Tree t = (A,B);\nEndBlock;",                                vec![("t", vec!["- ", "A ", "B "])])]
  #[case::tab_after_begin(               "Begin\tTrees;\n Tree t = (A,B);\nEnd;",                                    vec![("t", vec!["- ", "A ", "B "])])]
  #[case::mrbayes_note(                  "Begin Trees;\n [Note: This tree contains information]\n Tree con_50 = (A,B);\nEnd;", vec![("con_50", vec!["- ", "A ", "B "])])]
  #[case::paup_star(                     "Begin Trees;\n Tree * PAUP_1 = [&U] (A,B);\nEnd;",                         vec![("PAUP_1", vec!["- ", "A ", "B "])])]
  #[case::utree(                         "Begin Trees;\n UTree t = (A,B);\nEnd;",                                    vec![("t", vec!["- ", "A ", "B "])])]
  #[case::quoted_name_with_equals(       "Begin Trees;\n Tree 'a=b' = (A,B);\nEnd;",                                 vec![("a=b", vec!["- ", "A ", "B "])])]
  #[case::data_block_skipped(            "Begin Data;\n Matrix\n A ACGT\n B 'AC;GT'\n ;\nEnd;\nBegin Trees;\n Tree t = (A,B);\nEnd;", vec![("t", vec!["- ", "A ", "B "])])]
  #[case::lowercase(                     "begin trees;\n tree t = (A,B);\nend;",                                     vec![("t", vec!["- ", "A ", "B "])])]
  #[trace]
  fn test_nexus_reads_trees(#[case] body: &str, #[case] expected: Vec<(&str, Vec<&str>)>) {
    let expected: Vec<(String, Vec<String>)> = expected
      .into_iter()
      .map(|(name, nodes)| (name.to_owned(), nodes.into_iter().map(str::to_owned).collect()))
      .collect();

    assert_eq!(expected, trees(&format!("#NEXUS\n{body}\n"), &NewickReadOptions::default()));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_header(           "Begin Trees;\n Tree t = (A,B);\nEnd;",                        "line 1, column 1: The input does not start with a '#NEXUS' header")]
  #[case::no_end(              "#NEXUS\nBegin Trees;\n Tree t = (A,B);\n",                    "line 4, column 1: The Trees block never ends")]
  #[case::no_semicolon(        "#NEXUS\nBegin Trees;\n Tree t = (A,B)\nEnd",                  "line 4, column 4: The command does not end with ';'")]
  #[case::tree_name_spaces(    "#NEXUS\nBegin Trees;\n tree my tree = (A,B);\nEnd;",           "line 3, column 10: expected comment, '='")]
  #[case::translate_spaces(    "#NEXUS\nBegin Trees;\n Translate 1 Homo sapiens, 2 B;\nEnd;", "line 3, column 19: expected ',', ';'")]
  #[case::unknown_taxon(       "#NEXUS\nBegin Taxa;\n TaxLabels A B;\nEnd;\nBegin Trees;\n Tree t = (A,C);\nEnd;", r#"line 6, column 2: The tree label "C" is neither a taxon label nor a 'Translate' key"#)]
  #[case::taxa_count(          "#NEXUS\nBegin Taxa;\n Dimensions NTax=3;\n TaxLabels A B;\nEnd;", "line 4, column 2: The 'Dimensions' command declares 3 taxa, but 'TaxLabels' lists 2")]
  #[case::outside_block(       "#NEXUS\nTree t = (A,B);",                                       "line 2, column 1: A command appears outside of a block")]
  #[case::nested_begin(        "#NEXUS\nBegin Trees;\nBegin Taxa;\nEnd;",                      "line 3, column 1: A block begins before the Trees block ends")]
  #[case::duplicate_key(       "#NEXUS\nBegin Trees;\n Translate 1 A, 1 B;\nEnd;",              r#"line 3, column 2: The 'Translate' command lists the key "1" more than once"#)]
  #[case::bad_tree(            "#NEXUS\nBegin Trees;\n Tree t = (A,B));\nEnd;",                 "line 3, column 16: unexpected ')' without a matching '('")]
  #[trace]
  fn test_nexus_strict_errors(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(expected, read_error(input, &NewickReadOptions::default()));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::tree_name_spaces(  "#NEXUS\nBegin Trees;\n tree my tree = (A,B);\nEnd;",            vec![("my tree", vec!["- ", "A ", "B "])], 0)]
  #[case::translate_spaces(  "#NEXUS\nBegin Trees;\n Translate 1 Homo sapiens, 2 B;\n Tree t = (1,2);\nEnd;", vec![("t", vec!["- ", "Homo sapiens ", "B "])], 0)]
  #[case::unknown_taxon(     "#NEXUS\nBegin Taxa;\n TaxLabels A B;\nEnd;\nBegin Trees;\n Tree t = (A,C);\nEnd;", vec![("t", vec!["- ", "A ", "C "])], 1)]
  #[case::no_end(            "#NEXUS\nBegin Trees;\n Tree t = (A,B);\n",                     vec![("t", vec!["- ", "A ", "B "])], 1)]
  #[trace]
  fn test_nexus_tolerant(#[case] input: &str, #[case] expected: Vec<(&str, Vec<&str>)>, #[case] warning_count: usize) {
    let file = nexus_from_str(input, &tolerant()).unwrap();

    let actual: Vec<(String, Vec<String>)> = file.trees.iter().map(|tree| (tree.name.clone(), summary(&tree.tree.graph))).collect();
    let warnings = file.warnings.len() + file.trees.iter().map(|tree| tree.tree.warnings.len()).sum::<usize>();
    let expected: Vec<(String, Vec<String>)> = expected
      .into_iter()
      .map(|(name, nodes)| (name.to_owned(), nodes.into_iter().map(str::to_owned).collect()))
      .collect();
    assert_eq!((expected, warning_count), (actual, warnings));
  }

  #[test]
  fn test_nexus_taxon_numbers_refer_to_taxa() {
    let input = "#NEXUS\nBegin Taxa;\n TaxLabels A B;\nEnd;\nBegin Trees;\n Tree t = (2,1);\nEnd;\n";

    assert_eq!(
      vec![("t".to_owned(), vec!["- ".to_owned(), "B ".to_owned(), "A ".to_owned()])],
      trees(input, &NewickReadOptions::default())
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::properties_rooted(  "Begin Trees;\n PROPERTIES rooted=yes;\n Tree t = (A,B);\nEnd;",          Some(true))]
  #[case::properties_reset(   "Begin Trees;\n PROPERTIES rooted=yes;\nEnd;\nBegin Trees;\n Tree t = (A,B);\nEnd;", None)]
  #[case::tree_comment_wins(  "Begin Trees;\n PROPERTIES rooted=no;\n Tree t = [&R] (A,B);\nEnd;",     Some(true))]
  #[trace]
  fn test_nexus_rooting(#[case] body: &str, #[case] expected: Option<bool>) {
    let options = NewickReadOptions {
      dialects: vec![NewickDialect::Beast],
      ..NewickReadOptions::default()
    };

    let file = nexus_from_str(&format!("#NEXUS\n{body}\n"), &options).unwrap();

    assert_eq!(expected, file.trees[0].tree.graph.rooted());
  }

  #[test]
  fn test_nexus_tree_command_comments_use_tree_dialect() {
    let input = "#NEXUS\nBegin trees;\n tree STATE_0 [&lnP=-3195.24] = [&R] (A:[&rate=1.0]1.0,B:1.0);\nEnd;\n";

    let file = nexus_from_str(input, &NewickReadOptions::all_dialects()).unwrap();

    let tree = &file.trees[0];
    assert_eq!(
      (
        NewickDialect::Beast,
        vec![NewickComment::Beast(vec![(
          "lnP".to_owned(),
          NewickValue::Number(-3195.24)
        )])]
      ),
      (tree.tree.dialect, tree.comments.clone())
    );
  }

  #[test]
  fn test_nexus_mrbayes_sample_file() {
    let input = indoc! {"
      #NEXUS
      [ID: 1234]
      begin trees;
         translate
            1 Homo_sapiens,
            2 Pan,
            3 Gorilla;
         tree gen.0 = [&U] ((1:0.1[&B TK02Brlens 0.1],2:0.1):0.2,3:0.3);
      end;
    "};
    let options = NewickReadOptions::all_dialects();

    let file = nexus_from_str(input, &options).unwrap();

    let tree = &file.trees[0];
    let leaf_edge = tree
      .tree
      .graph
      .child_edges(tree.tree.graph.children(tree.tree.graph.root()).next().unwrap())[0];
    let expected_comment = NewickComment::MrBayesMcmc(MrBayesComment {
      kind: MrBayesKind::B,
      name: "TK02Brlens".to_owned(),
      values: vec!["0.1".to_owned()],
    });
    assert_eq!(
      (
        "gen.0",
        NewickDialect::MrBayes,
        vec!["- ", "- :0.2", "Homo_sapiens :0.1", "Pan :0.1", "Gorilla :0.3"]
          .into_iter()
          .map(str::to_owned)
          .collect::<Vec<_>>(),
        vec![expected_comment]
      ),
      (
        tree.name.as_str(),
        tree.tree.dialect,
        summary(&tree.tree.graph),
        tree
          .tree
          .graph
          .edge(leaf_edge)
          .data()
          .comments()
          .iter()
          .map(|comment| comment.comment.clone())
          .collect::<Vec<_>>()
      )
    );
  }

  #[test]
  fn test_nexus_reports_skipped_commands() {
    let input = "#NEXUS\nBegin Data;\n Format datatype=dna;\nEnd;\nBegin Trees;\n Title my_trees;\n Tree t = (A,B);\nEnd;\nBegin figtree;\n set appearance.foregroundColour=#-16777216;\nEnd;\n";

    let file = nexus_from_str(input, &NewickReadOptions::default()).unwrap();

    let expected = vec![
      NexusCommand {
        block: Some("Data".to_owned()),
        command: "Format".to_owned(),
        line: 3,
      },
      NexusCommand {
        block: Some("Trees".to_owned()),
        command: "Title".to_owned(),
        line: 6,
      },
      NexusCommand {
        block: Some("figtree".to_owned()),
        command: "set".to_owned(),
        line: 10,
      },
    ];
    assert_eq!(expected, file.skipped);
  }

  #[test]
  fn test_nexus_trees_one_by_one() {
    let input = "#NEXUS\nBegin Trees;\n Tree t1 = (A,B);\n Tree t2 = (C,D);\nEnd;\n";

    let mut reading = nexus_trees(input.as_bytes(), NewickReadOptions::default());
    let names: Vec<String> = reading.by_ref().map(|tree| tree.unwrap().name).collect();

    assert_eq!(
      (vec!["t1".to_owned(), "t2".to_owned()], 0),
      (names, reading.skipped().len())
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::plain(       "\u{feff}  #nexus\n",  true)]
  #[case::with_suffix( "#NEXUSX\n",           false)]
  #[case::newick(      "(A,B);",              false)]
  #[trace]
  fn test_nexus_header_detection(#[case] input: &str, #[case] expected: bool) {
    assert_eq!(expected, is_nexus(input));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::without_translate(false, indoc! {"
    #NEXUS

    Begin Taxa;
      Dimensions NTax=3;
      TaxLabels
        'hCoV-19/USA/1'
        B
        'a b'
      ;
    End;

    Begin Trees;
      Tree 'tree one' [&lnP=-1] = [&R](hCoV-19/USA/1:0.1,(B,'a b')95:0.2);
      Tree t2 = [&U](B,'a b');
    End;
  "})]
  #[case::with_translate(true, indoc! {"
    #NEXUS

    Begin Taxa;
      Dimensions NTax=3;
      TaxLabels
        'hCoV-19/USA/1'
        B
        'a b'
      ;
    End;

    Begin Trees;
      Translate
        1 'hCoV-19/USA/1',
        2 B,
        3 'a b'
      ;
      Tree 'tree one' [&lnP=-1] = [&R](1:0.1,(2,3)95:0.2);
      Tree t2 = [&U](2,3);
    End;
  "})]
  #[trace]
  fn test_nexus_write_and_read_back(#[case] translate: bool, #[case] expected: &str) {
    let source = indoc! {"
      #NEXUS
      Begin Trees;
        Tree 'tree one' [&lnP=-1] = [&R] ('hCoV-19/USA/1':0.1,(B,'a b')95:0.2);
        Tree t2 = [&U] (B,'a b');
      End;
    "};
    let read_options = NewickReadOptions {
      dialects: vec![NewickDialect::Beast],
      ..NewickReadOptions::default()
    };
    let file = nexus_from_str(source, &read_options).unwrap();
    let refs: Vec<NexusTreeRef<'_>> = file.trees.iter().map(|tree| tree.as_ref()).collect();
    let options = NexusWriteOptions {
      newick: NewickWriteOptions::new(NewickDialect::Beast),
      translate,
    };

    let written = nexus_to_string(&refs, &options).unwrap();
    let reread = nexus_from_str(&written, &read_options).unwrap();

    let same = file.trees.iter().zip(&reread.trees).all(|(a, b)| a.name == b.name && a.comments == b.comments && a.tree.graph.eq_ordered(&b.tree.graph));
    assert_eq!((expected, true), (written.as_str(), same));
  }

  mod helpers {
    use crate::__tests__::test_read_basic::tests::helpers::summary;
    use crate::nexus::read::nexus_from_str;
    use crate::read::options::{NewickReadOptions, ReadMode};

    pub(super) fn trees(input: &str, options: &NewickReadOptions) -> Vec<(String, Vec<String>)> {
      match nexus_from_str(input, options) {
        Ok(file) => file
          .trees
          .iter()
          .map(|tree| (tree.name.clone(), summary(&tree.tree.graph)))
          .collect(),
        Err(error) => panic!("{input:?} does not read: {error}"),
      }
    }

    pub(super) fn read_error(input: &str, options: &NewickReadOptions) -> String {
      match nexus_from_str(input, options) {
        Ok(file) => panic!("{input:?} reads as {} trees", file.trees.len()),
        Err(error) => error.to_string(),
      }
    }

    pub(super) fn tolerant() -> NewickReadOptions {
      NewickReadOptions {
        mode: ReadMode::Tolerant,
        ..NewickReadOptions::default()
      }
    }
  }
}

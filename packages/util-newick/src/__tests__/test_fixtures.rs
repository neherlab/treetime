#[cfg(test)]
mod tests {
  use crate::__tests__::test_read_basic::tests::helpers::summary;
  use crate::dialect::NewickDialect;
  use crate::read::error::NewickErrorKind;
  use crate::read::options::NewickReadOptions;
  use crate::read::stream::newick_from_str;
  use crate::write::newick::newick_to_string;
  use crate::write::options::NewickWriteOptions;
  use helpers::{dendropy_view, ete_view, normalize_lines, phylonet_view, read_in};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  const BEAST1_TREE_LOG: &str = "[&R] ((1:[&rate=1.0]0.25,2:[&rate=0.9]0.25):[&rate=1.1]0.5,3:[&rate=1.0]0.75);";
  const BEAST2_NODE_AND_BRANCH_META: &str = "((A[&height=0.25]:[&rate=1.2]0.5,B[&height=0.25]:[&rate=0.8]0.5)[&height=0.75]:0.25,C[&height=0.0]:1.0)[&height=1.0];";
  const TREEANNOTATOR_MCC: &str = "((A[&height=0.0,height_95%_HPD={0.0,0.0},posterior=1.0]:1.5,B[&height=0.0,height_95%_HPD={0.0,0.0},posterior=1.0]:1.5)[&height=1.5,height_95%_HPD={1.2,1.8},posterior=0.9876]:0.5,C[&height=0.0,posterior=1.0]:2.0)[&height=2.0,posterior=1.0];";
  const MRBAYES_CONSENSUS: &str = "[&U] (A[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)=\"100\",prob+-sd=\"100+-0\"]:1.0[&length_mean=1.0e-01,length_median=1.0e-01,length_95%HPD={9.0e-02,1.1e-01}],B[&prob=1.00000000e+00]:0.2[&length_mean=2.0e-01],C[&prob=1.00000000e+00]:0.3[&length_mean=3.0e-01]);";
  const MRBAYES_SAMPLE: &str = "[&U] ((1:1.000000e-01[&B TK02Brlens 1.000000e-01],2:1.000000e-01[&B TK02Brlens 2.000000e-01]):2.000000e-01[&E ibr 2: 1.000000e-02 3.000000e-02],3:3.000000e-01);";
  const IQTREE_CONCORDANCE: &str = "(A:0.1,B:0.2,(C:0.3,D:0.4)[&gCF=\"33.33\",gDF1/gDF2=\"0/33.33\"]:0.5);";
  const FASTTREE_SH_SUPPORT: &str = "(A:0.1,B:0.2,(C:0.3,D:0.4)0.950:0.05);";
  const RAXML_NG_BOOTSTRAP: &str = "(A:0.1,B:0.2,(C:0.3,D:0.4)100:0.05);";
  const IQTREE_TWO_SUPPORTS: &str = "(A:0.1,B:0.2,(C:0.3,D:0.4)80.5/95:0.05);";
  const FORESTER_NHX: &str = "(((ADH2:0.1[&&NHX:S=human:E=1.1.1.1],ADH1:0.11[&&NHX:S=human:E=1.1.1.1]):0.05[&&NHX:S=Primates:E=1.1.1.1:D=Y:B=100],ADHY:0.1[&&NHX:S=nematode:E=1.1.1.1],ADHX:0.12[&&NHX:S=insect:E=1.1.1.1]):0.1[&&NHX:S=Metazoa:E=1.1.1.1:D=N],(ADH4:0.09[&&NHX:S=yeast:E=1.1.1.1],ADH3:0.13[&&NHX:S=yeast:E=1.1.1.1],ADH2:0.12[&&NHX:S=yeast:E=1.1.1.1],ADH1:0.11[&&NHX:S=yeast:E=1.1.1.1]):0.1[&&NHX:S=Fungi])[&&NHX:E=1.1.1.1:D=N];";
  const CARDONA_ENEWICK: &str = "(A,B,((C,(Y)x#H1)c,(x#H1,D)d)e)f;";
  const DENDROSCOPE_LGT_ACCEPTOR: &str = "((A,(B)x##LGT1)c,(x#LGT1,C)d)r;";
  const SPLITSTREE_HYBRID: &str = "((a:1,(b:1)#H1:1)e:1,(#H1:1,c:1)f:1)g;";
  const PHYLONET_RICH_GAMMA: &str = "[&U]((A:1,(B:1)#H1:1::0.3):1,(#H1:1::0.7,C:1):1);";
  const PHYLONET_SUPPORT_FIELDS: &str = "((A:1:90,B:1:80):1:100,C:2);";
  const PHYLONET_RICH_WEIGHT: &str = "[&R][&W 0.5]((A:1:90,B:1:80):1:100,C:2);";

  #[rustfmt::skip]
  #[rstest]
  #[case::dendropy_5_0_8_beast1_tree_log(              BEAST1_TREE_LOG,              NewickDialect::Beast, indoc! {"
    -
    - :0.5 rate=1.1
    1 :0.25 rate=1.0
    2 :0.25 rate=0.9
    3 :0.75 rate=1.0"})]
  #[case::dendropy_5_0_8_beast2_node_and_branch_meta(  BEAST2_NODE_AND_BRANCH_META,  NewickDialect::Beast, indoc! {"
    - height=1.0
    - :0.25 height=0.75
    A :0.5 height=0.25 rate=1.2
    B :0.5 height=0.25 rate=0.8
    C :1.0 height=0.0"})]
  #[case::dendropy_5_0_8_treeannotator_mcc(            TREEANNOTATOR_MCC,            NewickDialect::Beast, indoc! {"
    - height=2.0 posterior=1.0
    - :0.5 height=1.5 height_95%_HPD={1.2,1.8} posterior=0.9876
    A :1.5 height=0.0 height_95%_HPD={0.0,0.0} posterior=1.0
    B :1.5 height=0.0 height_95%_HPD={0.0,0.0} posterior=1.0
    C :2.0 height=0.0 posterior=1.0"})]
  #[case::dendropy_5_0_8_mrbayes_consensus(            MRBAYES_CONSENSUS,            NewickDialect::Beast, indoc! {"
    -
    A :1.0 length_95%HPD={9.0e-02,1.1e-01} length_mean=1.0e-01 length_median=1.0e-01 prob(percent)=100 prob+-sd=100+-0 prob=1.00000000e+00 prob_range={1.00000000e+00,1.00000000e+00} prob_stddev=0.00000000e+00
    B :0.2 length_mean=2.0e-01 prob=1.00000000e+00
    C :0.3 length_mean=3.0e-01 prob=1.00000000e+00"})]
  #[case::dendropy_5_0_8_iqtree_concordance(           IQTREE_CONCORDANCE,           NewickDialect::Beast, indoc! {"
    -
    A :0.1
    B :0.2
    - :0.5 gCF=33.33 gDF1/gDF2=0/33.33
    C :0.3
    D :0.4"})]
  #[case::dendropy_5_0_8_fasttree_sh_support(          FASTTREE_SH_SUPPORT,          NewickDialect::Classic, indoc! {"
    -
    A :0.1
    B :0.2
    0.950 :0.05
    C :0.3
    D :0.4"})]
  #[case::dendropy_5_0_8_raxml_ng_bootstrap(           RAXML_NG_BOOTSTRAP,           NewickDialect::Classic, indoc! {"
    -
    A :0.1
    B :0.2
    100 :0.05
    C :0.3
    D :0.4"})]
  #[case::dendropy_5_0_8_iqtree_two_supports(          IQTREE_TWO_SUPPORTS,          NewickDialect::Classic, indoc! {"
    -
    A :0.1
    B :0.2
    80.5/95 :0.05
    C :0.3
    D :0.4"})]
  #[trace]
  fn test_fixtures_match_dendropy(#[case] input: &str, #[case] dialect: NewickDialect, #[case] expected: &str) {
    let graph = read_in(input, dialect);

    assert_eq!(normalize_lines(expected), normalize_lines(&dendropy_view(&graph)));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::ete3_3_1_3_forester_nhx(FORESTER_NHX, indoc! {"
    - :0.0 D=N E=1.1.1.1
    - :0.1 D=N E=1.1.1.1 S=Metazoa
    - :0.05 B=100 D=Y E=1.1.1.1 S=Primates
    ADH2 :0.1 E=1.1.1.1 S=human
    ADH1 :0.11 E=1.1.1.1 S=human
    ADHY :0.1 E=1.1.1.1 S=nematode
    ADHX :0.12 E=1.1.1.1 S=insect
    - :0.1 S=Fungi
    ADH4 :0.09 E=1.1.1.1 S=yeast
    ADH3 :0.13 E=1.1.1.1 S=yeast
    ADH2 :0.12 E=1.1.1.1 S=yeast
    ADH1 :0.11 E=1.1.1.1 S=yeast"})]
  #[trace]
  fn test_fixtures_match_ete3(#[case] input: &str, #[case] expected: &str) {
    let graph = read_in(input, NewickDialect::Nhx);

    assert_eq!(expected, ete_view(&graph));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::phylonet_3_8_5_cardona_enewick(      CARDONA_ENEWICK,          NewickDialect::ENewick, indoc! {"
    f
    A | from f length=-Infinity support=-Infinity probability=-Infinity
    B | from f length=-Infinity support=-Infinity probability=-Infinity
    e | from f length=-Infinity support=-Infinity probability=-Infinity
    c | from e length=-Infinity support=-Infinity probability=-Infinity
    C | from c length=-Infinity support=-Infinity probability=-Infinity
    x network | from c length=-Infinity support=-Infinity probability=-Infinity | from d length=-Infinity support=-Infinity probability=-Infinity
    Y | from x length=-Infinity support=-Infinity probability=-Infinity
    d | from e length=-Infinity support=-Infinity probability=-Infinity
    D | from d length=-Infinity support=-Infinity probability=-Infinity"})]
  #[case::phylonet_3_8_5_splitstree_hybrid(    SPLITSTREE_HYBRID,        NewickDialect::ENewick, indoc! {"
    g
    e | from g length=1.0 support=-Infinity probability=-Infinity
    a | from e length=1.0 support=-Infinity probability=-Infinity
    - network | from e length=1.0 support=-Infinity probability=-Infinity | from f length=1.0 support=-Infinity probability=-Infinity
    b | from - length=1.0 support=-Infinity probability=-Infinity
    f | from g length=1.0 support=-Infinity probability=-Infinity
    c | from f length=1.0 support=-Infinity probability=-Infinity"})]
  #[case::phylonet_3_8_5_rich_gamma(           PHYLONET_RICH_GAMMA,      NewickDialect::Rich,    indoc! {"
    -
    - | from - length=1.0 support=-Infinity probability=-Infinity
    A | from - length=1.0 support=-Infinity probability=-Infinity
    - network | from - length=1.0 support=-Infinity probability=0.3 | from - length=1.0 support=-Infinity probability=0.7
    B | from - length=1.0 support=-Infinity probability=-Infinity
    - | from - length=1.0 support=-Infinity probability=-Infinity
    C | from - length=1.0 support=-Infinity probability=-Infinity"})]
  #[case::phylonet_3_8_5_support_fields(       PHYLONET_SUPPORT_FIELDS,  NewickDialect::Rich,    indoc! {"
    -
    - | from - length=1.0 support=100.0 probability=-Infinity
    A | from - length=1.0 support=90.0 probability=-Infinity
    B | from - length=1.0 support=80.0 probability=-Infinity
    C | from - length=2.0 support=-Infinity probability=-Infinity"})]
  #[trace]
  fn test_fixtures_match_phylonet(#[case] input: &str, #[case] dialect: NewickDialect, #[case] expected: &str) {
    let graph = read_in(input, dialect);

    assert_eq!(expected, phylonet_view(&graph));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::spec_mrbayes_sample(          MRBAYES_SAMPLE,           NewickDialect::MrBayes, vec!["- ", "- :0.2", "1 :0.1", "2 :0.1", "3 :0.3"])]
  #[case::spec_dendroscope_acceptor(    DENDROSCOPE_LGT_ACCEPTOR, NewickDialect::ENewick, vec!["r ", "c ", "A ", "x#LGT1 acceptor | ", "B ", "d ", "C "])]
  #[case::spec_phylonet_rich_weight(    PHYLONET_RICH_WEIGHT,     NewickDialect::Rich,    vec!["- ", "- :1 support=100(field)", "A :1 support=90(field)", "B :1 support=80(field)", "C :2"])]
  #[trace]
  fn test_fixtures_spec_expectations(#[case] input: &str, #[case] dialect: NewickDialect, #[case] expected: Vec<&str>) {
    assert_eq!(expected, summary(&read_in(input, dialect)));
  }

  #[test]
  fn test_fixtures_spec_rich_weight_and_rooting() {
    let graph = read_in(PHYLONET_RICH_WEIGHT, NewickDialect::Rich);

    assert_eq!((Some(true), Some(0.5)), (graph.rooted(), graph.weight()));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::beast1(          BEAST1_TREE_LOG,              NewickDialect::Beast)]
  #[case::beast2(          BEAST2_NODE_AND_BRANCH_META,  NewickDialect::Beast)]
  #[case::treeannotator(   TREEANNOTATOR_MCC,            NewickDialect::Beast)]
  #[case::mrbayes_con(     MRBAYES_CONSENSUS,            NewickDialect::Beast)]
  #[case::mrbayes_sample(  MRBAYES_SAMPLE,               NewickDialect::MrBayes)]
  #[case::iqtree_cf(       IQTREE_CONCORDANCE,           NewickDialect::Beast)]
  #[case::fasttree(        FASTTREE_SH_SUPPORT,          NewickDialect::Classic)]
  #[case::raxml_ng(        RAXML_NG_BOOTSTRAP,           NewickDialect::Classic)]
  #[case::iqtree_supports( IQTREE_TWO_SUPPORTS,          NewickDialect::Classic)]
  #[case::forester(        FORESTER_NHX,                 NewickDialect::Nhx)]
  #[case::cardona(         CARDONA_ENEWICK,              NewickDialect::ENewick)]
  #[case::dendroscope(     DENDROSCOPE_LGT_ACCEPTOR,     NewickDialect::ENewick)]
  #[case::splitstree(      SPLITSTREE_HYBRID,            NewickDialect::ENewick)]
  #[case::phylonet_gamma(  PHYLONET_RICH_GAMMA,          NewickDialect::Rich)]
  #[case::phylonet_weight( PHYLONET_RICH_WEIGHT,         NewickDialect::Rich)]
  #[trace]
  fn test_fixtures_roundtrip_in_own_dialect(#[case] input: &str, #[case] dialect: NewickDialect) {
    let graph = read_in(input, dialect);

    let written = newick_to_string(&graph, &NewickWriteOptions::new(dialect)).unwrap();

    assert!(graph.eq_ordered(&read_in(&written, dialect)), "{written} reads back differently");
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::beast1_all(          BEAST1_TREE_LOG,          NewickDialect::ALL.to_vec(),        Ok(NewickDialect::Beast))]
  #[case::beast1_classic(      BEAST1_TREE_LOG,          vec![NewickDialect::Classic],       Ok(NewickDialect::Classic))]
  #[case::beast1_nhx(          BEAST1_TREE_LOG,          vec![NewickDialect::Nhx],           Err(NewickErrorKind::Annotation))]
  #[case::mrbayes_con_all(     MRBAYES_CONSENSUS,        NewickDialect::ALL.to_vec(),        Ok(NewickDialect::Beast))]
  #[case::mrbayes_sample_all(  MRBAYES_SAMPLE,           NewickDialect::ALL.to_vec(),        Ok(NewickDialect::MrBayes))]
  #[case::mrbayes_sample_beast(MRBAYES_SAMPLE,           vec![NewickDialect::Beast],         Err(NewickErrorKind::Annotation))]
  #[case::forester_all(        FORESTER_NHX,             NewickDialect::ALL.to_vec(),        Ok(NewickDialect::Nhx))]
  #[case::forester_beast(      FORESTER_NHX,             vec![NewickDialect::Beast],         Err(NewickErrorKind::Annotation))]
  #[case::fasttree_all(        FASTTREE_SH_SUPPORT,      NewickDialect::ALL.to_vec(),        Ok(NewickDialect::Rich))]
  #[case::splitstree_all(      SPLITSTREE_HYBRID,        NewickDialect::ALL.to_vec(),        Ok(NewickDialect::Rich))]
  #[case::splitstree_classic(  SPLITSTREE_HYBRID,        vec![NewickDialect::Classic],       Ok(NewickDialect::Classic))]
  #[case::phylonet_all(        PHYLONET_RICH_GAMMA,      NewickDialect::ALL.to_vec(),        Ok(NewickDialect::Rich))]
  #[case::phylonet_enewick(    PHYLONET_RICH_GAMMA,      vec![NewickDialect::ENewick],       Err(NewickErrorKind::Syntax))]
  #[case::quoted_bracket_classic("(A[&a=\"x]y\"],B);",   vec![NewickDialect::Classic],       Err(NewickErrorKind::Syntax))]
  #[trace]
  fn test_fixtures_dialect_selection(#[case] input: &str, #[case] dialects: Vec<NewickDialect>, #[case] expected: Result<NewickDialect, NewickErrorKind>) {
    let options = NewickReadOptions {
      dialects,
      ..NewickReadOptions::default()
    };

    let actual = newick_from_str(input, &options).map(|tree| tree.dialect).map_err(|error| error.kind);

    assert_eq!(expected, actual);
  }

  mod helpers {
    use crate::__tests__::test_read_basic::tests::helpers::read_with;
    use crate::dialect::NewickDialect;
    use crate::model::data::NewickEdgeData;
    use crate::model::graph::NewickGraph;
    use crate::model::value::NewickValue;
    use crate::number::format_shortest;
    use crate::read::options::NewickReadOptions;
    use std::fmt::Write;

    pub(super) fn read_in(text: &str, dialect: NewickDialect) -> NewickGraph {
      let options = NewickReadOptions {
        dialects: vec![dialect],
        ..NewickReadOptions::default()
      };
      read_with(text, &options).graph
    }

    pub(super) fn dendropy_view(graph: &NewickGraph) -> String {
      graph
        .preorder()
        .map(|node| {
          let edge = edge_above(graph, node);
          let label = match (graph.node(node).name(), edge.support()) {
            (Some(name), _) => name.to_owned(),
            (None, []) => "-".to_owned(),
            (None, support) => support
              .iter()
              .map(|value| format_shortest(*value).unwrap())
              .collect::<Vec<_>>()
              .join("/"),
          };
          let mut parts = vec![label];
          if let Some(length) = edge.branch_length() {
            parts.push(format!(":{length:?}"));
          }
          let mut notes = annotations(graph, node);
          notes.sort();
          parts.extend(notes);
          parts.join(" ")
        })
        .collect::<Vec<_>>()
        .join("\n")
    }

    pub(super) fn ete_view(graph: &NewickGraph) -> String {
      graph
        .preorder()
        .map(|node| {
          let default_length = if node == graph.root() { 0.0 } else { 1.0 };
          let length = edge_above(graph, node).branch_length().unwrap_or(default_length);
          let mut parts = vec![
            graph.node(node).name().unwrap_or("-").to_owned(),
            format!(":{length:?}"),
          ];
          let mut notes = annotations(graph, node);
          notes.sort();
          parts.extend(notes);
          parts.join(" ")
        })
        .collect::<Vec<_>>()
        .join("\n")
    }

    pub(super) fn phylonet_view(graph: &NewickGraph) -> String {
      let java_double = |value: Option<f64>| value.map_or_else(|| "-Infinity".to_owned(), |value| format!("{value:?}"));
      graph
        .preorder()
        .map(|node| {
          let data = graph.node(node);
          let mut text = data.name().unwrap_or("-").to_owned();
          if data.hybrid().is_some() {
            text.push_str(" network");
          }
          for &edge in graph.parent_edges(node) {
            let entry = graph.edge(edge);
            let edge_data = entry.data();
            write!(
              text,
              " | from {} length={} support={} probability={}",
              graph.node(entry.parent()).name().unwrap_or("-"),
              java_double(edge_data.branch_length()),
              java_double(edge_data.support().first().copied()),
              java_double(edge_data.probability()),
            )
            .unwrap();
          }
          text
        })
        .collect::<Vec<_>>()
        .join("\n")
    }

    pub(super) fn normalize_lines(text: &str) -> Vec<Vec<String>> {
      text
        .lines()
        .map(|line| {
          line
            .split(' ')
            .enumerate()
            .map(|(idx, token)| {
              if idx == 0 {
                normalize_label(token)
              } else {
                token.to_owned()
              }
            })
            .collect()
        })
        .collect()
    }

    fn normalize_label(label: &str) -> String {
      let numbers: Result<Vec<f64>, _> = label.split('/').map(str::parse::<f64>).collect();
      match numbers {
        Ok(values) if label != "-" => values
          .iter()
          .map(|value| format_shortest(*value).unwrap())
          .collect::<Vec<_>>()
          .join("/"),
        Ok(_) | Err(_) => label.to_owned(),
      }
    }

    fn edge_above(graph: &NewickGraph, node: usize) -> &NewickEdgeData {
      match graph.parent_edges(node).first() {
        Some(&edge) => graph.edge(edge).data(),
        None => graph.root_edge(),
      }
    }

    fn annotations(graph: &NewickGraph, node: usize) -> Vec<String> {
      graph
        .node(node)
        .annotations()
        .chain(edge_above(graph, node).annotations())
        .map(|(key, value)| format!("{key}={}", value_text(value)))
        .collect()
    }

    fn value_text(value: &NewickValue) -> String {
      match value {
        NewickValue::Boolean(value) => value.to_string(),
        NewickValue::Number(number) => format_shortest(*number).unwrap(),
        NewickValue::NumberText(text) | NewickValue::String(text) => text.clone(),
        NewickValue::Color([red, green, blue]) => format!("{red}.{green}.{blue}"),
        NewickValue::Array(values) => format!(
          "{{{}}}",
          values.as_slice().iter().map(value_text).collect::<Vec<_>>().join(",")
        ),
      }
    }
  }
}

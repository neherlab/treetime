#[cfg(test)]
mod tests {
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;

  use helpers::{DatesCase, dates_table};

  const TREE: &str = "((A:0.1,B:0.2)C:0.1,D:0.3)root;";

  #[test]
  fn test_dates_table_writes_date_and_numeric_date_of_every_node() -> Result<(), Report> {
    let actual = dates_table(&DatesCase::default())?;
    assert_eq!(
      indoc! {"
        #node\tdate\tnumeric date
        root\t2010-01-01\t2010.001370
        C\t2012-07-02\t2012.500000
        A\t2016-06-13\t2016.450000
        B\t2015-12-31\t2015.999000
        D\t2016-02-01\t2016.085000
      "},
      actual
    );
    Ok(())
  }

  #[test]
  fn test_dates_table_adds_comment_and_bounds_with_confidence() -> Result<(), Report> {
    let actual = dates_table(&DatesCase {
      confidence: true,
      ..DatesCase::default()
    })?;
    assert_eq!(
      indoc! {"
        #Lower and upper bound delineate the 90% max posterior region
        #node\tdate\tnumeric date\tlower bound\tupper bound
        root\t2010-01-01\t2010.001370\t2008.500000\t2011.250000
        C\t2012-07-02\t2012.500000\t2012.000000\t2013.000000
        A\t2016-06-13\t2016.450000\t2016.450000\t2016.450000
        B\t2015-12-31\t2015.999000\t2015.999000\t2015.999000
        D\t2016-02-01\t2016.085000\t2016.085000\t2016.085000
      "},
      actual
    );
    Ok(())
  }

  #[test]
  fn test_dates_table_writes_dashes_for_nodes_on_bad_branches() -> Result<(), Report> {
    let case = DatesCase {
      excluded: &["B", "D"],
      ..DatesCase::default()
    };
    let actual = (
      dates_table(&case)?,
      dates_table(&DatesCase {
        confidence: true,
        ..case
      })?,
    );
    assert_eq!(
      (
        indoc! {"
          #node\tdate\tnumeric date
          root\t2010-01-01\t2010.001370
          C\t2012-07-02\t2012.500000
          A\t2016-06-13\t2016.450000
          B\t--\t--
          D\t--\t--
        "}
        .to_owned(),
        indoc! {"
          #Lower and upper bound delineate the 90% max posterior region
          #node\tdate\tnumeric date\tlower bound\tupper bound
          root\t2010-01-01\t2010.001370\t2008.500000\t2011.250000
          C\t2012-07-02\t2012.500000\t2012.000000\t2013.000000
          A\t2016-06-13\t2016.450000\t2016.450000\t2016.450000
          B\t--\t--\t--\t--
          D\t--\t--\t--\t--
        "}
        .to_owned()
      ),
      actual
    );
    Ok(())
  }

  #[test]
  fn test_dates_table_writes_dashes_for_nodes_without_a_date() -> Result<(), Report> {
    let actual = dates_table(&DatesCase {
      undated: &["C"],
      ..DatesCase::default()
    })?;
    assert_eq!(
      indoc! {"
        #node\tdate\tnumeric date
        root\t2010-01-01\t2010.001370
        C\t--\t--
        A\t2016-06-13\t2016.450000
        B\t2015-12-31\t2015.999000
        D\t2016-02-01\t2016.085000
      "},
      actual
    );
    Ok(())
  }

  #[test]
  fn test_dates_table_writes_the_date_of_a_node_without_bounds() -> Result<(), Report> {
    let actual = dates_table(&DatesCase {
      confidence: true,
      unbounded: &["C"],
      ..DatesCase::default()
    })?;
    assert_eq!(
      indoc! {"
        #Lower and upper bound delineate the 90% max posterior region
        #node\tdate\tnumeric date\tlower bound\tupper bound
        root\t2010-01-01\t2010.001370\t2008.500000\t2011.250000
        C\t2012-07-02\t2012.500000\t--\t--
        A\t2016-06-13\t2016.450000\t2016.450000\t2016.450000
        B\t2015-12-31\t2015.999000\t2015.999000\t2015.999000
        D\t2016-02-01\t2016.085000\t2016.085000\t2016.085000
      "},
      actual
    );
    Ok(())
  }

  #[test]
  fn test_dates_table_writes_commas_for_a_csv_path() -> Result<(), Report> {
    let actual = dates_table(&DatesCase {
      confidence: true,
      commas: true,
      ..DatesCase::default()
    })?;
    assert_eq!(
      indoc! {"
        #Lower and upper bound delineate the 90% max posterior region
        #node,date,numeric date,lower bound,upper bound
        root,2010-01-01,2010.001370,2008.500000,2011.250000
        C,2012-07-02,2012.500000,2012.000000,2013.000000
        A,2016-06-13,2016.450000,2016.450000,2016.450000
        B,2015-12-31,2015.999000,2015.999000,2015.999000
        D,2016-02-01,2016.085000,2016.085000,2016.085000
      "},
      actual
    );
    Ok(())
  }

  mod helpers {
    use super::TREE;
    use crate::__tests__::test_support::tests::{Topology, annotated, node_key, topology_from};
    use crate::annotated_graph::{AnnotatedGraph, AnnotatedTreeView, TreeDates};
    use crate::dates_table::write_dates_table;
    use eyre::Report;
    use maplit::btreemap;
    use std::collections::{BTreeMap, BTreeSet};
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::csv::TableFormat;

    #[derive(Default)]
    pub(super) struct DatesCase {
      pub(super) confidence: bool,
      pub(super) excluded: &'static [&'static str],
      pub(super) undated: &'static [&'static str],
      pub(super) unbounded: &'static [&'static str],
      pub(super) commas: bool,
    }

    pub(super) fn dates_table(case: &DatesCase) -> Result<String, Report> {
      let topology = topology_from(TREE)?;
      let key = |name: &str| node_key(&topology, name);
      let num_date: BTreeMap<GraphNodeKey, Option<f64>> = btreemap! {
        key("root") => 2010.00137,
        key("C") => 2012.5,
        key("A") => 2016.45,
        key("B") => 2015.999,
        key("D") => 2016.085,
      }
      .into_iter()
      .map(|(node, date)| {
        (
          node,
          (!case.undated.iter().any(|&name| key(name) == node)).then_some(date),
        )
      })
      .collect();
      let confidence: BTreeMap<GraphNodeKey, [f64; 2]> = btreemap! {
        key("root") => [2008.5, 2011.25],
        key("C") => [2012.0, 2013.0],
        key("A") => [2016.45, 2016.45],
        key("B") => [2015.999, 2015.999],
        key("D") => [2016.085, 2016.085],
      }
      .into_iter()
      .filter(|(node, _)| !case.unbounded.iter().any(|&name| key(name) == *node))
      .collect();
      let excluded: BTreeSet<GraphNodeKey> = case.excluded.iter().map(|&name| key(name)).collect();
      let format = if case.commas {
        TableFormat::Csv
      } else {
        TableFormat::Tsv
      };
      write_table(
        &topology,
        &num_date,
        case.confidence.then_some(&confidence),
        &excluded,
        format,
      )
    }

    fn write_table(
      topology: &Topology,
      num_date: &BTreeMap<GraphNodeKey, Option<f64>>,
      confidence: Option<&BTreeMap<GraphNodeKey, [f64; 2]>>,
      excluded: &BTreeSet<GraphNodeKey>,
      format: TableFormat,
    ) -> Result<String, Report> {
      let graph = AnnotatedGraph {
        dates: Some(TreeDates {
          num_date,
          confidence,
          excluded,
          input_dates: None,
        }),
        ..annotated(topology)
      };
      let mut buf = Vec::new();
      write_dates_table(&AnnotatedTreeView::new(&graph)?, &mut buf, format)?;
      Ok(String::from_utf8(buf)?)
    }
  }
}

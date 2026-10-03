#[cfg(test)]
mod tests {
  use crate::clock::clock_model::ClockModel;
  use crate::test_utils::{RecordingLog, find_node_key_by_name};
  use crate::timetree::optimization::outliers::report_outliers;
  use eyre::Report;
  use helpers::{header, row};
  use maplit::{btreemap, btreeset};
  use pretty_assertions::assert_eq;
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_report_outliers_lists_dated_outliers_by_descending_residual() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:1,B:2,C:3,D:4,E:5,F:6)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("named node exists");
    let divergences = btreemap! {
      key("root") => 0.0,
      key("A") => 1.0,
      key("B") => 2.0,
      key("C") => 3.0,
      key("D") => 4.0,
      key("E") => 5.0,
      key("F") => 6.0,
    };
    let given_dates = btreemap! {
      key("A") => Some(2010.0),
      key("B") => Some(2000.0),
      key("C") => Some(2006.0),
      key("D") => Some(1990.0),
      key("E") => Some(2018.0),
      key("F") => None,
    };
    let outliers = btreeset! { key("A"), key("B"), key("D"), key("E"), key("F") };
    let log = RecordingLog::default();

    report_outliers(
      &graph,
      &outliers,
      &divergences,
      &ClockModel::for_testing(0.5, -1000.0),
      2.0,
      &given_dates,
      &names,
      &log,
    );

    let expected = vec![
      "Clock filter marked 4 outliers:".to_owned(),
      header(),
      row("D", 1990.0, 2008.0, -4.5),
      row("E", 2018.0, 2010.0, 2.0),
      row("A", 2010.0, 2002.0, 2.0),
      row("B", 2000.0, 2004.0, -1.0),
    ];
    assert_eq!(expected, log.warnings());
    Ok(())
  }

  #[test]
  fn test_report_outliers_is_silent_without_outliers() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:1,B:2)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let log = RecordingLog::default();

    report_outliers(
      &graph,
      &btreeset! {},
      &btreemap! {},
      &ClockModel::for_testing(0.5, -1000.0),
      2.0,
      &btreemap! {},
      &names,
      &log,
    );

    assert_eq!(Vec::<String>::new(), log.warnings());
    Ok(())
  }

  mod helpers {
    pub(super) fn header() -> String {
      format!(
        "{:>20} {:>12} {:>14} {:>10}",
        "name", "given_date", "apparent_date", "residual"
      )
    }

    pub(super) fn row(name: &str, given_date: f64, apparent_date: f64, residual: f64) -> String {
      format!("{name:>20} {given_date:>12.2} {apparent_date:>14.2} {residual:>10.2}")
    }
  }
}

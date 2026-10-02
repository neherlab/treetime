#[cfg(test)]
mod tests {
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_state::ClockInputs;
  use crate::clock::rtt::{ClockRegressionResult, gather_clock_regression_results};
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use maplit::{btreemap, btreeset};
  use std::collections::BTreeMap;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::pretty_assert_eq;

  use helpers::row;

  #[test]
  fn test_gather_clock_regression_results_reports_every_node_against_the_clock_line() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:1,B:2)AB:1,C:4)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("named node exists");
    let divergences = btreemap! {
      key("root") => 0.0,
      key("AB") => 1.0,
      key("A") => 2.0,
      key("B") => 3.0,
      key("C") => 4.0,
    };
    let times = btreemap! { key("A") => Some(2005.0), key("B") => None, key("C") => Some(2007.0) };
    let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
    let clock_model = ClockModel::for_testing(0.5, -1000.0);

    let mut actual = gather_clock_regression_results(
      &graph,
      &inputs,
      &divergences,
      &btreeset! { key("C") },
      &clock_model,
      &names,
    );
    actual.sort_by(|lhs, rhs| lhs.name.cmp(&rhs.name));

    let expected: Vec<ClockRegressionResult> = vec![
      row("A", 2.0, Some(2005.0), 2004.0, Some(0.5), false, true),
      row("AB", 1.0, None, 2002.0, None, false, false),
      row("B", 3.0, None, 2006.0, None, false, true),
      row("C", 4.0, Some(2007.0), 2008.0, Some(-0.5), true, true),
      row("root", 0.0, None, 2000.0, None, false, false),
    ];
    pretty_assert_eq!(expected, actual);
    Ok(())
  }

  mod helpers {
    use crate::clock::rtt::ClockRegressionResult;
    use crate::o;

    pub(super) fn row(
      name: &str,
      div: f64,
      date: Option<f64>,
      predicted_date: f64,
      clock_deviation: Option<f64>,
      is_outlier: bool,
      is_leaf: bool,
    ) -> ClockRegressionResult {
      ClockRegressionResult {
        name: Some(o!(name)),
        div,
        date,
        predicted_date,
        clock_deviation,
        is_outlier,
        is_leaf,
        date_source: None,
      }
    }
  }
}

#[cfg(test)]
mod tests {
  use crate::progress::{RunWarning, RunWarningKind};
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::dates_outside_input::dates_outside_input_dates;
  use eyre::Report;
  use treetime_graph::keymap;

  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::date::{DateConstraint, DateRange, DateValue};
  use treetime_utils::{o, vec_of_owned};

  const TREE: &str = "((A:1,B:1)X:1,(C:1,D:1)Y:1)root;";

  #[rustfmt::skip]
  #[rstest]
  #[case::inside(              [2010.0, 2011.2, 2012.5, 2012.0, 2005.0, 2008.0, 2000.0], None)]
  #[case::tip_before_exact(    [2008.5, 2011.2, 2012.5, 2012.0, 2005.0, 2008.0, 2000.0], Some(("Some final dates lie outside the input dates: 1 tip dated outside the input date range of the tip, by up to 1.50 years (A). The clock alone dates these nodes, as it does for the tips that the clock filter flagged and for their ancestors", vec_of_owned!["A"])))]
  #[case::tip_after_range(     [2010.0, 2012.0, 2012.5, 2012.0, 2005.0, 2008.0, 2000.0], Some(("Some final dates lie outside the input dates: 1 tip dated outside the input date range of the tip, by up to 0.50 years (B). The clock alone dates these nodes, as it does for the tips that the clock filter flagged and for their ancestors", vec_of_owned!["B"])))]
  #[case::internal_after_last( [2010.0, 2011.2, 2012.5, 2012.0, 2013.25, 2008.0, 2000.0], Some(("Some final dates lie outside the input dates: 1 internal node dated after the latest input date 2013.00, by up to 0.25 years (X). The clock alone dates these nodes, as it does for the tips that the clock filter flagged and for their ancestors", vec_of_owned!["X"])))]
  #[case::both(                [2014.0, 2015.0, 2012.5, 2016.0, 2017.0, 2018.0, 2000.0], Some(("Some final dates lie outside the input dates: 3 tips dated outside the input date range of the tip, by up to 4.00 years (A, B, D); 2 internal nodes dated after the latest input date 2013.00, by up to 5.00 years (X, Y). The clock alone dates these nodes, as it does for the tips that the clock filter flagged and for their ancestors", vec_of_owned!["A", "B", "D", "X", "Y"])))]
  #[trace]
  fn test_dates_outside_input_dates(
    #[case] node_dates: [f64; 7],
    #[case] expected: Option<(&str, Vec<String>)>,
  ) -> Result<(), Report> {
    let [a, b, c, d, x, y, root] = node_dates;
    let parsed = nwk_read(TREE.as_bytes())?;
    let names = parsed.names.clone();
    let key = |name: &str| find_node_key_by_name(&parsed.graph, &names, name).expect("node must exist");
    let range = |start: f64, end: f64| DateConstraint {
      raw: o!("range"),
      value: DateValue::Range(DateRange { start, end }),
    };
    let dates = keymap! {
      key("A") => DateConstraint::exact(2010.0),
      key("B") => range(2011.0, 2011.5),
      key("C") => range(2012.0, 2013.0),
      key("D") => DateConstraint::exact(2012.0),
    };
    let node_dates = keymap! {
      key("A") => Some(a),
      key("B") => Some(b),
      key("C") => Some(c),
      key("D") => Some(d),
      key("X") => Some(x),
      key("Y") => Some(y),
      key("root") => Some(root),
    };

    let actual = dates_outside_input_dates(&parsed.graph, &names, &dates, &node_dates);

    let expected = expected.map(|(message, names)| RunWarning {
      kind: RunWarningKind::DatesOutsideInputDates,
      message: message.to_owned(),
      names,
    });
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_dates_outside_input_dates_skips_undated_tips_and_nodes_without_dates() -> Result<(), Report> {
    let parsed = nwk_read(TREE.as_bytes())?;
    let names = parsed.names.clone();
    let key = |name: &str| find_node_key_by_name(&parsed.graph, &names, name).expect("node must exist");
    let dates = keymap! {
      key("A") => DateConstraint::exact(2010.0),
    };
    let node_dates = keymap! {
      key("A") => Some(2010.0),
      key("B") => Some(2009.0),
      key("C") => Some(2009.0),
      key("D") => None,
      key("X") => None,
      key("Y") => Some(2005.0),
      key("root") => Some(2000.0),
    };

    assert_eq!(
      None,
      dates_outside_input_dates(&parsed.graph, &names, &dates, &node_dates)
    );
    Ok(())
  }
}

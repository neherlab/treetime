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
  #[case::tip_before_exact(    [2008.5, 2011.2, 2012.5, 2012.0, 2005.0, 2008.0, 2000.0], Some(("Some final dates lie outside the input dates: 1 tip dated outside the input date range of the tip, by up to 1.50 years (A at 2008.50). The clock alone dates these nodes, as it does for the tips that the clock filter flagged and for their ancestors", vec_of_owned!["A"])))]
  #[case::tip_after_range(     [2010.0, 2012.0, 2012.5, 2012.0, 2005.0, 2008.0, 2000.0], Some(("Some final dates lie outside the input dates: 1 tip dated outside the input date range of the tip, by up to 0.50 years (B at 2012.00). The clock alone dates these nodes, as it does for the tips that the clock filter flagged and for their ancestors", vec_of_owned!["B"])))]
  #[case::internal_after_last( [2010.0, 2011.2, 2012.5, 2012.0, 2013.25, 2008.0, 2000.0], Some(("Some final dates lie outside the input dates: 1 internal node dated after the latest input date 2013.00, by up to 0.25 years (X at 2013.25). The clock alone dates these nodes, as it does for the tips that the clock filter flagged and for their ancestors", vec_of_owned!["X"])))]
  #[case::both(                [2014.0, 2015.0, 2012.5, 2016.0, 2017.0, 2018.0, 2000.0], Some(("Some final dates lie outside the input dates: 3 tips dated outside the input date range of the tip, by up to 4.00 years (A at 2014.00, B at 2015.00, D at 2016.00); 2 internal nodes dated after the latest input date 2013.00, by up to 5.00 years (X at 2017.00, Y at 2018.00). The clock alone dates these nodes, as it does for the tips that the clock filter flagged and for their ancestors", vec_of_owned!["A", "B", "D", "X", "Y"])))]
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

  #[rustfmt::skip]
  #[rstest]
  #[case::inside(    [2011.0, 2012.0], None)]
  #[case::at_bounds( [2010.0, 2012.0], None)]
  #[case::after(     [2254.4, 2011.0], Some(("Some final dates lie outside the input dates: 1 undated tip dated outside the sampling period 2010.00 to 2012.00, by up to 242.40 years (B at 2254.40). The clock alone dates these nodes, as it does for the tips that the clock filter flagged and for their ancestors", vec_of_owned!["B"])))]
  #[case::before(    [2011.0, 2009.5], Some(("Some final dates lie outside the input dates: 1 undated tip dated outside the sampling period 2010.00 to 2012.00, by up to 0.50 years (C at 2009.50). The clock alone dates these nodes, as it does for the tips that the clock filter flagged and for their ancestors", vec_of_owned!["C"])))]
  #[trace]
  fn test_dates_outside_input_dates_checks_undated_tips_against_the_sampling_period(
    #[case] undated_dates: [f64; 2],
    #[case] expected: Option<(&str, Vec<String>)>,
  ) -> Result<(), Report> {
    let [b, c] = undated_dates;
    let parsed = nwk_read(TREE.as_bytes())?;
    let names = parsed.names.clone();
    let key = |name: &str| find_node_key_by_name(&parsed.graph, &names, name).expect("node must exist");
    let dates = keymap! {
      key("A") => DateConstraint::exact(2010.0),
      key("D") => DateConstraint::exact(2012.0),
    };
    let node_dates = keymap! {
      key("A") => Some(2010.0),
      key("B") => Some(b),
      key("C") => Some(c),
      key("D") => Some(2012.0),
      key("X") => Some(2005.0),
      key("Y") => Some(2008.0),
      key("root") => Some(2000.0),
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
  fn test_dates_outside_input_dates_skips_nodes_without_dates() -> Result<(), Report> {
    let parsed = nwk_read(TREE.as_bytes())?;
    let names = parsed.names.clone();
    let key = |name: &str| find_node_key_by_name(&parsed.graph, &names, name).expect("node must exist");
    let dates = keymap! {
      key("A") => DateConstraint::exact(2010.0),
      key("C") => DateConstraint::exact(2012.0),
    };
    let node_dates = keymap! {
      key("A") => Some(2010.0),
      key("B") => None,
      key("C") => Some(2012.0),
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

  #[test]
  fn test_dates_outside_input_dates_without_input_dates_gives_no_warning() -> Result<(), Report> {
    let parsed = nwk_read(TREE.as_bytes())?;
    let names = parsed.names.clone();
    let key = |name: &str| find_node_key_by_name(&parsed.graph, &names, name).expect("node must exist");
    let node_dates = keymap! {
      key("A") => Some(2010.0),
      key("B") => Some(2254.4),
      key("C") => Some(1900.0),
      key("D") => Some(2012.0),
      key("X") => Some(2005.0),
      key("Y") => Some(2008.0),
      key("root") => Some(2000.0),
    };

    assert_eq!(
      None,
      dates_outside_input_dates(&parsed.graph, &names, &keymap! {}, &node_dates)
    );
    Ok(())
  }
}

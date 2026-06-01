#[cfg(test)]
mod tests {
  use crate::dates_csv::*;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::datetime::options::DateParserOptions;
  use treetime_utils::o;
  use treetime_utils::pretty_assert_ulps_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::empty("")]
  #[case::whitespace("   ")]
  #[case::nan("NaN")]
  #[case::null("null")]
  #[case::all_x_dash("XXXX-XX-XX")]
  #[case::all_x_slash("XXXX/XX/XX")]
  #[case::all_x_dot("XXXX.XX.XX")]
  #[trace]
  fn test_date_read_empty(#[case] input: &str) -> Result<(), Report>  {
    let actual = read_date(input, &DateParserOptions::default())?;
    assert_eq!(None, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::iso_dash("2024-07-23",   2024.5587431693)]
  #[case::iso_slash("2024/07/23",   2024.5587431693)]
  #[case::iso_dot("2024.07.23",   2024.5587431693)]
  //
  #[case::year_fraction_short("2024.558743",  2024.558743)]
  #[case::year_fraction_long("2024.5587431693",  2024.5587431693)]
  #[case::compact("20240723",     2024.5587431693)]
  #[trace]
  fn test_date_read_exact(#[case] input: &str, #[case] expected: f64) -> Result<(), Report> {
    let constraint = read_date(input, &DateParserOptions::default())?.unwrap();
    assert!(constraint.is_exact(), "expected Exact variant for input {input:?}");
    assert_eq!(input, constraint.raw);
    pretty_assert_ulps_eq!(expected, constraint.mean(), epsilon = 1e-8);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::uncertain_day_dash(  "2024-07-XX",  2024.5396174863)]
  #[case::uncertain_day_slash( "2024/07/XX",  2024.5396174863)]
  #[case::uncertain_day_dot(   "2024.07.XX",  2024.5396174863)]
  #[case::uncertain_month_day( "2024-XX-XX",  2024.5000000000)]
  #[trace]
  fn test_date_read_uncertain(#[case] input: &str, #[case] expected_mean: f64) -> Result<(), Report> {
    let constraint = read_date(input, &DateParserOptions::default())?.unwrap();
    assert!(matches!(constraint.value, DateValue::Uncertain(_)), "expected Uncertain variant for input {input:?}");
    assert_eq!(input, constraint.raw);
    pretty_assert_ulps_eq!(expected_mean, constraint.mean(), epsilon = 1e-8);
    Ok(())
  }

  #[test]
  fn test_date_read_range() -> Result<(), Report> {
    let constraint = read_date("2020-01-01/2020-06-30", &DateParserOptions::default())?.unwrap();
    assert!(
      matches!(constraint.value, DateValue::Range(_)),
      "expected Range variant"
    );
    assert_eq!("2020-01-01/2020-06-30", constraint.raw);
    if let DateValue::Range(r) = constraint.value {
      assert!(r.start < r.end);
      assert!(r.contains(constraint.mean()));
      assert!(r.width() > 0.0);
    }
    Ok(())
  }

  #[test]
  fn test_date_constraint_exact_constructor() {
    let c = DateConstraint::exact(2024.5);
    assert!(c.is_exact());
    #[allow(clippy::float_cmp, reason = "exact literal round-trip, no arithmetic")]
    {
      assert_eq!(2024.5, c.mean());
    }
    assert_eq!("2024.5", c.raw);
  }

  #[test]
  fn test_date_range_contains() {
    let r = DateRange {
      start: 2020.0,
      end: 2021.0,
    };
    assert!(r.contains(2020.0));
    assert!(r.contains(2020.5));
    assert!(r.contains(2021.0));
    assert!(!r.contains(2019.9));
    assert!(!r.contains(2021.1));
  }

  #[test]
  fn test_date_range_width() {
    let r = DateRange {
      start: 2020.0,
      end: 2021.0,
    };
    #[allow(clippy::float_cmp, reason = "exact subtraction of integer-valued literals")]
    {
      assert_eq!(1.0, r.width());
    }
  }

  #[test]
  fn test_read_dates_from_str() -> Result<(), Report> {
    let content = r#"name	 date
A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409	2013.40520192
A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409	2012.83778234
A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409	2009.48186174
A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416	2009.5229295
A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409	2000.13415469
A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416	2000.68172485
A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428	2011.98015058
A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409	2011.95550992
A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409	2012.8569473
A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409	2003.84052019
A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409	2007.48733744
A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409	2008.15058179
A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412	2008.86516085
A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423	2003.00273785
A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409	2012.25735797
A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409	2011.98562628
A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409	2003.00273785
A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416	2011.65160849
A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409	2013.11225188
"#;

    let actual = read_dates_from_str(content, b'\t', &[], &Some(o!("name")), &Some(o!("date")))?;

    let expected = btreemap! {
      o!("A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409") => Some(DateConstraint::exact(2013.40520192)),
      o!("A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409") => Some(DateConstraint::exact(2012.83778234)),
      o!("A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409") => Some(DateConstraint::exact(2009.48186174)),
      o!("A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416") => Some(DateConstraint::exact(2009.5229295)),
      o!("A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409") => Some(DateConstraint::exact(2000.13415469)),
      o!("A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416") => Some(DateConstraint::exact(2000.68172485)),
      o!("A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428") => Some(DateConstraint::exact(2011.98015058)),
      o!("A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409") => Some(DateConstraint::exact(2011.95550992)),
      o!("A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409") => Some(DateConstraint::exact(2012.8569473)),
      o!("A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409") => Some(DateConstraint::exact(2003.84052019)),
      o!("A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409") => Some(DateConstraint::exact(2007.48733744)),
      o!("A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409") => Some(DateConstraint::exact(2008.15058179)),
      o!("A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412") => Some(DateConstraint::exact(2008.86516085)),
      o!("A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423") => Some(DateConstraint::exact(2003.00273785)),
      o!("A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409") => Some(DateConstraint::exact(2012.25735797)),
      o!("A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409") => Some(DateConstraint::exact(2011.98562628)),
      o!("A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409") => Some(DateConstraint::exact(2003.00273785)),
      o!("A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416") => Some(DateConstraint::exact(2011.65160849)),
      o!("A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409") => Some(DateConstraint::exact(2013.11225188)),
    };

    assert_eq!(actual, expected);

    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use crate::csv::delimiter_from_path;
  use crate::dates_csv::*;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::datetime::options::DateParserOptions;
  use treetime_utils::pretty_assert_ulps_eq;
  use treetime_utils::{assert_error, o};

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
  #[case::negative_whole_year(   "-500",      -500.0 )]
  #[case::negative_fraction(     "-500.5",    -500.5 )]
  #[case::year_below_1000(       "500.25",     500.25)]
  #[case::far_future(            "1e9",        1e9   )]
  #[trace]
  fn test_date_read_number_as_written(#[case] input: &str, #[case] expected: f64) -> Result<(), Report> {
    let constraint = read_date(input, &DateParserOptions::default())?.unwrap();
    assert_eq!(DateValue::Exact(DateExact { value: expected }), constraint.value);
    Ok(())
  }

  #[test]
  fn test_date_read_compact_date_before_number() -> Result<(), Report> {
    let constraint = read_date("20200315", &DateParserOptions::default())?.unwrap();
    assert!(constraint.is_exact());
    pretty_assert_ulps_eq!(2020.0 + (75.0 - 0.5) / 366.0, constraint.mean(), max_ulps = 4);
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

  #[rustfmt::skip]
  #[rstest]
  #[case::built_in_reads_day_first(   "03/06/2016",  None,              Some("2016-06-03"))]
  #[case::format_reads_month_first(   "03/06/2016",  Some("%m/%d/%Y"),  Some("2016-03-06"))]
  #[case::built_in_lacks_month_name(  "13 Jun 2016", None,              None)]
  #[case::format_reads_month_name(    "13 Jun 2016", Some("%d %b %Y"),  Some("2016-06-13"))]
  #[case::format_keeps_iso(           "2016-06-13",  Some("%d %b %Y"),  Some("2016-06-13"))]
  #[case::format_keeps_year_fraction( "2016.45",     Some("%d %b %Y"),  Some("2016.45"))]
  #[case::format_keeps_uncertain(     "2016-06-XX",  Some("%d %b %Y"),  Some("2016-06-XX"))]
  #[trace]
  fn test_date_read_with_date_format(
    #[case] input: &str,
    #[case] date_format: Option<&str>,
    #[case] expected_as: Option<&str>,
  ) -> Result<(), Report> {
    let options = DateParserOptions {
      date_format: date_format.map(str::to_owned),
      ..DateParserOptions::default()
    };
    let expected = expected_as
      .map(|reference| read_date(reference, &DateParserOptions::default()))
      .transpose()?
      .flatten()
      .map(|constraint| constraint.value);
    let actual = read_date(input, &options)?.map(|constraint| constraint.value);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_date_read_table_dates_with_date_format() -> Result<(), Report> {
    let table = metadata_read(
      &b"name\tdate\nA\t03/06/2016\nB\t13 Jun 2016\n"[..],
      delimiter_from_path("metadata.tsv"),
      &['\t'],
      &[],
      Some("name"),
      Some("date"),
    )?;
    let options = DateParserOptions {
      date_format: Some(o!("%m/%d/%Y")),
      ..DateParserOptions::default()
    };
    let expected = vec![
      (
        o!("A"),
        read_date("2016-03-06", &DateParserOptions::default())?.map(|date| date.value),
      ),
      (o!("B"), None),
    ];
    let actual = table
      .dates(&options)?
      .into_iter()
      .map(|(name, date)| (name, date.map(|date| date.value)))
      .collect::<Vec<_>>();
    assert_eq!(expected, actual);
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
    #[expect(clippy::float_cmp, reason = "exact literal round-trip, no arithmetic")]
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
    #[expect(clippy::float_cmp, reason = "exact subtraction of integer-valued literals")]
    {
      assert_eq!(1.0, r.width());
    }
  }

  #[test]
  fn test_metadata_read_dates() -> Result<(), Report> {
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

    let actual = metadata_read(
      content.as_bytes(),
      delimiter_from_path("dates.tsv"),
      &['\t'],
      &[],
      Some("name"),
      Some("date"),
    )?
    .dates(&DateParserOptions::default())?;

    let expected = vec![
      (
        o!("A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409"),
        Some(DateConstraint::exact(2013.40520192)),
      ),
      (
        o!("A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409"),
        Some(DateConstraint::exact(2012.83778234)),
      ),
      (
        o!("A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409"),
        Some(DateConstraint::exact(2009.48186174)),
      ),
      (
        o!("A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416"),
        Some(DateConstraint::exact(2009.5229295)),
      ),
      (
        o!("A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409"),
        Some(DateConstraint::exact(2000.13415469)),
      ),
      (
        o!("A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416"),
        Some(DateConstraint::exact(2000.68172485)),
      ),
      (
        o!("A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428"),
        Some(DateConstraint::exact(2011.98015058)),
      ),
      (
        o!("A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409"),
        Some(DateConstraint::exact(2011.95550992)),
      ),
      (
        o!("A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409"),
        Some(DateConstraint::exact(2012.8569473)),
      ),
      (
        o!("A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409"),
        Some(DateConstraint::exact(2003.84052019)),
      ),
      (
        o!("A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409"),
        Some(DateConstraint::exact(2007.48733744)),
      ),
      (
        o!("A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409"),
        Some(DateConstraint::exact(2008.15058179)),
      ),
      (
        o!("A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412"),
        Some(DateConstraint::exact(2008.86516085)),
      ),
      (
        o!("A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423"),
        Some(DateConstraint::exact(2003.00273785)),
      ),
      (
        o!("A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409"),
        Some(DateConstraint::exact(2012.25735797)),
      ),
      (
        o!("A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409"),
        Some(DateConstraint::exact(2011.98562628)),
      ),
      (
        o!("A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409"),
        Some(DateConstraint::exact(2003.00273785)),
      ),
      (
        o!("A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416"),
        Some(DateConstraint::exact(2011.65160849)),
      ),
      (
        o!("A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409"),
        Some(DateConstraint::exact(2013.11225188)),
      ),
    ];

    assert_eq!(actual, expected);

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::csv(                 "metadata.csv",    "strain,date,country\nA,2020-01-15,usa\nB,2020-02-03,peru\n",          None,           (',',  Some("date")))]
  #[case::tsv(                 "metadata.tsv",    "strain\tdate\tcountry\nA\t2020-01-15\tusa\nB\t2020-02-03\tperu\n",    None,           ('\t', Some("date")))]
  #[case::tsv_named_as_csv(    "metadata.csv",    "strain\tdate\tcountry\nA\t2020-01-15\tusa\nB\t2020-02-03\tperu\n",    None,           ('\t', Some("date")))]
  #[case::csv_named_as_tsv(    "metadata.tsv.xz", "strain,date,country\nA,2020-01-15,usa\nB,2020-02-03,peru\n",          None,           (',',  Some("date")))]
  #[case::without_dates(       "metadata.tsv",    "strain\tcountry\nA\tusa\nB\tperu\n",                                  None,           ('\t', None))]
  #[case::explicit_date_column("metadata.tsv",    "strain\tdate,x\nA\t2020,1\nB\t2021,2\n",                             Some("date,x"), ('\t', Some("date,x")))]
  #[case::tie_broken_by_tsv(   "metadata.tsv",    "strain,x\tstrain\tdate\nq\tA\t2020-01-15\nq\tB\t2020-02-03\n",       None,           ('\t', Some("date")))]
  #[case::tie_broken_by_csv(   "metadata.csv",    "strain,date,x\tstrain\nA,2020-01-15,q\nB,2020-02-03,q\n",              None,           (',',  Some("date")))]
  #[case::num_date(            "metadata.tsv",    "strain\tnum_date\nA\t2020.1\nB\t2020.2\n",                            None,           ('\t', Some("num_date")))]
  #[case::leftmost_date_like(  "metadata.tsv",    "strain\tcollection_date\tdate\nA\t2020\t2021\nB\t2020\t2021\n",       None,           ('\t', Some("collection_date")))]
  #[case::explicit_hash_date(  "metadata.tsv",    "strain\t#date\nA\t2020\nB\t2021\n",                                  Some("#date"),  ('\t', Some("date")))]
  #[trace]
  fn test_metadata_read_picks_columns(
    #[case] file_name: &str,
    #[case] content: &str,
    #[case] date_column: Option<&str>,
    #[case] (delimiter, expected_date_column): (char, Option<&str>),
  ) {
    let table = metadata_read(
      content.as_bytes(),
      delimiter_from_path(file_name),
      &[',', '\t', ';'],
      &[o!("strain"), o!("name")],
      None,
      date_column,
    )
    .unwrap();
    assert_eq!(
      (delimiter, o!("strain"), expected_date_column.map(str::to_owned), vec![o!("A"), o!("B")]),
      (
        table.delimiter,
        table.id_column.clone(),
        table.date_column.as_ref().ok().cloned(),
        table.rows.iter().map(|row| row.name.clone()).collect::<Vec<_>>(),
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::candidate(vec![o!("#name")], None)]
  #[case::provided( vec![],            Some("#name"))]
  #[trace]
  fn test_metadata_read_hash_name_column(#[case] name_candidates: Vec<String>, #[case] name_column: Option<&str>) {
    let table = metadata_read(
      &b"#name\tdate\tcountry\nA\t2020\tusa\n"[..],
      delimiter_from_path("metadata.tsv"),
      &[',', '\t'],
      &name_candidates,
      name_column,
      None,
    )
    .unwrap();
    assert_eq!(
      (o!("name"), vec![o!("A")]),
      (table.id_column.clone(), table.rows.iter().map(|row| row.name.clone()).collect::<Vec<_>>())
    );
  }

  #[test]
  fn test_metadata_read_explicit_date_column_as_written() {
    let table = metadata_read(
      &b"name\tdate\t#date\nA\t2020\t2021\n"[..],
      delimiter_from_path("metadata.tsv"),
      &['\t'],
      &[o!("name")],
      None,
      Some("#date"),
    )
    .unwrap();
    assert_eq!(
      vec![Some(o!("2021"))],
      table.rows.iter().map(|row| row.date.clone()).collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_metadata_read_explicit_date_column_is_case_sensitive() {
    let table = metadata_read(
      &b"strain\tdate\nA\t2020\n"[..],
      delimiter_from_path("metadata.tsv"),
      &['\t'],
      &[o!("strain")],
      None,
      Some("Date"),
    )
    .unwrap();
    assert_error!(
      table.dates(&DateParserOptions::default()),
      "Unable to find column 'Date'. Available columns are: strain, date"
    );
  }

  #[test]
  fn test_metadata_read_dates_require_a_date_column() {
    let table = metadata_read(
      &b"strain\tcountry\nA\tusa\n"[..],
      delimiter_from_path("metadata.tsv"),
      &[',', '\t'],
      &[o!("strain")],
      None,
      None,
    )
    .unwrap();
    assert_error!(
      table.dates(&DateParserOptions::default()),
      "Unable to find column:\n  Looking for: a name that contains 'date'\n  Available columns are: strain, country"
    );
  }
}

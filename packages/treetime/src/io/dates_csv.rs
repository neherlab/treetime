use crate::io::csv::{get_col_name, guess_csv_delimiter};
use crate::{make_internal_report, vec_of_owned};
use csv::{ReaderBuilder as CsvReaderBuilder, StringRecord, Trim};
use eyre::{Report, WrapErr, eyre};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::io::Read;
use std::path::Path;
use treetime_utils::datetime::options::DateParserOptions;
use treetime_utils::datetime::parse_date::{parse_date, parse_date_range};
use treetime_utils::datetime::parse_uncertain_date::parse_date_uncertain;
use treetime_utils::datetime::year_frac::{date_range_to_year_fraction_range, date_to_year_fraction};
use treetime_utils::file::open_file_or_stdin;

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum DateOrRange {
  YearFraction(f64),
  YearFractionRange((f64, f64)),
}

impl DateOrRange {
  #[inline]
  pub fn mean(&self) -> f64 {
    match self {
      DateOrRange::YearFraction(date) => *date,
      DateOrRange::YearFractionRange((begin, end)) => (begin + end) / 2.0,
    }
  }
}

pub type DatesMap = BTreeMap<String, Option<DateOrRange>>;
pub type DateRecord = (String, Option<DateOrRange>);

pub fn read_dates_from_reader(
  reader: impl Read,
  delimiter: u8,
  name_column: &Option<String>,
  date_column: &Option<String>,
) -> Result<DatesMap, Report> {
  let mut reader = CsvReaderBuilder::new()
    .trim(Trim::All)
    .delimiter(delimiter)
    .from_reader(reader);

  let headers = reader
    .headers()
    .map(|record| record.iter().map(str::to_owned).collect_vec())
    .map_err(|err| eyre!("{err}"))?
    .iter()
    .map(|header| header.trim_start_matches('#').trim_end_matches('#').trim().to_owned())
    .collect_vec();

  let name_column_idx = get_col_name(&headers, &vec_of_owned!["name", "strain", "accession"], name_column)?;
  let date_column_idx = get_col_name(&headers, &vec_of_owned!["date"], date_column)?;

  reader
    .records()
    .enumerate()
    .map(|(index, record)| {
      let record = record?;
      convert_record(index, &record, name_column_idx, date_column_idx)
        .wrap_err_with(|| format!("When reading row {index}, column '{date_column_idx}'"))
    })
    .collect::<Result<DatesMap, Report>>()
}

pub fn read_dates_from_str(
  content: &str,
  delimiter: u8,
  name_column: &Option<String>,
  date_column: &Option<String>,
) -> Result<DatesMap, Report> {
  let reader = content.as_bytes();
  read_dates_from_reader(reader, delimiter, name_column, date_column)
}

pub fn read_dates(
  filepath: impl AsRef<Path>,
  name_column: &Option<String>,
  date_column: &Option<String>,
) -> Result<DatesMap, Report> {
  let filepath = filepath.as_ref();
  let file =
    open_file_or_stdin(&Some(filepath)).wrap_err_with(|| format!("When reading file: '{}'", filepath.display()))?;
  let delimiter = guess_csv_delimiter(filepath)
    .wrap_err_with(|| format!("When guessing CSV delimiter for '{}'", filepath.display()))?;
  read_dates_from_reader(file, delimiter, name_column, date_column)
    .wrap_err_with(|| format!("When reading dates from file: '{}'", filepath.display()))
}

pub fn convert_record(
  index: usize,
  record: &StringRecord,
  name_column_idx: usize,
  date_column_idx: usize,
) -> Result<DateRecord, Report> {
  let name = record
    .get(name_column_idx)
    .ok_or_else(|| make_internal_report!("Row '{index}': Unable to get column with index '{name_column_idx}'"))?
    .to_owned();

  let date = record
    .get(date_column_idx)
    .ok_or_else(|| make_internal_report!("Row '{index}': Unable to get column with index '{date_column_idx}'"))?;

  let date = read_date(date, &DateParserOptions::default())?;

  Ok((name, date))
}

pub fn read_date(date_str: &str, options: &DateParserOptions) -> Result<Option<DateOrRange>, Report> {
  let trimmed = date_str.trim();

  // Try parsing as fractional year first (e.g. "2020.5", "2003.84052019")
  if let Ok(year_fraction) = trimmed.parse::<f64>() {
    if year_fraction.is_finite() {
      return Ok(Some(DateOrRange::YearFraction(year_fraction)));
    }
  }

  // Try parsing as formatted date (e.g. "2020-07-15", "20200715", "2020/07/15")
  if let Ok(date) = parse_date(trimmed, options) {
    return Ok(Some(DateOrRange::YearFraction(date_to_year_fraction(&date))));
  }

  // Try parsing as uncertain date (e.g. "2020-XX-XX")
  if let Ok(date_range) = parse_date_uncertain(trimmed, options) {
    return Ok(Some(DateOrRange::YearFractionRange(date_range_to_year_fraction_range(
      &date_range,
    ))));
  }

  // Try parsing as date range (e.g. "2020-01-01/2020-12-31")
  if let Ok(date_range) = parse_date_range(trimmed, options) {
    return Ok(Some(DateOrRange::YearFractionRange(date_range_to_year_fraction_range(
      &date_range,
    ))));
  }

  // Unable to parse
  Ok(None)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::o;

  #[rustfmt::skip]
  #[rstest]
  #[case("")]
  #[case("   ")]
  #[case("NaN")]
  #[case("null")]
  #[case("XXXX-XX-XX")]
  #[case("XXXX/XX/XX")]
  #[case("XXXX.XX.XX")]
  #[trace]
  fn test_date_read_empty(#[case] input: &str) -> Result<(), Report>  {
    let actual = read_date(input, &DateParserOptions::default())?;
    assert_eq!(None, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case("2024-07-23",   2024.5587431693)]
  #[case("2024/07/23",   2024.5587431693)]
  #[case("2024.07.23",   2024.5587431693)]
  //
  #[case("2024-07-XX",   2024.5396174863)]
  #[case("2024/07/XX",   2024.5396174863)]
  #[case("2024.07.XX",   2024.5396174863)]
  //
  #[case("2024-XX-XX",   2024.5000000000)]
  #[case("2024/XX/XX",   2024.5000000000)]
  #[case("2024.XX.XX",   2024.5000000000)]
  //
  #[case("2024.558743",  2024.558743)]
  #[case("2024.5587431693",  2024.5587431693)]
  #[case("20240723",     2024.5587431693)]
  #[trace]
  fn test_date_read_ok(#[case] input: &str, #[case] expected: f64) -> Result<(), Report>  {
    let actual = read_date(input, &DateParserOptions::default())?.unwrap().mean();
    pretty_assert_ulps_eq!(expected, actual, epsilon = 1e-8);
    Ok(())
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

    let actual = read_dates_from_str(content, b'\t', &Some(o!("name")), &Some(o!("date")))?;

    let expected = btreemap! {
      o!("A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2013.40520192)),
      o!("A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2012.83778234)),
      o!("A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409") => Some(DateOrRange::YearFraction(2009.48186174)),
      o!("A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416") => Some(DateOrRange::YearFraction(2009.5229295)),
      o!("A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409") => Some(DateOrRange::YearFraction(2000.13415469)),
      o!("A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416") => Some(DateOrRange::YearFraction(2000.68172485)),
      o!("A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428") => Some(DateOrRange::YearFraction(2011.98015058)),
      o!("A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409") => Some(DateOrRange::YearFraction(2011.95550992)),
      o!("A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2012.8569473)),
      o!("A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409") => Some(DateOrRange::YearFraction(2003.84052019)),
      o!("A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409") => Some(DateOrRange::YearFraction(2007.48733744)),
      o!("A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409") => Some(DateOrRange::YearFraction(2008.15058179)),
      o!("A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412") => Some(DateOrRange::YearFraction(2008.86516085)),
      o!("A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423") => Some(DateOrRange::YearFraction(2003.00273785)),
      o!("A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409") => Some(DateOrRange::YearFraction(2012.25735797)),
      o!("A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409") => Some(DateOrRange::YearFraction(2011.98562628)),
      o!("A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409") => Some(DateOrRange::YearFraction(2003.00273785)),
      o!("A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416") => Some(DateOrRange::YearFraction(2011.65160849)),
      o!("A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2013.11225188)),
    };

    assert_eq!(actual, expected);

    Ok(())
  }
}

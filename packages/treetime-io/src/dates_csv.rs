use crate::csv::{delimiter_from_path, detect_csv_delimiter, get_col_name, normalize_csv_headers, table_reader};
use eyre::{Report, WrapErr};
use std::io::BufRead;
use std::path::Path;
pub use treetime_primitives::date::{DateConstraint, DateExact, DateRange, DateValue, DatesMap};
use treetime_utils::datetime::options::DateParserOptions;
use treetime_utils::datetime::parse_date::{parse_date, parse_date_range};
use treetime_utils::datetime::parse_uncertain_date::parse_date_uncertain;
use treetime_utils::datetime::year_fraction::{date_range_to_year_fraction_range, date_to_year_fraction};
use treetime_utils::io::file::read_file_with;
use treetime_utils::{make_internal_report, make_report, vec_of_owned};

#[derive(Debug)]
pub struct MetadataTable {
  pub delimiter: char,
  pub columns: Vec<String>,
  pub id_column: String,
  pub date_column: Result<String, Report>,
  pub rows: Vec<MetadataRow>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MetadataRow {
  pub name: String,
  pub date: Option<String>,
}

impl MetadataTable {
  pub fn dates(&self) -> Result<DatesMap, Report> {
    let column = self.date_column.as_ref().map_err(|report| make_report!("{report}"))?;
    let options = DateParserOptions::default();
    self
      .rows
      .iter()
      .enumerate()
      .map(|(index, row)| {
        let raw = row
          .date
          .as_deref()
          .ok_or_else(|| make_internal_report!("Row '{index}' has no value in date column '{column}'"))?;
        let date = read_date(raw, &options).wrap_err_with(|| format!("When reading row {index}, column '{column}'"))?;
        Ok((row.name.clone(), date))
      })
      .collect()
  }
}

pub fn metadata_read_file(
  filepath: impl AsRef<Path>,
  delimiters: &[char],
  name_candidates: &[String],
  name_column: Option<&str>,
  date_column: Option<&str>,
) -> Result<MetadataTable, Report> {
  let filepath = filepath.as_ref();
  read_file_with(filepath, |reader| {
    metadata_read(
      reader,
      delimiter_from_path(filepath),
      delimiters,
      name_candidates,
      name_column,
      date_column,
    )
  })
}

pub fn metadata_read(
  mut reader: impl BufRead,
  path_delimiter: Option<u8>,
  delimiters: &[char],
  name_candidates: &[String],
  name_column: Option<&str>,
  date_column: Option<&str>,
) -> Result<MetadataTable, Report> {
  let date_candidates = vec_of_owned!["date"];
  let delimiter = detect_csv_delimiter(&mut reader, path_delimiter, delimiters, |headers| {
    get_col_name(headers, name_candidates, name_column).is_ok()
      && get_col_name(headers, &date_candidates, date_column).is_ok()
  })
  .or_else(|_without_dates| {
    detect_csv_delimiter(&mut reader, path_delimiter, delimiters, |headers| {
      get_col_name(headers, name_candidates, name_column).is_ok()
    })
  })
  .wrap_err("When detecting the metadata delimiter")?;

  let mut reader = table_reader(reader, delimiter);
  let columns = reader
    .headers()
    .map(normalize_csv_headers)
    .map_err(|err| make_report!("{err}"))?;
  let id_index = get_col_name(&columns, name_candidates, name_column)?;
  let date_index = get_col_name(&columns, &date_candidates, date_column);

  let rows = reader
    .records()
    .enumerate()
    .map(|(index, record)| {
      let record = record?;
      let name = record
        .get(id_index)
        .ok_or_else(|| make_internal_report!("Row '{index}': Unable to get column with index '{id_index}'"))?
        .to_owned();
      let date = date_index
        .as_ref()
        .ok()
        .map(|&date_index| {
          record
            .get(date_index)
            .map(str::to_owned)
            .ok_or_else(|| make_internal_report!("Row '{index}': Unable to get column with index '{date_index}'"))
        })
        .transpose()?;
      Ok(MetadataRow { name, date })
    })
    .collect::<Result<Vec<_>, Report>>()?;

  Ok(MetadataTable {
    delimiter: char::from(delimiter),
    id_column: columns[id_index].clone(),
    date_column: date_index.map(|index| columns[index].clone()),
    columns,
    rows,
  })
}

#[cfg_attr(
  dylint_lib = "treetime_lints",
  expect(
    error_dropped_by_pattern,
    reason = "each date notation is tried in turn; a value that matches none is reported by the caller as missing"
  )
)]
pub(crate) fn read_date(date_str: &str, options: &DateParserOptions) -> Result<Option<DateConstraint>, Report> {
  let trimmed = date_str.trim();

  if let Ok(year_fraction) = trimmed.parse::<f64>()
    && year_fraction.is_finite()
    && (1000.0..=3000.0).contains(&year_fraction)
  {
    return Ok(Some(DateConstraint {
      raw: trimmed.to_owned(),
      value: DateValue::Exact(DateExact { value: year_fraction }),
    }));
  }

  if let Ok(date) = parse_date(trimmed, options) {
    return Ok(Some(DateConstraint {
      raw: trimmed.to_owned(),
      value: DateValue::Exact(DateExact {
        value: date_to_year_fraction(&date),
      }),
    }));
  }

  if let Ok(date_range) = parse_date_uncertain(trimmed, options) {
    let (start, end) = date_range_to_year_fraction_range(&date_range);
    return Ok(Some(DateConstraint {
      raw: trimmed.to_owned(),
      value: DateValue::Uncertain(DateRange { start, end }),
    }));
  }

  if let Ok(date_range) = parse_date_range(trimmed, options) {
    let (start, end) = date_range_to_year_fraction_range(&date_range);
    return Ok(Some(DateConstraint {
      raw: trimmed.to_owned(),
      value: DateValue::Range(DateRange { start, end }),
    }));
  }

  Ok(None)
}

use deser::{Deserialize, Serialize};
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::io::BufRead;
use std::ops::RangeInclusive;
use std::path::Path;
use treetime_utils::io::file::read_file_with;
use treetime_utils::make_error;

const LATITUDES: RangeInclusive<f64> = -90.0..=90.0;

const LONGITUDES: RangeInclusive<f64> = -180.0..=180.0;

pub fn lat_longs_read_file(filepath: impl AsRef<Path>) -> Result<LatLongs, Report> {
  read_file_with(filepath, lat_longs_read)
}

pub fn lat_longs_read(reader: impl BufRead) -> Result<LatLongs, Report> {
  let mut locations = BTreeMap::new();
  for (index, line) in reader.lines().enumerate() {
    let number = index + 1;
    let line = line.wrap_err_with(|| format!("When reading line {number}"))?;
    if let Some((key, coordinates)) = location(&line).wrap_err_with(|| format!("Line {number}"))? {
      locations.insert(key, coordinates);
    }
  }
  Ok(LatLongs { locations })
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct LatLongs {
  locations: BTreeMap<(String, String), Coordinates>,
}

impl LatLongs {
  pub fn get(&self, resolution: &str, name: &str) -> Option<Coordinates> {
    self
      .locations
      .get(&(resolution.to_lowercase(), name.to_lowercase()))
      .copied()
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Coordinates {
  pub latitude: f64,
  pub longitude: f64,
}

fn location(line: &str) -> Result<Option<((String, String), Coordinates)>, Report> {
  let trimmed = line.trim();
  if trimmed.is_empty() || trimmed.starts_with('#') {
    return Ok(None);
  }
  let fields: Vec<&str> = if line.contains('\t') {
    trimmed.split('\t').collect()
  } else {
    trimmed.split_whitespace().collect()
  };
  let [resolution, name, latitude, longitude] = fields[..] else {
    return make_error!(
      "expected 4 fields (resolution, name, latitude and longitude), found {}: '{trimmed}'",
      fields.len()
    );
  };
  let coordinates = Coordinates {
    latitude: degrees("latitude", latitude, &LATITUDES)?,
    longitude: degrees("longitude", longitude, &LONGITUDES)?,
  };
  Ok(Some(((resolution.to_lowercase(), name.to_lowercase()), coordinates)))
}

fn degrees(what: &str, text: &str, range: &RangeInclusive<f64>) -> Result<f64, Report> {
  let value = text
    .trim()
    .parse::<f64>()
    .map_err(Report::from)
    .wrap_err_with(|| format!("{what} '{text}' is not a number"))?;
  if !range.contains(&value) {
    return make_error!(
      "{what} {text} is outside the range {} to {}",
      range.start(),
      range.end()
    );
  }
  Ok(value)
}

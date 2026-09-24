use csv::{ReaderBuilder, Trim, Writer, WriterBuilder};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use serde::Serialize;
use std::io::{BufRead, Write};
use std::path::Path;
use treetime_utils::io::compression::remove_compression_ext;
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::io::fs::extension;
use treetime_utils::make_error;
use treetime_utils::make_report;

pub struct CsvStructFileWriter {
  writer: CsvStructWriter<Box<dyn Write + Send>>,
}

impl CsvStructFileWriter {
  pub fn new(filepath: impl AsRef<Path>, delimiter: u8) -> Result<Self, Report> {
    let filepath = filepath.as_ref();
    let file = create_file_or_stdout(filepath)?;
    let writer = CsvStructWriter::new(file, delimiter)?;
    Ok(Self { writer })
  }

  pub fn write<T: Serialize>(&mut self, record: &T) -> Result<(), Report> {
    self.writer.write(record)?;
    Ok(())
  }
}

pub struct CsvStructWriter<W: Write + Send> {
  pub writer: Writer<W>,
}

impl<W: Write + Send> CsvStructWriter<W> {
  pub fn new(writer: W, delimiter: u8) -> Result<Self, Report> {
    let writer = WriterBuilder::new().delimiter(delimiter).from_writer(writer);
    Ok(Self { writer })
  }

  pub fn write<T: Serialize>(&mut self, record: &T) -> Result<(), Report> {
    self.writer.serialize(record)?;
    Ok(())
  }
}

pub fn default_name_candidates() -> Vec<String> {
  vec!["strain".to_owned(), "name".to_owned(), "accession".to_owned()]
}

pub fn default_metadata_delimiters() -> Vec<char> {
  vec![',', '\t', ';']
}

pub(crate) fn get_col_name(
  headers: &[String],
  possible_names: &[String],
  provided_name: Option<&str>,
) -> Result<usize, Report> {
  if let Some(provided_name) = provided_name {
    match headers.iter().position(|header| header == provided_name) {
      Some(idx) => Ok(idx),
      None => make_error!(
        "Unable to find column '{provided_name}'. Available columns are: {}",
        headers.join(", ")
      ),
    }
  } else {
    let candidates_lower: Vec<String> = possible_names.iter().map(|c| c.to_lowercase()).collect();
    headers
      .iter()
      .enumerate()
      .find_map(|(idx, header)| {
        let header_lower = header.to_lowercase();
        candidates_lower.contains(&header_lower).then_some(idx)
      })
      .ok_or_else(|| {
        make_report!(
          "Unable to find column:\n  Looking for: {}\n  Available columns are: {}",
          possible_names.join(", "),
          headers.join(", ")
        )
      })
  }
}

pub(crate) fn detect_csv_delimiter<R: BufRead + ?Sized>(
  reader: &mut R,
  filepath: impl AsRef<Path>,
  delimiters: &[char],
  header_matches: impl Fn(&[String]) -> bool,
) -> Result<u8, Report> {
  const SAMPLE_SIZE: usize = 64 * 1024;

  let delimiters = delimiters
    .iter()
    .copied()
    .unique()
    .map(delimiter_to_byte)
    .collect::<Result<Vec<_>, _>>()?;

  if delimiters.is_empty() {
    return make_error!("At least one metadata delimiter is required");
  }
  if let [delimiter] = delimiters.as_slice() {
    return Ok(*delimiter);
  }

  let filepath = filepath.as_ref();
  let sample = reader.fill_buf().wrap_err_with(|| {
    format!(
      "When reading the start of '{}' to detect its delimiter",
      filepath.display()
    )
  })?;
  let sample = &sample[..sample.len().min(SAMPLE_SIZE)];
  let matches = delimiters
    .iter()
    .copied()
    .filter(|delimiter| csv_headers(sample, *delimiter).is_ok_and(|headers| header_matches(&headers)))
    .collect_vec();

  if let [delimiter] = matches.as_slice() {
    Ok(*delimiter)
  } else {
    let path_delimiter = delimiter_from_path(filepath);
    if let Some(delimiter) = path_delimiter
      && delimiters.contains(&delimiter)
      && (matches.is_empty() || matches.contains(&delimiter))
    {
      return Ok(delimiter);
    }
    make_error!(
      "Unable to detect metadata delimiter from candidates: {}",
      delimiters
        .iter()
        .map(|delimiter| format!("{:?}", char::from(*delimiter)))
        .join(", ")
    )
  }
}

fn delimiter_to_byte(delimiter: char) -> Result<u8, Report> {
  u8::try_from(u32::from(delimiter))
    .map_err(Report::from)
    .wrap_err_with(|| format!("Metadata delimiter {delimiter:?} must fit in one byte"))
}

fn csv_headers(sample: &[u8], delimiter: u8) -> Result<Vec<String>, csv::Error> {
  let mut reader = ReaderBuilder::new()
    .trim(Trim::All)
    .delimiter(delimiter)
    .from_reader(sample);
  reader.headers().map(normalize_csv_headers)
}

pub(crate) fn normalize_csv_headers(headers: &csv::StringRecord) -> Vec<String> {
  headers
    .iter()
    .map(|header| header.trim_start_matches('#').trim_end_matches('#').trim().to_owned())
    .collect()
}

fn delimiter_from_path(filepath: impl AsRef<Path>) -> Option<u8> {
  let filepath = remove_compression_ext(filepath);
  match extension(filepath)?.to_lowercase().as_str() {
    "csv" => Some(b','),
    "tsv" => Some(b'\t'),
    "ssv" => Some(b';'),
    _ => None,
  }
}

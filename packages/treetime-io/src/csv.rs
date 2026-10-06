use deser::Serialize;
use deser::de::DeserializeOwned;
use deser::io::{Reader, Writer};
use deser_csv::{DeserializerConfig, Headers, Serializer, SerializerConfig, StreamDeserializer, Trim};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::io::{BufRead, Read, Write};
use std::path::{Path, PathBuf};
use treetime_utils::io::compression::remove_compression_ext;
use treetime_utils::io::file::{FileWriter, create_file_or_stdout, read_file_with, write_file_with};
use treetime_utils::io::fs::extension;
use treetime_utils::make_error;
use treetime_utils::make_report;

pub const DELIMITED_EXTENSIONS: [(&str, u8); 3] = [("csv", b','), ("tsv", b'\t'), ("ssv", b';')];

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TableFormat {
  Csv,
  Tsv,
}

impl TableFormat {
  pub const fn delimiter(self) -> u8 {
    match self {
      Self::Csv => b',',
      Self::Tsv => b'\t',
    }
  }
}

pub fn csv_read_file<T: DeserializeOwned>(filepath: impl AsRef<Path>, format: TableFormat) -> Result<Vec<T>, Report> {
  read_file_with(filepath, |reader| csv_read(reader, format))
}

pub fn csv_read<T: DeserializeOwned>(reader: impl Read, format: TableFormat) -> Result<Vec<T>, Report> {
  table_config(format.delimiter())
    .reader(reader)
    .iter()
    .enumerate()
    .map(|(index, record)| record.wrap_err_with(|| format!("When parsing row {}", index + 1)))
    .collect()
}

pub fn csv_write_file<T: Serialize>(
  filepath: impl AsRef<Path>,
  rows: impl IntoIterator<Item = T>,
  format: TableFormat,
) -> Result<(), Report> {
  write_file_with(filepath, |writer| {
    let mut csv = CsvWriter::new(writer, format);
    rows.into_iter().try_for_each(|row| csv.write_row(&row))?;
    csv.into_inner()?;
    Ok(())
  })
}

pub struct CsvWriter<W: Write> {
  writer: Writer<W, Serializer>,
  filepath: Option<PathBuf>,
}

impl<W: Write> CsvWriter<W> {
  pub fn new(writer: W, format: TableFormat) -> Self {
    Self::with_filepath(writer, format, None)
  }

  pub fn write_row<T: Serialize>(&mut self, row: &T) -> Result<(), Report> {
    let result = self.writer.write(row).wrap_err("When writing a table row");
    self.wrap_filepath(result)
  }

  pub fn write_record<I, T>(&mut self, record: I) -> Result<(), Report>
  where
    I: IntoIterator<Item = T>,
    T: AsRef<str>,
  {
    let record = record.into_iter().collect_vec();
    let fields = record.iter().map(AsRef::as_ref).collect_vec();
    let result = self.writer.write(&fields).wrap_err("When writing a table row");
    self.wrap_filepath(result)
  }

  pub fn into_inner(self) -> Result<W, Report> {
    let Self { mut writer, filepath } = self;
    let result = writer
      .flush()
      .map(|()| writer.into_inner())
      .wrap_err("When flushing the table");
    match filepath {
      Some(filepath) => result.wrap_err_with(|| format!("When writing file '{}'", filepath.display())),
      None => result,
    }
  }

  fn with_filepath(writer: W, format: TableFormat, filepath: Option<PathBuf>) -> Self {
    Self {
      writer: SerializerConfig::new().delimiter(format.delimiter()).writer(writer),
      filepath,
    }
  }

  fn wrap_filepath<T>(&self, result: Result<T, Report>) -> Result<T, Report> {
    match &self.filepath {
      Some(filepath) => result.wrap_err_with(|| format!("When writing file '{}'", filepath.display())),
      None => result,
    }
  }
}

impl CsvWriter<FileWriter> {
  pub fn create(filepath: impl AsRef<Path>, format: TableFormat) -> Result<Self, Report> {
    let filepath = filepath.as_ref();
    Ok(Self::with_filepath(
      create_file_or_stdout(filepath)?,
      format,
      Some(filepath.to_owned()),
    ))
  }

  pub fn finish(self) -> Result<(), Report> {
    self.into_inner()?.finish()
  }
}

pub(crate) fn table_records<R: Read>(reader: R, delimiter: u8) -> Reader<R, StreamDeserializer> {
  DeserializerConfig::new()
    .delimiter(delimiter)
    .trim(Trim::All)
    .headers(Headers::None)
    .reader(reader)
}

const fn table_config(delimiter: u8) -> DeserializerConfig {
  DeserializerConfig::new().delimiter(delimiter).trim(Trim::All)
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
  path_delimiter: Option<u8>,
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

  let sample = reader
    .fill_buf()
    .wrap_err("When reading the start of the table to detect its delimiter")?;
  let sample = &sample[..sample.len().min(SAMPLE_SIZE)];
  let matches = delimiters
    .iter()
    .copied()
    .filter(|delimiter| csv_headers(sample, *delimiter).is_ok_and(|headers| header_matches(&headers)))
    .collect_vec();

  if let [delimiter] = matches.as_slice() {
    Ok(*delimiter)
  } else {
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

fn csv_headers(sample: &[u8], delimiter: u8) -> Result<Vec<String>, Report> {
  table_records(sample, delimiter)
    .read::<Vec<String>>()?
    .map(|headers| normalize_csv_headers(&headers))
    .ok_or_else(|| make_report!("The table is empty"))
}

pub(crate) fn normalize_csv_headers(headers: &[String]) -> Vec<String> {
  headers
    .iter()
    .map(|header| header.trim_start_matches('#').trim_end_matches('#').trim().to_owned())
    .collect()
}

pub(crate) fn delimiter_from_path(filepath: impl AsRef<Path>) -> Option<u8> {
  let filepath = remove_compression_ext(filepath);
  let ext = extension(filepath)?.to_lowercase();
  DELIMITED_EXTENSIONS
    .iter()
    .find(|(known, _)| *known == ext)
    .map(|(_, delimiter)| *delimiter)
}

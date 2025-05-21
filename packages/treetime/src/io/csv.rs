use crate::io::file::create_file_or_stdout;
use crate::io::fs::{extension, read_file_to_string};
use crate::make_error;
use crate::utils::error::to_eyre_error;
use csv::{ReaderBuilder as CsvReaderBuilder, Writer as CsvWriterImpl, WriterBuilder as CsvWriterBuilder};
use eyre::{Report, eyre};
use serde::{Deserialize, Serialize};
use std::io::Write;
use std::path::{Path, PathBuf};

/// Writes CSV. Each row is a serde-annotated struct.
pub struct CsvStructWriter<W: Write + Send> {
  pub writer: CsvWriterImpl<W>,
}

impl<W: Write + Send> CsvStructWriter<W> {
  pub fn new(writer: W, delimiter: u8) -> Result<Self, Report> {
    let writer = CsvWriterBuilder::new().delimiter(delimiter).from_writer(writer);
    Ok(Self { writer })
  }

  pub fn write<T: Serialize>(&mut self, record: &T) -> Result<(), Report> {
    self.writer.serialize(record)?;
    Ok(())
  }
}

/// Writes CSV files. Each row is a serde-annotated struct.
pub struct CsvStructFileWriter {
  pub filepath: PathBuf,
  pub writer: CsvStructWriter<Box<dyn Write + Send>>,
}

impl CsvStructFileWriter {
  pub fn new(filepath: impl AsRef<Path>, delimiter: u8) -> Result<Self, Report> {
    let filepath = filepath.as_ref();
    let file = create_file_or_stdout(filepath)?;
    let writer = CsvStructWriter::new(file, delimiter)?;
    Ok(Self {
      filepath: filepath.to_owned(),
      writer,
    })
  }

  pub fn write<T: Serialize>(&mut self, record: &T) -> Result<(), Report> {
    self.writer.write(record)?;
    Ok(())
  }
}

pub trait VecWriter {
  fn write<I: IntoIterator<Item = T>, T: AsRef<[u8]>>(&mut self, values: I) -> Result<(), Report>;
}

/// Writes CSV. Each row is a vec of strings.
pub struct CsvVecWriter<W: Write + Send> {
  pub headers: Vec<String>,
  pub writer: CsvWriterImpl<W>,
}

impl<W: Write + Send> CsvVecWriter<W> {
  pub fn new(writer: W, delimiter: u8, headers: &[String]) -> Result<Self, Report> {
    let mut writer = CsvWriterBuilder::new().delimiter(delimiter).from_writer(writer);
    writer.write_record(headers)?;
    Ok(Self {
      headers: headers.to_owned(),
      writer,
    })
  }
}

impl<W: Write + Send> VecWriter for CsvVecWriter<W> {
  fn write<I: IntoIterator<Item = T>, T: AsRef<[u8]>>(&mut self, values: I) -> Result<(), Report> {
    self.writer.write_record(values)?;
    Ok(())
  }
}

/// Writes CSV files. Each row is a vec of strings.
pub struct CsvVecFileWriter {
  pub filepath: PathBuf,
  pub headers: Vec<String>,
  pub writer: CsvVecWriter<Box<dyn Write + Send>>,
}

impl CsvVecFileWriter {
  pub fn new(filepath: impl AsRef<Path>, delimiter: u8, headers: &[String]) -> Result<Self, Report> {
    let filepath = filepath.as_ref();
    let file = create_file_or_stdout(filepath)?;
    let writer = CsvVecWriter::new(file, delimiter, headers)?;
    Ok(Self {
      filepath: filepath.to_owned(),
      headers: headers.to_owned(),
      writer,
    })
  }
}

impl VecWriter for CsvVecFileWriter {
  fn write<I: IntoIterator<Item = T>, T: AsRef<[u8]>>(&mut self, values: I) -> Result<(), Report> {
    self.writer.write(values)?;
    Ok(())
  }
}

/// Parse entire CSV file
pub fn csv_read_file<T: for<'de> Deserialize<'de>>(filepath: impl AsRef<Path>) -> Result<Vec<T>, Report> {
  let filepath = filepath.as_ref();
  let data = read_file_to_string(filepath)?;
  csv_read_str(data)
}

/// Parse entire CSV string
pub fn csv_read_str<T: for<'de> Deserialize<'de>, S: AsRef<str>>(data: S) -> Result<Vec<T>, Report> {
  let reader = CsvReaderBuilder::new()
    .has_headers(true)
    .from_reader(data.as_ref().as_bytes());
  reader
    .into_deserialize::<T>()
    .map(to_eyre_error)
    .collect::<Result<Vec<T>, Report>>()
}

pub fn get_col_name(
  headers: &[String],
  possible_names: &[String],
  provided_name: &Option<String>,
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
    headers
      .iter()
      .position(|header| possible_names.contains(header))
      .ok_or_else(|| {
        eyre!(
          "Unable to find column:\n  Looking for: {}\n  Available columns are: {}",
          possible_names.join(", "),
          headers.join(", ")
        )
      })
  }
}

pub fn guess_csv_delimiter(filepath: impl AsRef<Path>) -> Result<u8, Report> {
  let filepath = filepath.as_ref();
  let ext = extension(filepath)
    .ok_or_else(|| eyre!("Unable to detect file extension: '{filepath:?}': "))?
    .to_lowercase();
  match ext.as_str() {
    "csv" => Ok(b','),
    "tsv" => Ok(b'\t'),
    "ssv" => Ok(b';'),
    _ => make_error!("Unknown file extension: '{ext}'"),
  }
}

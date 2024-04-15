use crate::io::compression::Decompressor;
use crate::io::concat::concat;
use crate::io::file::{create_file, open_file_or_stdin, open_stdin};
use crate::make_error;
use eyre::Report;
use log::info;
use serde::{Deserialize, Serialize};
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Default, Debug, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct FastaRecord {
  pub seq_name: String,
  pub desc: String,
  pub seq: String,
  pub index: usize,
}

impl FastaRecord {
  pub fn new() -> Self {
    Self::default()
  }

  pub fn clear(&mut self) {
    self.seq_name.clear();
    self.seq.clear();
    self.index = 0;
  }

  pub fn is_empty(&self) -> bool {
    self.seq_name.is_empty() && self.seq_name.is_empty() && self.index == 0
  }
}

pub struct FastaReader<'a> {
  reader: Box<dyn BufRead + 'a>,
  line: String,
  index: usize,
}

impl<'a> FastaReader<'a> {
  pub fn new(reader: Box<dyn BufRead + 'a>) -> Self {
    Self {
      reader,
      line: String::new(),
      index: 0,
    }
  }

  pub fn from_str(contents: &'a str) -> Result<Self, Report> {
    let reader = contents.as_bytes();
    Ok(Self::new(Box::new(reader)))
  }

  pub fn from_str_and_path(contents: &'static str, filepath: impl AsRef<Path>) -> Result<Self, Report> {
    let decompressor = Decompressor::from_str_and_path(contents, filepath)?;
    let reader = BufReader::new(decompressor);
    Ok(Self::new(Box::new(reader)))
  }

  /// Reads multiple files sequentially given a set of paths
  pub fn from_paths<P: AsRef<Path>>(filepaths: &[P]) -> Result<Self, Report> {
    if filepaths.is_empty() {
      info!("Reading input fasta from standard input");
      return Ok(Self::new(open_stdin()?));
    }

    let readers: Vec<Box<dyn BufRead + 'a>> = filepaths
      .iter()
      .map(|filepath| -> Result<Box<dyn BufRead + 'a>, Report> { open_file_or_stdin(&Some(filepath)) })
      .collect::<Result<Vec<Box<dyn BufRead + 'a>>, Report>>()?;

    let concat = concat(readers.into_iter());
    let concat_buf = BufReader::new(concat);

    Ok(Self::new(Box::new(concat_buf)))
  }

  #[allow(clippy::string_slice)]
  pub fn read(&mut self, record: &mut FastaRecord) -> Result<(), Report> {
    record.clear();

    if self.line.is_empty() {
      self.reader.read_line(&mut self.line)?;
      if self.line.is_empty() {
        return Ok(());
      }
    }

    if !self.line.starts_with('>') {
      return make_error!("Expected character '>' at record start.");
    }

    let (seq_name, desc) = self.line[1..]
      .split_once(char::is_whitespace)
      .unwrap_or_else(|| (&self.line[1..], ""));

    record.seq_name = seq_name.trim().to_owned();
    record.desc = desc.trim().to_owned();

    loop {
      self.line.clear();
      self.reader.read_line(&mut self.line)?;
      if self.line.is_empty() || self.line.starts_with('>') {
        break;
      }

      let fragment = self.line.trim_end().chars().map(|c| c.to_ascii_uppercase());

      record.seq.extend(fragment);
    }

    record.index = self.index;
    self.index += 1;

    Ok(())
  }
}

pub fn read_one_fasta(filepath: impl AsRef<Path>) -> Result<FastaRecord, Report> {
  let filepath = filepath.as_ref();
  let mut reader = FastaReader::from_paths(&[filepath])?;
  let mut record = FastaRecord::default();
  reader.read(&mut record)?;
  Ok(record)
}

pub fn read_many_fasta<P: AsRef<Path>>(filepaths: &[P]) -> Result<Vec<FastaRecord>, Report> {
  let mut reader = FastaReader::from_paths(filepaths)?;
  let mut fasta_records = Vec::<FastaRecord>::new();

  loop {
    let mut record = FastaRecord::default();
    reader.read(&mut record).unwrap();
    if record.is_empty() {
      break;
    }
    fasta_records.push(record);
  }

  Ok(fasta_records)
}

pub fn read_one_fasta_str(contents: &str) -> Result<FastaRecord, Report> {
  let mut reader = FastaReader::from_str(contents)?;
  let mut record = FastaRecord::default();
  reader.read(&mut record)?;
  Ok(record)
}

// Writes sequences into given fasta file
pub struct FastaWriter {
  writer: Box<dyn std::io::Write>,
}

impl FastaWriter {
  pub fn new(writer: Box<dyn std::io::Write>) -> Self {
    Self { writer }
  }

  pub fn from_path(filepath: impl AsRef<Path>) -> Result<Self, Report> {
    Ok(Self::new(create_file(filepath)?))
  }

  pub fn write(&mut self, seq_name: &str, seq: &str) -> Result<(), Report> {
    write!(self.writer, ">{seq_name}\n{seq}\n")?;
    Ok(())
  }

  pub fn flush(&mut self) -> Result<(), Report> {
    self.writer.flush()?;
    Ok(())
  }
}

#[derive(Clone, Debug, Serialize)]
struct OutputTranslationsTemplateContext<'a> {
  gene: &'a str,
}

use crate::io::compression::Decompressor;
use crate::io::concat::concat;
use crate::make_error;
use eyre::{Report, WrapErr};
use log::{info, warn};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{stdin, BufRead, BufReader, Read};
use std::path::Path;
use std::str::FromStr;

use crate::io::file::create_file;
#[cfg(not(target_arch = "wasm32"))]
use atty::{is as is_tty, Stream};

const TTY_WARNING: &str = r#"Reading from standard input which is a TTY (e.g. an interactive terminal). This is likely not what you meant. Instead:

 - if you want to read fasta from the output of another program, try:

    cat sequences.fasta | treetime <your other flags>

 - if you want to read fasta from file(s), try:

    treetime sequences.fasta sequences2.fasta <your other flags>
"#;

#[derive(Clone, Default, Debug, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct FastaRecord {
  pub seq_name: String,
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

      #[cfg(not(target_arch = "wasm32"))]
      if is_tty(Stream::Stdin) {
        warn!("{TTY_WARNING}");
      }

      return Ok(Self::new(Box::new(BufReader::new(stdin()))));
    }

    let readers: Vec<Box<dyn BufRead + 'a>> = filepaths
      .iter()
      .map(|filepath| -> Result<Box<dyn BufRead + 'a>, Report> {
        let filepath = filepath.as_ref();
        let file = File::open(&filepath).wrap_err_with(|| format!("When opening file: {filepath:?}"))?;
        let decompressor = Decompressor::from_path(file, &filepath)?;
        let reader = BufReader::with_capacity(32 * 1024, decompressor);
        Ok(Box::new(reader))
      })
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

    record.seq_name = self.line[1..].trim().to_owned();

    loop {
      self.line.clear();
      self.reader.read_line(&mut self.line)?;
      if self.line.is_empty() || self.line.starts_with('>') {
        break;
      }

      let fragment = self.line.trim_end().chars().into_iter().map(|c| c.to_ascii_uppercase());

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

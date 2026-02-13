use crate::alphabet::alphabet::Alphabet;
use crate::io::concat::Concat;
use crate::make_error;
use eyre::{Context, Report};
use itertools::Itertools;
use log::warn;
use serde::{Deserialize, Serialize};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use treetime_primitives::{AsciiChar, Seq};
use treetime_utils::compression::Decompressor;
use treetime_utils::file::{create_file_or_stdout, open_file_or_stdin};
use treetime_utils::string::quote_single;

#[derive(Clone, Default, Debug, Deserialize, Serialize, Eq, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct FastaRecord {
  pub seq_name: String,
  pub desc: Option<String>,
  pub seq: Seq,
  pub index: usize,
}

impl FastaRecord {
  pub fn new() -> Self {
    Self::default()
  }

  pub fn clear(&mut self) {
    self.seq_name.clear();
    self.desc = None;
    self.seq.clear();
    self.index = 0;
  }

  pub fn is_empty(&self) -> bool {
    self.seq_name.is_empty() && self.seq_name.is_empty() && self.desc.is_none() && self.index == 0
  }

  pub fn header(&self) -> String {
    match &self.desc {
      Some(desc) => format!(">{} {}", self.seq_name, desc),
      None => format!(">{}", self.seq_name),
    }
  }
}

pub struct FastaReader<'a, 'b> {
  reader: Box<dyn BufRead + 'a>,
  alphabet: &'b Alphabet,
  line: String,
  n_lines: usize,
  n_chars: usize,
  index: usize,
}

impl<'a, 'b> FastaReader<'a, 'b> {
  pub fn new(reader: Box<dyn BufRead + 'a>, alphabet: &'b Alphabet) -> Self {
    Self {
      reader,
      alphabet,
      line: String::new(),
      n_lines: 0,
      n_chars: 0,
      index: 0,
    }
  }

  pub fn from_str(contents: &'a impl AsRef<str>, alphabet: &'b Alphabet) -> Result<Self, Report> {
    let reader = contents.as_ref().as_bytes();
    Ok(Self::new(Box::new(reader), alphabet))
  }

  pub fn from_str_and_path(
    contents: &'static str,
    filepath: impl AsRef<Path>,
    alphabet: &'b Alphabet,
  ) -> Result<Self, Report> {
    let decompressor = Decompressor::from_str_and_path(contents, filepath)?;
    let reader = BufReader::new(decompressor);
    Ok(Self::new(Box::new(reader), alphabet))
  }

  pub fn from_path(filepath: impl AsRef<Path>, alphabet: &'b Alphabet) -> Result<Self, Report> {
    Self::from_paths(&[filepath], alphabet)
  }

  /// Reads multiple files sequentially given a set of paths
  pub fn from_paths<P: AsRef<Path>>(filepaths: &[P], alphabet: &'b Alphabet) -> Result<Self, Report> {
    let readers: Vec<Box<dyn BufRead + 'a>> = filepaths
      .iter()
      .map(|filepath| -> Result<Box<dyn BufRead + 'a>, Report> { open_file_or_stdin(&Some(filepath)) })
      .collect::<Result<Vec<Box<dyn BufRead + 'a>>, Report>>()?;

    let concat = Concat::with_delimiter(readers.into_iter(), Some(b"\n".to_vec()));
    let concat_buf = BufReader::new(concat);

    Ok(Self::new(Box::new(concat_buf), alphabet))
  }

  #[allow(clippy::string_slice)]
  pub fn read(&mut self, record: &mut FastaRecord) -> Result<(), Report> {
    record.clear();

    if self.line.is_empty() {
      // Read lines until we find the next record or EOF
      loop {
        self.line.clear();
        if self.reader.read_line(&mut self.line)? == 0 {
          if self.index > 0 {
            // We have read at least one record  by the end of input - this is normal operation mode
            return Ok(());
          }

          if self.index == 0 && self.n_chars == 0 {
            // We have read no records and no non-whitespace characters by the end of input
            warn!(
              "FASTA input is empty or consists entirely from whitespace: this is allowed but might not be what's intended"
            );
            return Ok(());
          }

          // We have read some characters, but no records detected by the end of input
          return make_error!(
            "FASTA input is incorrectly formatted: expected at least one FASTA record starting with character '>', but none found"
          );
        }

        let trimmed = self.line.trim();
        self.n_lines += 1;
        self.n_chars += trimmed.len();
        if trimmed.starts_with('>') {
          break; // Found the header of the next record
        }
      }
    }

    let header_line = self.line.trim();
    let (name, desc) = header_line[1..]
      .split_once(' ')
      .unwrap_or_else(|| (&header_line[1..], ""));
    record.seq_name = name.to_owned();
    record.desc = if desc.is_empty() { None } else { Some(desc.to_owned()) };
    record.index = self.index;
    self.index += 1;

    // Read sequence lines until the next header or EOF
    self.line.clear();
    while self.reader.read_line(&mut self.line)? > 0 {
      let trimmed = self.line.trim();
      self.n_lines += 1;
      self.n_chars += trimmed.len();
      if trimmed.starts_with('>') {
        // We have reached the next record
        break;
      }

      record.seq.reserve(trimmed.len());
      for c in trimmed.chars() {
        let uc = AsciiChar::from(c.to_ascii_uppercase());
        if self.alphabet.contains(uc) {
          record.seq.push(uc);
        } else {
          return make_error!(
            "FASTA input is incorrect: character \"{c}\" is not in the alphabet. Expected characters: {}",
            self.alphabet.chars().map(char::from).map(quote_single).join(", ")
          )
          .wrap_err_with(|| format!("When processing sequence #{}: \"{}\"", self.index, record.header()));
        }
      }

      self.line.clear();
    }

    Ok(())
  }
}

pub fn read_one_fasta(filepath: impl AsRef<Path>, alphabet: &Alphabet) -> Result<FastaRecord, Report> {
  let filepath = filepath.as_ref();
  let mut reader = FastaReader::from_path(filepath, alphabet)?;
  let mut record = FastaRecord::default();
  reader.read(&mut record)?;
  Ok(record)
}

pub fn read_many_fasta<P: AsRef<Path>>(filepaths: &[P], alphabet: &Alphabet) -> Result<Vec<FastaRecord>, Report> {
  let mut reader = FastaReader::from_paths(filepaths, alphabet)?;
  let mut fasta_records = Vec::<FastaRecord>::new();

  loop {
    let mut record = FastaRecord::default();
    reader.read(&mut record)?;
    if record.is_empty() {
      break;
    }
    fasta_records.push(record);
  }

  Ok(fasta_records)
}

pub fn read_one_fasta_str(contents: impl AsRef<str>, alphabet: &Alphabet) -> Result<FastaRecord, Report> {
  let mut reader = FastaReader::from_str(&contents, alphabet)?;
  let mut record = FastaRecord::default();
  reader.read(&mut record)?;
  Ok(record)
}

pub fn read_many_fasta_str(contents: impl AsRef<str>, alphabet: &Alphabet) -> Result<Vec<FastaRecord>, Report> {
  let mut reader = FastaReader::from_str(&contents, alphabet)?;
  let mut fasta_records = Vec::<FastaRecord>::new();

  loop {
    let mut record = FastaRecord::default();
    reader.read(&mut record)?;
    if record.is_empty() {
      break;
    }
    fasta_records.push(record);
  }

  Ok(fasta_records)
}

// Writes sequences into given fasta file
pub struct FastaWriter {
  writer: Box<dyn Write>,
}

impl FastaWriter {
  pub fn new(writer: Box<dyn Write>) -> Self {
    Self { writer }
  }

  pub fn from_path(filepath: impl AsRef<Path>) -> Result<Self, Report> {
    Ok(Self::new(create_file_or_stdout(filepath)?))
  }

  pub fn write(&mut self, seq_name: impl AsRef<str>, desc: &Option<String>, seq: &Seq) -> Result<(), Report> {
    self.writer.write_all(b">")?;
    self.writer.write_all(seq_name.as_ref().as_bytes())?;

    if let Some(desc) = desc {
      self.writer.write_all(b" ")?;
      self.writer.write_all(desc.as_bytes())?;
    }

    self.writer.write_all(b"\n")?;
    self.writer.write_all(seq.as_ref())?;
    self.writer.write_all(b"\n")?;
    Ok(())
  }

  pub fn flush(&mut self) -> Result<(), Report> {
    self.writer.flush()?;
    Ok(())
  }
}

pub fn write_one_fasta(
  filepath: impl AsRef<Path>,
  seq_name: impl AsRef<str>,
  desc: &Option<String>,
  seq: &Seq,
) -> Result<(), Report> {
  let mut writer = FastaWriter::from_path(&filepath)?;
  writer.write(seq_name, desc, seq)
}

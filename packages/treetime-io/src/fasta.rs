use crate::concat::Concat;
use eyre::{Context, Report};
use itertools::Itertools;
use log::warn;
use serde::{Deserialize, Serialize};
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use treetime_primitives::{AlignmentRecord, AlphabetLike, AsciiChar, Seq};
use treetime_utils::fmt::string::quote_single;
use treetime_utils::io::file::open_file_or_stdin;
use treetime_utils::make_error;

pub fn read_many_fasta_path<P: AsRef<Path>, A: AlphabetLike>(
  filepaths: &[P],
  alphabet: &A,
) -> Result<Vec<FastaRecord>, Report> {
  let reader = FastaReader::from_paths(filepaths, alphabet)?;
  read_many_fasta(reader)
}

#[cfg_attr(
  dylint_lib = "treetime_lints",
  allow(pub_unused_in_workspace, reason = "used only by tests of other workspace crates, which a cfg(test) item cannot reach")
)]
pub fn read_many_fasta_str<A: AlphabetLike>(
  contents: impl AsRef<str>,
  alphabet: &A,
) -> Result<Vec<FastaRecord>, Report> {
  let reader = FastaReader::from_str(&contents, alphabet)?;
  read_many_fasta(reader)
}

pub fn read_many_fasta<A: AlphabetLike>(mut reader: FastaReader<'_, '_, A>) -> Result<Vec<FastaRecord>, Report> {
  let mut records = Vec::new();

  loop {
    let mut record = FastaRecord::default();
    reader.read(&mut record)?;
    if record.is_empty() {
      break;
    }
    records.push(record);
  }

  Ok(records)
}

pub struct FastaReader<'a, 'b, A: AlphabetLike> {
  reader: Box<dyn BufRead + 'a>,
  alphabet: &'b A,
  line: String,
  n_lines: usize,
  n_chars: usize,
  index: usize,
}

impl<'a, 'b, A: AlphabetLike> FastaReader<'a, 'b, A> {
  pub fn new(reader: Box<dyn BufRead + 'a>, alphabet: &'b A) -> Self {
    Self {
      reader,
      alphabet,
      line: String::new(),
      n_lines: 0,
      n_chars: 0,
      index: 0,
    }
  }

  fn from_str(contents: &'a impl AsRef<str>, alphabet: &'b A) -> Result<Self, Report> {
    let reader = contents.as_ref().as_bytes();
    Ok(Self::new(Box::new(reader), alphabet))
  }

  fn from_paths<P: AsRef<Path>>(filepaths: &[P], alphabet: &'b A) -> Result<Self, Report> {
    let readers: Vec<Box<dyn BufRead + 'a>> = filepaths
      .iter()
      .map(|filepath| -> Result<Box<dyn BufRead + 'a>, Report> { open_file_or_stdin(&Some(filepath)) })
      .collect::<Result<Vec<Box<dyn BufRead + 'a>>, Report>>()?;

    let concat = Concat::with_delimiter(readers.into_iter(), Some(b"\n".to_vec()));
    let concat_buf = BufReader::new(concat);

    Ok(Self::new(Box::new(concat_buf), alphabet))
  }

  #[expect(
    clippy::string_slice,
    reason = "index 1 follows the ASCII '>' marker, a char boundary"
  )]
  pub fn read(&mut self, record: &mut FastaRecord) -> Result<(), Report> {
    record.clear();

    if self.line.is_empty() {
      loop {
        self.line.clear();
        let n_read = self
          .reader
          .read_line(&mut self.line)
          .wrap_err_with(|| format!("When reading line {} of FASTA input", self.n_lines + 1))?;
        if n_read == 0 {
          if self.index > 0 {
            return Ok(());
          }

          if self.index == 0 && self.n_chars == 0 {
            warn!(
              "FASTA input is empty or consists entirely from whitespace: this is allowed but might not be what's intended"
            );
            return Ok(());
          }

          return make_error!(
            "FASTA input is incorrectly formatted: expected at least one FASTA record starting with character '>', but none found"
          );
        }

        let trimmed = self.line.trim();
        self.n_lines += 1;
        self.n_chars += trimmed.len();
        if trimmed.starts_with('>') {
          break;
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

    self.line.clear();
    while self
      .reader
      .read_line(&mut self.line)
      .wrap_err_with(|| format!("When reading line {} of FASTA input", self.n_lines + 1))?
      > 0
    {
      let trimmed = self.line.trim();
      self.n_lines += 1;
      self.n_chars += trimmed.len();
      if trimmed.starts_with('>') {
        break;
      }

      record.seq.reserve(trimmed.len());
      for c in trimmed.chars() {
        let uc = AsciiChar::try_from_char(c.to_ascii_uppercase())
          .wrap_err_with(|| format!("When processing sequence #{}: \"{}\"", self.index, record.header()))?;
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

#[derive(Clone, Default, Debug, Deserialize, Serialize, Eq, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct FastaRecord {
  pub seq_name: String,
  pub desc: Option<String>,
  pub seq: Seq,
  pub index: usize,
}

impl From<FastaRecord> for AlignmentRecord {
  fn from(record: FastaRecord) -> Self {
    Self {
      name: record.seq_name,
      seq: record.seq,
    }
  }
}

impl FastaRecord {
  #[cfg_attr(
    dylint_lib = "treetime_lints",
    allow(pub_unused_in_workspace, reason = "used only by tests of other workspace crates, which a cfg(test) item cannot reach")
  )]
  pub fn new() -> Self {
    Self::default()
  }

  fn clear(&mut self) {
    self.seq_name.clear();
    self.desc = None;
    self.seq.clear();
    self.index = 0;
  }

  pub fn is_empty(&self) -> bool {
    self.seq_name.is_empty() && self.seq.is_empty() && self.desc.is_none() && self.index == 0
  }

  fn header(&self) -> String {
    match &self.desc {
      Some(desc) => format!(">{} {}", self.seq_name, desc),
      None => format!(">{}", self.seq_name),
    }
  }
}

pub struct FastaWriter {
  writer: Box<dyn Write>,
}

impl FastaWriter {
  pub fn new(writer: Box<dyn Write>) -> Self {
    Self { writer }
  }

  pub fn write(&mut self, seq_name: impl AsRef<str>, desc: &Option<String>, seq: &Seq) -> Result<(), Report> {
    let seq_name = seq_name.as_ref();
    write_fasta_record(&mut self.writer, seq_name, desc.as_deref(), seq)
      .wrap_err_with(|| format!("When writing FASTA record '{seq_name}'"))
  }
}

fn write_fasta_record(writer: &mut impl Write, seq_name: &str, desc: Option<&str>, seq: &Seq) -> io::Result<()> {
  writer.write_all(b">")?;
  writer.write_all(seq_name.as_bytes())?;
  if let Some(desc) = desc {
    writer.write_all(b" ")?;
    writer.write_all(desc.as_bytes())?;
  }
  writer.write_all(b"\n")?;
  writer.write_all(seq.as_ref())?;
  writer.write_all(b"\n")
}

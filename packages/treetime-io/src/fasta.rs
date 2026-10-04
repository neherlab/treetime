use eyre::{Context, Report};
use itertools::Itertools;
use log::warn;
use std::io::{self, BufRead, Write};
use std::path::Path;
use treetime_primitives::{AlignmentRecord, AlphabetLike, AsciiChar, Seq};
use treetime_utils::fmt::string::quote_single;
use treetime_utils::io::file::{FileWriter, create_file_or_stdout, read_file_with};
use treetime_utils::make_error;

pub const FASTA_EXTENSIONS: [&str; 4] = ["fasta", "fa", "fas", "aln"];

pub fn fasta_read_file<A: AlphabetLike>(filepath: impl AsRef<Path>, alphabet: &A) -> Result<Vec<FastaRecord>, Report> {
  read_file_with(filepath, |reader| fasta_read(reader, alphabet))
}

pub fn fasta_read<A: AlphabetLike>(reader: impl BufRead, alphabet: &A) -> Result<Vec<FastaRecord>, Report> {
  FastaRecords::new(reader, alphabet).collect()
}

#[derive(Clone, Default, Debug, Eq, PartialEq)]
pub struct FastaRecord {
  pub seq_name: String,
  pub desc: Option<String>,
  pub seq: Seq,
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
  fn header(&self) -> String {
    match &self.desc {
      Some(desc) => format!(">{} {}", self.seq_name, desc),
      None => format!(">{}", self.seq_name),
    }
  }
}

pub struct FastaWriter {
  writer: FileWriter,
}

impl FastaWriter {
  pub fn create(filepath: impl AsRef<Path>) -> Result<Self, Report> {
    Ok(Self {
      writer: create_file_or_stdout(filepath)?,
    })
  }

  pub fn write(&mut self, seq_name: &str, desc: Option<&str>, seq: &Seq) -> Result<(), Report> {
    fasta_write_record(&mut self.writer, seq_name, desc, seq)
      .wrap_err_with(|| format!("When writing FASTA record '{seq_name}'"))
      .wrap_err_with(|| format!("When writing file '{}'", self.writer.filepath().display()))
  }

  pub fn finish(self) -> Result<(), Report> {
    self.writer.finish()
  }
}

struct FastaRecords<'a, R: BufRead, A: AlphabetLike> {
  reader: R,
  alphabet: &'a A,
  accepted: [bool; 128],
  line: String,
  n_lines: usize,
  n_chars: usize,
  n_records: usize,
}

impl<'a, R: BufRead, A: AlphabetLike> FastaRecords<'a, R, A> {
  fn new(reader: R, alphabet: &'a A) -> Self {
    Self {
      reader,
      alphabet,
      accepted: accepted_bytes(alphabet),
      line: String::new(),
      n_lines: 0,
      n_chars: 0,
      n_records: 0,
    }
  }

  fn read_line(&mut self) -> Result<bool, Report> {
    self.line.clear();
    let n_read = self
      .reader
      .read_line(&mut self.line)
      .wrap_err_with(|| format!("When reading line {} of FASTA input", self.n_lines + 1))?;
    if n_read > 0 {
      self.n_lines += 1;
      self.n_chars += self.line.trim().len();
    }
    Ok(n_read > 0)
  }

  fn find_header(&mut self) -> Result<bool, Report> {
    if self.line.trim_start().starts_with('>') {
      return Ok(true);
    }
    while self.read_line()? {
      if self.line.trim_start().starts_with('>') {
        return Ok(true);
      }
    }
    if self.n_records > 0 {
      return Ok(false);
    }
    if self.n_chars == 0 {
      warn!(
        "FASTA input is empty or consists entirely from whitespace: this is allowed but might not be what's intended"
      );
      return Ok(false);
    }
    make_error!(
      "FASTA input is incorrectly formatted: expected at least one FASTA record starting with character '>', but none found"
    )
  }

  #[expect(
    clippy::string_slice,
    reason = "index 1 follows the ASCII '>' marker, a char boundary"
  )]
  fn read_record(&mut self) -> Result<Option<FastaRecord>, Report> {
    if !self.find_header()? {
      return Ok(None);
    }

    let header_line = self.line.trim();
    let (name, desc) = header_line[1..]
      .split_once(' ')
      .unwrap_or_else(|| (&header_line[1..], ""));
    let mut record = FastaRecord {
      seq_name: name.to_owned(),
      desc: (!desc.is_empty()).then(|| desc.to_owned()),
      seq: Seq::new(),
    };
    self.n_records += 1;

    while self.read_line()? {
      if self.line.trim_start().starts_with('>') {
        break;
      }
      self.append_sequence_line(&mut record)?;
    }

    Ok(Some(record))
  }

  fn append_sequence_line(&self, record: &mut FastaRecord) -> Result<(), Report> {
    let trimmed = self.line.trim();
    record.seq.reserve(trimmed.len());
    if trimmed.is_ascii() && trimmed.bytes().all(|byte| self.accepted[usize::from(byte)]) {
      record.seq.extend(
        trimmed
          .bytes()
          .map(|byte| AsciiChar::from_byte_unchecked(byte.to_ascii_uppercase())),
      );
      return Ok(());
    }
    for c in trimmed.chars() {
      let uc = AsciiChar::try_from_char(c.to_ascii_uppercase())
        .wrap_err_with(|| format!("When processing sequence #{}: \"{}\"", self.n_records, record.header()))?;
      if self.alphabet.contains(uc) {
        record.seq.push(uc);
      } else {
        return make_error!(
          "FASTA input is incorrect: character \"{c}\" is not in the alphabet. Expected characters: {}",
          self.alphabet.chars().map(char::from).map(quote_single).join(", ")
        )
        .wrap_err_with(|| format!("When processing sequence #{}: \"{}\"", self.n_records, record.header()));
      }
    }
    Ok(())
  }
}

impl<R: BufRead, A: AlphabetLike> Iterator for FastaRecords<'_, R, A> {
  type Item = Result<FastaRecord, Report>;

  fn next(&mut self) -> Option<Self::Item> {
    self.read_record().transpose()
  }
}

#[allow(
  clippy::as_conversions,
  reason = "table indices below 128 convert exactly to ASCII bytes"
)]
fn accepted_bytes<A: AlphabetLike>(alphabet: &A) -> [bool; 128] {
  std::array::from_fn(|byte| alphabet.contains(AsciiChar::from_byte_unchecked((byte as u8).to_ascii_uppercase())))
}

pub(crate) fn fasta_write_record(
  writer: &mut impl Write,
  seq_name: &str,
  desc: Option<&str>,
  seq: &Seq,
) -> io::Result<()> {
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

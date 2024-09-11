use crate::io::compression::Decompressor;
use crate::io::concat::Concat;
use crate::io::file::{create_file_or_stdout, open_file_or_stdin, open_stdin};
use crate::make_error;
use eyre::Report;
use log::info;
use serde::{Deserialize, Serialize};
use std::io::{BufRead, BufReader};
use std::path::Path;

pub fn is_char_allowed(c: char) -> bool {
  c.is_ascii()
}

#[derive(Clone, Default, Debug, Deserialize, Serialize, Eq, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct FastaRecord {
  pub seq_name: String,
  pub desc: Option<String>,
  pub seq: String,
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

  pub fn from_str(contents: &'a impl AsRef<str>) -> Result<Self, Report> {
    let reader = contents.as_ref().as_bytes();
    Ok(Self::new(Box::new(reader)))
  }

  pub fn from_str_and_path(contents: &'static str, filepath: impl AsRef<Path>) -> Result<Self, Report> {
    let decompressor = Decompressor::from_str_and_path(contents, filepath)?;
    let reader = BufReader::new(decompressor);
    Ok(Self::new(Box::new(reader)))
  }

  pub fn from_path(filepath: impl AsRef<Path>) -> Result<Self, Report> {
    Self::from_paths(&[filepath])
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

    let concat = Concat::with_delimiter(readers.into_iter(), Some(b"\n".to_vec()));
    let concat_buf = BufReader::new(concat);

    Ok(Self::new(Box::new(concat_buf)))
  }

  #[allow(clippy::string_slice)]
  pub fn read(&mut self, record: &mut FastaRecord) -> Result<(), Report> {
    record.clear();

    if self.line.is_empty() {
      loop {
        self.line.clear();

        let n_bytes = self.reader.read_line(&mut self.line)?;
        if n_bytes == 0 {
          return Ok(());
        }

        let trimmed = self.line.trim();
        if !trimmed.is_empty() {
          self.line = trimmed.to_owned();
          break;
        }
      }
    }

    if !self.line.starts_with('>') {
      return make_error!("Expected character '>' at record start.");
    }

    let s = self.line[1..].trim();

    let mut parts = s.splitn(2, ' ');
    let name = parts.next().unwrap_or_default().to_owned();
    let desc = parts.next().map(ToOwned::to_owned);

    record.seq_name = name;
    record.desc = desc;

    loop {
      self.line.clear();

      let n_bytes = self.reader.read_line(&mut self.line)?;
      if n_bytes == 0 {
        record.index = self.index;
        self.index += 1;
        return Ok(());
      }

      let trimmed = self.line.trim();
      if !trimmed.is_empty() {
        self.line = trimmed.to_owned();
        break;
      }
    }

    if self.line.is_empty() || self.line.starts_with('>') {
      record.index = self.index;
      self.index += 1;
      return Ok(());
    }

    let fragment = self
      .line
      .chars()
      .filter(|c| is_char_allowed(*c))
      .map(|c| c.to_ascii_uppercase());

    record.seq.extend(fragment);

    loop {
      self.line.clear();
      self.reader.read_line(&mut self.line)?;
      self.line = self.line.trim().to_owned();
      if self.line.is_empty() || self.line.starts_with('>') {
        break;
      }

      let fragment = self
        .line
        .chars()
        .filter(|c| is_char_allowed(*c))
        .map(|c| c.to_ascii_uppercase());

      record.seq.extend(fragment);
    }

    record.index = self.index;
    self.index += 1;

    Ok(())
  }
}

pub fn read_one_fasta(filepath: impl AsRef<Path>) -> Result<FastaRecord, Report> {
  let filepath = filepath.as_ref();
  let mut reader = FastaReader::from_path(filepath)?;
  let mut record = FastaRecord::default();
  reader.read(&mut record)?;
  Ok(record)
}

pub fn read_many_fasta<P: AsRef<Path>>(filepaths: &[P]) -> Result<Vec<FastaRecord>, Report> {
  let mut reader = FastaReader::from_paths(filepaths)?;
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

pub fn read_one_fasta_str(contents: impl AsRef<str>) -> Result<FastaRecord, Report> {
  let mut reader = FastaReader::from_str(&contents)?;
  let mut record = FastaRecord::default();
  reader.read(&mut record)?;
  Ok(record)
}

pub fn read_many_fasta_str(contents: impl AsRef<str>) -> Result<Vec<FastaRecord>, Report> {
  let mut reader = FastaReader::from_str(&contents)?;
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
  writer: Box<dyn std::io::Write>,
}

impl FastaWriter {
  pub fn new(writer: Box<dyn std::io::Write>) -> Self {
    Self { writer }
  }

  pub fn from_path(filepath: impl AsRef<Path>) -> Result<Self, Report> {
    Ok(Self::new(create_file_or_stdout(filepath)?))
  }

  pub fn write(
    &mut self,
    seq_name: impl AsRef<str>,
    desc: &Option<String>,
    seq: impl AsRef<str>,
  ) -> Result<(), Report> {
    self.writer.write_all(b">")?;
    self.writer.write_all(seq_name.as_ref().as_bytes())?;

    if let Some(desc) = desc {
      self.writer.write_all(b" ")?;
      self.writer.write_all(desc.as_bytes())?;
    }

    self.writer.write_all(b"\n")?;
    self.writer.write_all(seq.as_ref().as_bytes())?;
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
  seq: impl AsRef<str>,
) -> Result<(), Report> {
  let mut writer = FastaWriter::from_path(&filepath)?;
  writer.write(seq_name, desc, seq)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::o;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::io::Cursor;

  #[test]
  fn test_fasta_reader_fail_on_non_fasta() {
    let data =
      b"This is not a valid FASTA string.\nIt is not empty, and not entirely whitespace\nbut does not contain 'greater than' character.\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));
    let mut record = FastaRecord::new();
    assert_eq!(
      reader.read(&mut record).unwrap_err().to_string(),
      "Expected character '>' at record start."
    );
  }

  #[test]
  fn test_fasta_reader_read_empty() {
    let data = b"";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert!(record.is_empty());
  }

  #[test]
  fn test_fasta_reader_read_whitespace_only() {
    let data = b"\n \n \n\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert!(record.is_empty());
  }

  #[test]
  fn test_fasta_reader_read_single_record() {
    let data = b">seq1\nATCG\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_read_single_record_with_leading_newline() {
    let data = b"\n>seq1\nATCG\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_read_single_record_with_multiple_leading_newlines() {
    let data = b"\n\n\n>seq1\nATCG\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_read_single_record_without_trailing_newline() {
    let data = b">seq1\nATCG";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_read_multiple_records() {
    let data = b">seq1\nATCG\n>seq2\nGCTA\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record1 = FastaRecord::new();
    reader.read(&mut record1).unwrap();

    let mut record2 = FastaRecord::new();
    reader.read(&mut record2).unwrap();

    assert_eq!(record1.seq_name, "seq1");
    assert_eq!(record1.seq, "ATCG");
    assert_eq!(record1.index, 0);

    assert_eq!(record2.seq_name, "seq2");
    assert_eq!(record2.seq, "GCTA");
    assert_eq!(record2.index, 1);
  }

  #[test]
  fn test_fasta_reader_read_empty_lines_between_records() {
    let data = b"\n>seq1\n\nATCG\n\n\n>seq2\nGCTA\n\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record1 = FastaRecord::new();
    reader.read(&mut record1).unwrap();

    let mut record2 = FastaRecord::new();
    reader.read(&mut record2).unwrap();

    assert_eq!(record1.seq_name, "seq1");
    assert_eq!(record1.seq, "ATCG");
    assert_eq!(record1.index, 0);

    assert_eq!(record2.seq_name, "seq2");
    assert_eq!(record2.seq, "GCTA");
    assert_eq!(record2.index, 1);
  }

  #[test]
  fn test_fasta_reader_read_with_trailing_newline() {
    let data = b">seq1\nATCG\n\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_example_1() {
    let data = b"\n\n>a\nACGCTCGATC\n\n>b\nCCGCGC";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("a"),
        desc: None,
        seq: o!("ACGCTCGATC"),
        index: 0,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("b"),
        desc: None,
        seq: o!("CCGCGC"),
        index: 1,
      }
    );
  }

  #[test]
  fn test_fasta_reader_example_2() {
    let data = b">a\nACGCTCGATC\n>b\nCCGCGC\n>c";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("a"),
        desc: None,
        seq: o!("ACGCTCGATC"),
        index: 0,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("b"),
        desc: None,
        seq: o!("CCGCGC"),
        index: 1,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("c"),
        desc: None,
        seq: o!(""),
        index: 2,
      }
    );
  }

  #[test]
  fn test_fasta_reader_example_3() {
    let data = b">a\nACGCTCGATC\n>b\n>c\nCCGCGC";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("a"),
        desc: None,
        seq: o!("ACGCTCGATC"),
        index: 0,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("b"),
        desc: None,
        seq: o!(""),
        index: 1,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("c"),
        desc: None,
        seq: o!("CCGCGC"),
        index: 2,
      }
    );
  }

  #[test]
  fn test_fasta_reader_name_desc() -> Result<(), Report> {
    let actual = read_many_fasta_str(indoc! {r#"
      >Identifier Description
      ACGT
      >Identifier Description with spaces
      ACGT
    "#})?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("Identifier"),
        desc: Some(o!("Description")),
        seq: o!("ACGT"),
        index: 0,
      },
      FastaRecord {
        seq_name: o!("Identifier"),
        desc: Some(o!("Description with spaces")),
        seq: o!("ACGT"),
        index: 1,
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_fasta_reader_dedent_nuc() -> Result<(), Report> {
    let actual = read_many_fasta_str(indoc! {r#"
      >FluBuster-001
      ACAGCCATGTATTG--
      >CommonCold-AB
      ACATCCCTGTA-TG--
      >Ecoli/Joke/2024|XD
      ACATCGCCNNA--GAC

      >Sniffles-B
      GCATCCCTGTA-NG--
      >StrawberryYogurtCulture|🍓
      CCGGCCATGTATTG--
      > SneezeC-19
      CCGGCGATGTRTTG--
        >MisindentedVirus|D-skew
        TCGGCCGTGTRTTG--
    "#})?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("FluBuster-001"),
        desc: None,
        seq: o!("ACAGCCATGTATTG--"),
        index: 0,
      },
      FastaRecord {
        seq_name: o!("CommonCold-AB"),
        desc: None,
        seq: o!("ACATCCCTGTA-TG--"),
        index: 1,
      },
      FastaRecord {
        seq_name: o!("Ecoli/Joke/2024|XD"),
        desc: None,
        seq: o!("ACATCGCCNNA--GAC"),
        index: 2,
      },
      FastaRecord {
        seq_name: o!("Sniffles-B"),
        desc: None,
        seq: o!("GCATCCCTGTA-NG--"),
        index: 3,
      },
      FastaRecord {
        seq_name: o!("StrawberryYogurtCulture|🍓"),
        desc: None,
        seq: o!("CCGGCCATGTATTG--"),
        index: 4,
      },
      FastaRecord {
        seq_name: o!("SneezeC-19"),
        desc: None,
        seq: o!("CCGGCGATGTRTTG--"),
        index: 5,
      },
      FastaRecord {
        seq_name: o!("MisindentedVirus|D-skew"),
        desc: None,
        seq: o!("TCGGCCGTGTRTTG--"),
        index: 6,
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_fasta_reader_dedent_aa() -> Result<(), Report> {
    let actual = read_many_fasta_str(indoc! {r#"
      >Prot/000|β-Napkinase
      MXDXXXTQ-B--
      >Enzyme/2024|LaughzymeFactor
      AX*XB-TQVWR*

      >😊-Gigglecatalyst
      MKXTQWX-B**
      >CellFunSignal
      MQXQXXBQRW**
      >Pathway/042|Doodlease
      MXQ-*XTQWBQR
    "#})?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("Prot/000|β-Napkinase"),
        desc: None,
        seq: o!("MXDXXXTQ-B--"),
        index: 0,
      },
      FastaRecord {
        seq_name: o!("Enzyme/2024|LaughzymeFactor"),
        desc: None,
        seq: o!("AX*XB-TQVWR*"),
        index: 1,
      },
      FastaRecord {
        seq_name: o!("😊-Gigglecatalyst"),
        desc: None,
        seq: o!("MKXTQWX-B**"),
        index: 2,
      },
      FastaRecord {
        seq_name: o!("CellFunSignal"),
        desc: None,
        seq: o!("MQXQXXBQRW**"),
        index: 3,
      },
      FastaRecord {
        seq_name: o!("Pathway/042|Doodlease"),
        desc: None,
        seq: o!("MXQ-*XTQWBQR"),
        index: 4,
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_fasta_reader_multiline_and_skewed_indentation() -> Result<(), Report> {
    let actual = read_many_fasta_str(indoc! {r#"
      >MixedCaseSeq
      aCaGcCAtGtAtTG--
      >LowercaseSeq
      acagccatgtattg--
      >UppercaseSeq
      ACAGCCATGTATTG--
      >MultilineSeq
      ACAGCC
      ATGT
      ATTG--
      >SkewedIndentSeq
        ACAGCC
      ATGTATTG
       ATTG--
    "#})?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("MixedCaseSeq"),
        desc: None,
        seq: o!("ACAGCCATGTATTG--"),
        index: 0,
      },
      FastaRecord {
        seq_name: o!("LowercaseSeq"),
        desc: None,
        seq: o!("ACAGCCATGTATTG--"),
        index: 1,
      },
      FastaRecord {
        seq_name: o!("UppercaseSeq"),
        desc: None,
        seq: o!("ACAGCCATGTATTG--"),
        index: 2,
      },
      FastaRecord {
        seq_name: o!("MultilineSeq"),
        desc: None,
        seq: o!("ACAGCCATGTATTG--"),
        index: 3,
      },
      FastaRecord {
        seq_name: o!("SkewedIndentSeq"),
        desc: None,
        seq: o!("ACAGCCATGTATTGATTG--"),
        index: 4,
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }
}

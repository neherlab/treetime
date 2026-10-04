#[cfg(test)]
mod tests {
  use crate::commands::shared::alignment::read_alignment;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::fs;
  use tempfile::tempdir;
  use treetime::alphabet::alphabet::Alphabet;
  use treetime_utils::assert_error;

  #[test]
  fn test_alignment_reads_the_records_of_every_file_in_order() -> Result<(), Report> {
    let dir = tempdir()?;
    let first = dir.path().join("first.fasta");
    let second = dir.path().join("second.fasta");
    fs::write(&first, ">A\nAC\n")?;
    fs::write(&second, ">B\nGT\n>C\nTT\n")?;

    let records = read_alignment(&[first, second], &Alphabet::default())?;

    let actual: Vec<(String, String)> = records
      .into_iter()
      .map(|record| (record.seq_name, record.seq.to_string()))
      .collect();
    let expected = vec![
      ("A".to_owned(), "AC".to_owned()),
      ("B".to_owned(), "GT".to_owned()),
      ("C".to_owned(), "TT".to_owned()),
    ];
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_alignment_rejects_a_file_without_a_record() -> Result<(), Report> {
    let dir = tempdir()?;
    let first = dir.path().join("first.fasta");
    let second = dir.path().join("second.fasta");
    fs::write(&first, ">A\nAC\n")?;
    fs::write(&second, "GTAC\n")?;

    assert_error!(
      read_alignment(&[first, second.clone()], &Alphabet::default()),
      format!(
        "When reading file '{}': FASTA input is incorrectly formatted: expected at least one FASTA record starting with character '>', but none found",
        second.display()
      )
    );
    Ok(())
  }

  #[test]
  fn test_alignment_keeps_text_before_the_first_header_out_of_the_previous_file() -> Result<(), Report> {
    let dir = tempdir()?;
    let first = dir.path().join("first.fasta");
    let second = dir.path().join("second.fasta");
    fs::write(&first, ">A\nAC\n")?;
    fs::write(&second, "GT\n>B\nACGT\n")?;

    let records = read_alignment(&[first, second], &Alphabet::default())?;

    let actual: Vec<(String, String)> = records
      .into_iter()
      .map(|record| (record.seq_name, record.seq.to_string()))
      .collect();
    let expected = vec![("A".to_owned(), "AC".to_owned()), ("B".to_owned(), "ACGT".to_owned())];
    assert_eq!(expected, actual);
    Ok(())
  }
}

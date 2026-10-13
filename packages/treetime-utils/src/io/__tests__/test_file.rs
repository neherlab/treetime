#[cfg(test)]
mod tests {
  use crate::io::file::{DEFAULT_FILE_BUF_SIZE, read_file_with, write_file_with};
  use eyre::{Report, WrapErr};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::io::{Read, Write};
  use tempfile::tempdir;

  #[rstest]
  #[case::plain("out.fasta")]
  #[case::gzip("out.fasta.gz")]
  #[trace]
  fn test_file_writer_write_all_keeps_small_and_oversized_pieces_in_order(
    #[case] filename: &str,
  ) -> Result<(), Report> {
    let dir = tempdir().wrap_err("When creating a temporary directory")?;
    let path = dir.path().join(filename);
    let small = b">sample\nACGT\n";
    let oversized = vec![b'N'; DEFAULT_FILE_BUF_SIZE + 1];

    write_file_with(&path, |writer| {
      for _ in 0..1000 {
        writer.write_all(small).wrap_err("When writing a small piece")?;
      }
      writer.write_all(&oversized).wrap_err("When writing the oversized piece")?;
      writer.write_all(small).wrap_err("When writing the last piece")?;
      Ok(())
    })?;

    let actual = read_file_with(&path, |mut reader| {
      let mut bytes = Vec::new();
      reader.read_to_end(&mut bytes).wrap_err("When reading the file back")?;
      Ok(bytes)
    })?;
    let expected = [small.repeat(1000), oversized, small.to_vec()].concat();
    assert_eq!(expected, actual);
    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use crate::assert_error;
  use crate::io::compression::{CompressionType, Compressor, Decompressor};
  use eyre::Report;
  use helpers::{FailingFlush, compress, decompress};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::io::Write;

  const CONTENT: &str = "ACGTACGTNNNN-ACGT\n";

  #[rstest]
  #[case::none(CompressionType::None)]
  #[case::gzip(CompressionType::Gzip)]
  #[case::bzip2(CompressionType::Bzip2)]
  #[case::xz(CompressionType::Xz)]
  #[case::zstd(CompressionType::Zstandard)]
  fn test_compression_round_trips_after_finish(#[case] compression_type: CompressionType) -> Result<(), Report> {
    let bytes = compress(&compression_type, true)?;
    assert_eq!(CONTENT.repeat(100), decompress(&bytes, &compression_type)?);
    Ok(())
  }

  #[rstest]
  #[case::gzip(CompressionType::Gzip)]
  #[case::bzip2(CompressionType::Bzip2)]
  #[case::xz(CompressionType::Xz)]
  #[case::zstd(CompressionType::Zstandard)]
  fn test_compression_drop_finishes_an_unfinished_stream(
    #[case] compression_type: CompressionType,
  ) -> Result<(), Report> {
    let bytes = compress(&compression_type, false)?;
    assert_eq!(CONTENT.repeat(100), decompress(&bytes, &compression_type)?);
    Ok(())
  }

  #[test]
  fn test_compression_finish_returns_the_flush_error() -> Result<(), Report> {
    let mut compressor = Compressor::new(FailingFlush, &CompressionType::None)?;
    compressor.write_all(CONTENT.as_bytes())?;
    assert_error!(compressor.finish(), "While finishing compressed file: disk full");
    Ok(())
  }

  mod helpers {
    use super::*;
    use std::io::{self, Read};

    pub(super) struct FailingFlush;

    impl Write for FailingFlush {
      fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        Ok(buf.len())
      }

      fn flush(&mut self) -> io::Result<()> {
        Err(io::Error::other("disk full"))
      }
    }

    pub(super) fn compress(compression_type: &CompressionType, finish: bool) -> Result<Vec<u8>, Report> {
      let mut bytes = vec![];
      let mut compressor = Compressor::new(&mut bytes, compression_type)?;
      for _ in 0..100 {
        compressor.write_all(CONTENT.as_bytes())?;
      }
      if finish {
        compressor.finish()?;
      } else {
        drop(compressor);
      }
      Ok(bytes)
    }

    pub(super) fn decompress(bytes: &[u8], compression_type: &CompressionType) -> Result<String, Report> {
      let mut content = String::new();
      Decompressor::new(bytes, compression_type)?.read_to_string(&mut content)?;
      Ok(content)
    }
  }
}

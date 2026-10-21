#[cfg(test)]
mod tests {
  use crate::assert_error;
  use crate::io::compression::{CompressionType, Compressor, Decompressor};
  use eyre::{Report, WrapErr};
  use helpers::{FailingFlush, alignment, compress, compress_with_threads, decompress};
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
    let bytes = compress(&CONTENT.repeat(100), &compression_type, true)?;
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
    let bytes = compress(&CONTENT.repeat(100), &compression_type, false)?;
    assert_eq!(CONTENT.repeat(100), decompress(&bytes, &compression_type)?);
    Ok(())
  }

  #[rstest]
  #[case::gzip(CompressionType::Gzip)]
  #[case::xz(CompressionType::Xz)]
  #[case::zstd(CompressionType::Zstandard)]
  fn test_compression_output_does_not_depend_on_thread_count(
    #[case] compression_type: CompressionType,
  ) -> Result<(), Report> {
    let content = alignment(500, 30_000);
    let one_thread = compress_with_threads(&content, &compression_type, 1)?;
    let four_threads = compress_with_threads(&content, &compression_type, 4)?;
    assert_eq!(content, decompress(&four_threads, &compression_type)?);
    assert!(one_thread == four_threads, "output differs between 1 and 4 threads");
    Ok(())
  }

  #[test]
  fn test_compression_finish_returns_the_flush_error() -> Result<(), Report> {
    let mut compressor = Compressor::new(FailingFlush, &CompressionType::None)?;
    compressor
      .write_all(CONTENT.as_bytes())
      .wrap_err("When writing the content")?;
    assert_error!(compressor.finish(), "While finishing compressed file: disk full");
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::sync::random::get_random_number_generator;
    use itertools::Itertools;
    use parking_lot::Mutex;
    use rand::RngExt;
    use rayon::ThreadPoolBuilder;
    use std::io::{self, Read};
    use std::iter::repeat_with;
    use std::sync::Arc;

    pub(super) struct FailingFlush;

    impl Write for FailingFlush {
      fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        Ok(buf.len())
      }

      fn flush(&mut self) -> io::Result<()> {
        Err(io::Error::other("disk full"))
      }
    }

    #[derive(Clone, Default)]
    struct SharedBuffer(Arc<Mutex<Vec<u8>>>);

    impl Write for SharedBuffer {
      fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.lock().extend_from_slice(buf);
        Ok(buf.len())
      }

      fn flush(&mut self) -> io::Result<()> {
        Ok(())
      }
    }

    pub(super) fn alignment(n_sequences: usize, length: usize) -> String {
      let mut rng = get_random_number_generator(42);
      let reference = repeat_with(|| b"ACGT"[rng.random_range(0..4)])
        .take(length)
        .collect_vec();
      (0..n_sequences)
        .map(|i| {
          let mut sequence = reference.clone();
          for _ in 0..20 {
            sequence[rng.random_range(0..length)] = b"ACGT"[rng.random_range(0..4)];
          }
          format!(">seq{i}\n{}\n", String::from_utf8_lossy(&sequence))
        })
        .join("")
    }

    pub(super) fn compress(content: &str, compression_type: &CompressionType, finish: bool) -> Result<Vec<u8>, Report> {
      let buffer = SharedBuffer::default();
      let mut compressor = Compressor::new(buffer.clone(), compression_type)?;
      for line in content.split_inclusive('\n') {
        compressor
          .write_all(line.as_bytes())
          .wrap_err("When writing the content")?;
      }
      if finish {
        compressor.finish()?;
      } else {
        drop(compressor);
      }
      Ok(buffer.0.lock().clone())
    }

    pub(super) fn compress_with_threads(
      content: &str,
      compression_type: &CompressionType,
      threads: usize,
    ) -> Result<Vec<u8>, Report> {
      ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .wrap_err("When building the thread pool")?
        .install(|| compress(content, compression_type, true))
    }

    pub(super) fn decompress(bytes: &[u8], compression_type: &CompressionType) -> Result<String, Report> {
      let mut content = String::new();
      Decompressor::new(bytes, compression_type)?
        .read_to_string(&mut content)
        .wrap_err("When decompressing the content")?;
      Ok(content)
    }
  }
}

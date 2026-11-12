#[cfg(test)]
mod tests {
  use crate::runs::files::{FileChunks, RUN_DOWNLOAD_CHUNK_SIZE};
  use helpers::patterned_bytes;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::fs;
  use treetime_utils::assert_error;

  #[test]
  fn test_files_file_chunks_split_a_file_into_blocks_of_the_chunk_size() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("large.fasta");
    let content = patterned_bytes(2 * RUN_DOWNLOAD_CHUNK_SIZE + 1234);
    fs::write(&path, &content).unwrap();
    let chunks: Vec<Vec<u8>> = FileChunks::open(&path).unwrap().try_collect().unwrap();
    assert_eq!(
      (vec![RUN_DOWNLOAD_CHUNK_SIZE, RUN_DOWNLOAD_CHUNK_SIZE, 1234], true),
      (chunks.iter().map(Vec::len).collect_vec(), chunks.concat() == content)
    );
  }

  #[test]
  fn test_files_file_chunks_of_an_empty_file_are_empty() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("empty.csv");
    fs::write(&path, "").unwrap();
    assert_eq!(0, FileChunks::open(&path).unwrap().count());
  }

  #[test]
  fn test_files_file_chunks_open_fails_for_a_missing_file() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("missing.nwk");
    assert_error!(
      FileChunks::open(&path),
      format!(
        "When opening '{}': No such file or directory (os error 2)",
        path.display()
      )
    );
  }

  mod helpers {
    pub(super) fn patterned_bytes(length: usize) -> Vec<u8> {
      b"ACGT".iter().copied().cycle().take(length).collect()
    }
  }
}

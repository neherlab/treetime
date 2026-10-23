#[cfg(test)]
mod tests {
  use crate::runs::files::{FileChunks, RUN_DOWNLOAD_CHUNK_SIZE, RunZipChunks, stream_run_zip, write_run_zip};
  use helpers::{archive_entries, out_dir_with_files, patterned_bytes};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::fs;
  use std::io::Cursor;
  use treetime_utils::assert_error;
  use zip::CompressionMethod;

  #[test]
  fn test_files_write_run_zip_deflates_every_output_file() {
    let out_dir = out_dir_with_files();
    let mut buffer = Cursor::new(vec![]);
    write_run_zip(out_dir.path(), "run-1", &mut buffer).unwrap();
    assert_eq!(
      vec![
        (
          "run-1/nested/tree.nwk".to_owned(),
          CompressionMethod::Deflated,
          "(A,B);".to_owned()
        ),
        (
          "run-1/rates.json".to_owned(),
          CompressionMethod::Deflated,
          "{\"rate\":1}".to_owned()
        ),
      ],
      archive_entries(buffer.into_inner())
    );
  }

  #[test]
  fn test_files_stream_run_zip_stores_every_output_file_without_seeking() {
    let out_dir = out_dir_with_files();
    let mut buffer = vec![];
    stream_run_zip(out_dir.path(), "run-1", &mut buffer).unwrap();
    assert_eq!(
      vec![
        (
          "run-1/nested/tree.nwk".to_owned(),
          CompressionMethod::Stored,
          "(A,B);".to_owned()
        ),
        (
          "run-1/rates.json".to_owned(),
          CompressionMethod::Stored,
          "{\"rate\":1}".to_owned()
        ),
      ],
      archive_entries(buffer)
    );
  }

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

  #[test]
  fn test_files_run_zip_chunks_stay_within_the_chunk_size_for_files_larger_than_a_chunk() {
    let out_dir = out_dir_with_files();
    let large = patterned_bytes(3 * RUN_DOWNLOAD_CHUNK_SIZE + 17);
    fs::write(out_dir.path().join("large.fasta"), &large).unwrap();
    let chunks: Vec<Vec<u8>> = RunZipChunks::new(out_dir.path(), "run-1")
      .unwrap()
      .try_collect()
      .unwrap();
    let largest = chunks.iter().map(Vec::len).max().unwrap();
    let entries = archive_entries(chunks.concat());
    assert_eq!(
      (
        true,
        vec![
          ("run-1/large.fasta".to_owned(), CompressionMethod::Stored, large.len()),
          ("run-1/nested/tree.nwk".to_owned(), CompressionMethod::Stored, 6),
          ("run-1/rates.json".to_owned(), CompressionMethod::Stored, 10),
        ],
        true,
      ),
      (
        largest <= RUN_DOWNLOAD_CHUNK_SIZE,
        entries
          .iter()
          .map(|(name, method, content)| (name.clone(), *method, content.len()))
          .collect_vec(),
        entries[0].2.as_bytes() == large.as_slice(),
      )
    );
  }

  #[test]
  fn test_files_run_zip_chunks_end_after_a_file_that_cannot_be_read() {
    let out_dir = out_dir_with_files();
    let mut chunks = RunZipChunks::new(out_dir.path(), "run-1").unwrap();
    let removed = out_dir.path().join("nested/tree.nwk");
    fs::remove_file(&removed).unwrap();
    let mut items = chunks.by_ref().collect_vec();
    let last = items.pop().unwrap();
    assert_eq!(
      (true, None),
      (
        items.iter().all(Result::is_ok),
        chunks.next().map(|chunk| chunk.is_ok())
      )
    );
    assert_error!(
      last,
      format!(
        "When opening '{}': No such file or directory (os error 2)",
        removed.display()
      )
    );
  }

  #[test]
  fn test_files_run_zip_chunks_of_an_empty_folder_form_an_empty_archive() {
    let dir = tempfile::tempdir().unwrap();
    let chunks: Vec<Vec<u8>> = RunZipChunks::new(dir.path(), "run-1").unwrap().try_collect().unwrap();
    assert_eq!(
      Vec::<(String, CompressionMethod, String)>::new(),
      archive_entries(chunks.concat())
    );
  }

  mod helpers {
    use itertools::Itertools;
    use std::fs;
    use std::io::{Cursor, Read};
    use tempfile::{TempDir, tempdir};
    use zip::{CompressionMethod, ZipArchive};

    pub(super) fn out_dir_with_files() -> TempDir {
      let dir = tempdir().unwrap();
      fs::create_dir_all(dir.path().join("nested")).unwrap();
      fs::write(dir.path().join("nested/tree.nwk"), "(A,B);").unwrap();
      fs::write(dir.path().join("rates.json"), "{\"rate\":1}").unwrap();
      dir
    }

    pub(super) fn patterned_bytes(length: usize) -> Vec<u8> {
      b"ACGT".iter().copied().cycle().take(length).collect()
    }

    pub(super) fn archive_entries(bytes: Vec<u8>) -> Vec<(String, CompressionMethod, String)> {
      let mut archive = ZipArchive::new(Cursor::new(bytes)).unwrap();
      (0..archive.len())
        .map(|index| {
          let mut file = archive.by_index(index).unwrap();
          let mut content = String::new();
          file.read_to_string(&mut content).unwrap();
          (file.name().unwrap().into_owned(), file.compression(), content)
        })
        .collect_vec()
    }
  }
}

#[cfg(test)]
mod tests {
  use crate::runs::files::{stream_run_zip, write_run_zip};
  use helpers::{archive_entries, out_dir_with_files};
  use pretty_assertions::assert_eq;
  use std::io::Cursor;
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

    pub(super) fn archive_entries(bytes: Vec<u8>) -> Vec<(String, CompressionMethod, String)> {
      let mut archive = ZipArchive::new(Cursor::new(bytes)).unwrap();
      (0..archive.len())
        .map(|index| {
          let mut file = archive.by_index(index).unwrap();
          let mut content = String::new();
          file.read_to_string(&mut content).unwrap();
          (file.name().to_owned(), file.compression(), content)
        })
        .collect_vec()
    }
  }
}

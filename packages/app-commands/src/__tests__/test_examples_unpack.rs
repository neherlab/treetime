#[cfg(test)]
mod tests {
  use crate::examples_unpack::UnpackLimits;
  use helpers::{Entry, archive, files_in, links_in, rename_entry, unpack};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::fs;
  #[cfg(unix)]
  use std::os::unix::fs::PermissionsExt;
  use tempfile::tempdir;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::parent_folder(     archive(&[Entry::File("../escaped.txt", "x")]),     "the example archive is invalid: entry '../escaped.txt' has a path outside the folder")]
  #[case::absolute(          archive(&[Entry::File("/abs/escaped.txt", "x")]),   "the example archive is invalid: entry '/abs/escaped.txt' has a path outside the folder")]
  #[case::backslash(         archive(&[Entry::File("a\\b.txt", "x")]),           "the example archive is invalid: entry 'a\\b.txt' has a path component with the character '\\\\'")]
  #[case::drive(             archive(&[Entry::File("C:/x.txt", "x")]),           "the example archive is invalid: entry 'C:/x.txt' has a path component with the character ':'")]
  #[case::device_with_dot(   archive(&[Entry::File("zika/CON.txt", "x")]),       "the example archive is invalid: entry 'zika/CON.txt' has a path component that is a device name on Windows")]
  #[case::device_superscript(archive(&[Entry::File("zika/COM¹.txt", "x")]),      "the example archive is invalid: entry 'zika/COM¹.txt' has a path component that is a device name on Windows")]
  #[case::device_space(      archive(&[Entry::File("zika/nul .txt", "x")]),      "the example archive is invalid: entry 'zika/nul .txt' has a path component that is a device name on Windows")]
  #[case::trailing_dot(      archive(&[Entry::File("zika/name.", "x")]),         "the example archive is invalid: entry 'zika/name.' has a path component that ends with '.' or a space")]
  #[case::symbolic_link(     archive(&[Entry::Link("zika/link", "../../outside")]), "the example archive is invalid: entry 'zika/link' is a symbolic link")]
  #[case::fifo(              archive(&[Entry::Mode("zika/fifo", "", 0o010_644)]), "the example archive is invalid: entry 'zika/fifo' is not a regular file or a folder")]
  #[case::folder_content(    archive(&[Entry::Mode("zika", "x", 0o040_755)]),    "the example archive is invalid: entry 'zika' is a folder with content")]
  #[case::duplicate_name(    rename_entry(archive(&[Entry::File("a.txt", "1"), Entry::File("b.txt", "2")]), "b.txt", "a.txt"), "the example archive is invalid: invalid Zip archive: header inconsistency detected: duplicate file names in central directory")]
  #[trace]
  fn test_examples_unpack_refuses_an_unsafe_archive(#[case] bytes: Vec<u8>, #[case] expected: &str) {
    let root = tempdir().unwrap();
    let folder = root.path().join("examples");
    let result = unpack(bytes, &folder, UnpackLimits::default());
    assert_error!(result, expected);
    assert_eq!(
      (Vec::<String>::new(), Vec::<String>::new()),
      (
        fs::read_dir(root.path())
          .unwrap()
          .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
          .filter(|name| name != "examples")
          .collect::<Vec<_>>(),
        links_in(&folder)
      )
    );
  }

  #[test]
  fn test_examples_unpack_refuses_more_entries_than_the_limit() {
    let root = tempdir().unwrap();
    let bytes = archive(&[Entry::File("a.txt", "1"), Entry::File("b.txt", "2")]);
    let limits = UnpackLimits {
      entries: 1,
      ..UnpackLimits::default()
    };
    assert_error!(
      unpack(bytes, &root.path().join("examples"), limits),
      "the example archive has 2 entries, more than the limit of 1"
    );
  }

  #[test]
  fn test_examples_unpack_refuses_more_bytes_than_the_limit() {
    let root = tempdir().unwrap();
    let bytes = archive(&[Entry::File("a.txt", "1234"), Entry::File("b.txt", "5678")]);
    let limits = UnpackLimits {
      bytes: 6,
      ..UnpackLimits::default()
    };
    assert_error!(
      unpack(bytes, &root.path().join("examples"), limits),
      "the example archive unpacks to more than the limit of 6 bytes (at entry 'b.txt')"
    );
  }

  #[test]
  fn test_examples_unpack_refuses_two_names_of_the_same_path() {
    let root = tempdir().unwrap();
    let bytes = archive(&[Entry::File("a/b.txt", "1"), Entry::File("a/./b.txt", "2")]);
    let Err(error) = unpack(bytes, &root.path().join("examples"), UnpackLimits::default()) else {
      panic!("expected Err, got Ok")
    };
    assert_eq!("When unpacking entry 'a/./b.txt'", error.to_string());
  }

  #[cfg(unix)]
  #[test]
  fn test_examples_unpack_drops_the_setuid_bit() {
    let root = tempdir().unwrap();
    let folder = root.path().join("examples");
    unpack(
      archive(&[Entry::Mode("zika/run", "content", 0o104_755)]),
      &folder,
      UnpackLimits::default(),
    )
    .unwrap();
    let mode = fs::metadata(folder.join("zika/run")).unwrap().permissions().mode();
    assert_eq!(
      (vec![("zika/run".to_owned(), "content".to_owned())], 0),
      (files_in(&folder), mode & 0o7000)
    );
  }

  #[test]
  fn test_examples_unpack_writes_a_nested_archive() {
    let root = tempdir().unwrap();
    let folder = root.path().join("examples");
    unpack(
      archive(&[
        Entry::Folder("flu/"),
        Entry::File("flu/h3n2/20/tree.nwk", "(A:1,B:1);\n"),
        Entry::File("./zika/20/metadata.tsv", "strain\tdate\n"),
        Entry::File("index.json", "{}\n"),
      ]),
      &folder,
      UnpackLimits::default(),
    )
    .unwrap();
    assert_eq!(
      vec![
        ("flu/h3n2/20/tree.nwk".to_owned(), "(A:1,B:1);\n".to_owned()),
        ("index.json".to_owned(), "{}\n".to_owned()),
        ("zika/20/metadata.tsv".to_owned(), "strain\tdate\n".to_owned()),
      ],
      files_in(&folder)
    );
  }

  mod helpers {
    use crate::examples_unpack::{UnpackLimits, unpack_examples};
    use eyre::Report;
    use std::fs;
    use std::io::{Cursor, Write};
    use std::path::Path;
    use zip::ZipArchive;
    use zip::ZipWriter;
    use zip::write::SimpleFileOptions;

    pub(super) enum Entry {
      File(&'static str, &'static str),
      Mode(&'static str, &'static str, u32),
      Link(&'static str, &'static str),
      Folder(&'static str),
    }

    pub(super) fn archive(entries: &[Entry]) -> Vec<u8> {
      let mut writer = ZipWriter::new(Cursor::new(Vec::new()));
      for entry in entries {
        match entry {
          Entry::File(name, content) => {
            writer.start_file(*name, SimpleFileOptions::default()).unwrap();
            writer.write_all(content.as_bytes()).unwrap();
          },
          Entry::Mode(name, content, mode) => {
            let options = SimpleFileOptions::default().external_attributes(mode << 16);
            writer.start_file(*name, options).unwrap();
            writer.write_all(content.as_bytes()).unwrap();
          },
          Entry::Link(name, target) => {
            writer
              .add_symlink(*name, *target, SimpleFileOptions::default())
              .unwrap();
          },
          Entry::Folder(name) => {
            writer.add_directory(*name, SimpleFileOptions::default()).unwrap();
          },
        }
      }
      writer.finish().unwrap().into_inner()
    }

    pub(super) fn rename_entry(bytes: Vec<u8>, from: &str, to: &str) -> Vec<u8> {
      assert_eq!(from.len(), to.len());
      let mut bytes = bytes;
      let positions = bytes
        .windows(from.len())
        .enumerate()
        .filter(|(_, window)| *window == from.as_bytes())
        .map(|(position, _)| position)
        .collect::<Vec<_>>();
      assert_eq!(2, positions.len());
      for position in positions {
        bytes[position..position + to.len()].copy_from_slice(to.as_bytes());
      }
      bytes
    }

    pub(super) fn unpack(bytes: Vec<u8>, folder: &Path, limits: UnpackLimits) -> Result<(), Report> {
      let mut zip = ZipArchive::new(Cursor::new(bytes))?;
      unpack_examples(&mut zip, folder, limits)
    }

    pub(super) fn files_in(folder: &Path) -> Vec<(String, String)> {
      let mut files = vec![];
      collect(folder, folder, &mut files);
      files.sort();
      files
    }

    pub(super) fn links_in(folder: &Path) -> Vec<String> {
      let mut links = vec![];
      collect_links(folder, folder, &mut links);
      links
    }

    fn collect(root: &Path, folder: &Path, files: &mut Vec<(String, String)>) {
      for entry in fs::read_dir(folder).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
          collect(root, &path, files);
        } else {
          let name = path.strip_prefix(root).unwrap().to_string_lossy().into_owned();
          files.push((name, fs::read_to_string(&path).unwrap()));
        }
      }
    }

    fn collect_links(root: &Path, folder: &Path, links: &mut Vec<String>) {
      if !folder.exists() {
        return;
      }
      for entry in fs::read_dir(folder).unwrap() {
        let path = entry.unwrap().path();
        let file_type = fs::symlink_metadata(&path).unwrap().file_type();
        if file_type.is_symlink() {
          links.push(path.strip_prefix(root).unwrap().to_string_lossy().into_owned());
        } else if file_type.is_dir() {
          collect_links(root, &path, links);
        }
      }
    }
  }
}

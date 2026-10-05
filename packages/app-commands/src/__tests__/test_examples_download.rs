#[cfg(test)]
mod tests {
  use crate::examples_download::{DownloadProgress, download_examples};
  use helpers::{archive, files_in, serve};
  use pretty_assertions::assert_eq;
  use std::fs;
  use tempfile::tempdir;
  use treetime_utils::assert_error;

  #[test]
  fn test_examples_download_unpacks_the_archive_into_a_missing_folder() {
    let bytes = archive(&[
      ("zika/20/tree.nwk", "(A:1,B:1);\n"),
      ("zika/20/ancestral.yaml", "tree: tree.nwk\n"),
    ]);
    let server = serve(bytes);
    let root = tempdir().unwrap();
    let target = root.path().join("examples");
    download_examples(&server.url("examples.zip"), &target, |_| {}).unwrap();
    assert_eq!(
      (
        vec![
          ("zika/20/ancestral.yaml".to_owned(), "tree: tree.nwk\n".to_owned()),
          ("zika/20/tree.nwk".to_owned(), "(A:1,B:1);\n".to_owned()),
        ],
        vec!["examples".to_owned()]
      ),
      (files_in(&target), files_in_top(root.path()))
    );
  }

  #[test]
  fn test_examples_download_replaces_an_empty_folder() {
    let server = serve(archive(&[("ebola/20/tree.nwk", "(A:1);\n")]));
    let root = tempdir().unwrap();
    let target = root.path().join("examples");
    fs::create_dir_all(&target).unwrap();
    download_examples(&server.url("examples.zip"), &target, |_| {}).unwrap();
    assert_eq!(
      vec![("ebola/20/tree.nwk".to_owned(), "(A:1);\n".to_owned())],
      files_in(&target)
    );
  }

  #[test]
  fn test_examples_download_treats_a_folder_with_finder_metadata_only_as_empty() {
    let server = serve(archive(&[("ebola/20/tree.nwk", "(A:1);\n")]));
    let root = tempdir().unwrap();
    let target = root.path().join("examples");
    fs::create_dir_all(&target).unwrap();
    fs::write(target.join(".DS_Store"), "finder").unwrap();
    download_examples(&server.url("examples.zip"), &target, |_| {}).unwrap();
    assert_eq!(
      vec![("ebola/20/tree.nwk".to_owned(), "(A:1);\n".to_owned())],
      files_in(&target)
    );
  }

  #[test]
  fn test_examples_download_refuses_a_folder_that_is_not_empty() {
    let server = serve(archive(&[("ebola/20/tree.nwk", "(A:1);\n")]));
    let root = tempdir().unwrap();
    let target = root.path().join("examples");
    fs::create_dir_all(&target).unwrap();
    fs::write(target.join("mine.nwk"), "(X:1);\n").unwrap();
    assert_error!(
      download_examples(&server.url("examples.zip"), &target, |_| {}),
      format!(
        "the folder '{}' is not empty; the example datasets go into an empty or missing folder",
        target.display()
      )
    );
    assert_eq!(vec![("mine.nwk".to_owned(), "(X:1);\n".to_owned())], files_in(&target));
  }

  #[test]
  fn test_examples_download_leaves_no_temporaries_after_a_failed_transfer() {
    let server = serve(archive(&[("ebola/20/tree.nwk", "(A:1);\n")]));
    let root = tempdir().unwrap();
    let target = root.path().join("examples");
    let url = server.url("missing.zip");
    assert_error!(
      download_examples(&url, &target, |_| {}),
      format!("When downloading '{url}': HTTP status client error (404 Not Found) for url ({url})")
    );
    assert_eq!(Vec::<String>::new(), files_in_top(root.path()));
  }

  #[test]
  fn test_examples_download_leaves_no_temporaries_after_an_invalid_archive() {
    let server = serve(b"not a zip archive".to_vec());
    let root = tempdir().unwrap();
    let target = root.path().join("examples");
    let url = server.url("examples.zip");
    assert_error!(
      download_examples(&url, &target, |_| {}),
      format!("When unpacking the example datasets from '{url}': invalid Zip archive: Could not find EOCD")
    );
    assert_eq!(Vec::<String>::new(), files_in_top(root.path()));
  }

  #[test]
  fn test_examples_download_removes_stale_temporaries_of_an_interrupted_download() {
    let server = serve(archive(&[("ebola/20/tree.nwk", "(A:1);\n")]));
    let root = tempdir().unwrap();
    fs::write(root.path().join(".treetime-examples-abc.zip"), "partial").unwrap();
    fs::create_dir_all(root.path().join(".treetime-examples-def/zika")).unwrap();
    fs::write(root.path().join("unrelated.txt"), "kept").unwrap();
    download_examples(&server.url("examples.zip"), &root.path().join("examples"), |_| {}).unwrap();
    assert_eq!(
      vec!["examples".to_owned(), "unrelated.txt".to_owned()],
      files_in_top(root.path())
    );
  }

  #[test]
  fn test_examples_download_refuses_an_entry_that_leaves_the_folder() {
    let server = serve(archive(&[("zika/tree.nwk", "(A:1);\n"), ("../escaped.txt", "outside")]));
    let root = tempdir().unwrap();
    let target = root.path().join("nested").join("examples");
    let result = download_examples(&server.url("examples.zip"), &target, |_| {});
    assert_eq!(
      (true, false, Vec::<String>::new()),
      (
        result.is_err(),
        root.path().join("nested").join("escaped.txt").exists(),
        files_in_top(&root.path().join("nested"))
      )
    );
  }

  #[test]
  fn test_examples_download_keeps_an_absolute_entry_inside_the_folder() {
    let outside = tempdir().unwrap();
    let absolute = outside.path().join("escaped.txt");
    let server = serve(archive(&[(absolute.to_str().unwrap(), "outside")]));
    let root = tempdir().unwrap();
    let target = root.path().join("examples");
    download_examples(&server.url("examples.zip"), &target, |_| {}).unwrap();
    let inside = absolute.strip_prefix("/").unwrap().to_string_lossy().into_owned();
    assert_eq!(
      (false, vec![(inside, "outside".to_owned())]),
      (absolute.exists(), files_in(&target))
    );
  }

  #[test]
  fn test_examples_download_reports_received_bytes_up_to_the_size_of_the_archive() {
    let bytes = archive(&[("zika/20/tree.nwk", &"(A:1,B:1);\n".repeat(20_000))]);
    let size = u64::try_from(bytes.len()).unwrap();
    let server = serve(bytes);
    let root = tempdir().unwrap();
    let mut reports = vec![];
    download_examples(&server.url("examples.zip"), &root.path().join("examples"), |progress| {
      reports.push(progress);
    })
    .unwrap();
    let increasing = reports
      .iter()
      .zip(reports.iter().skip(1))
      .all(|(earlier, later)| earlier.received <= later.received);
    assert_eq!(
      (
        Some(&DownloadProgress {
          received: 0,
          total: Some(size)
        }),
        Some(&DownloadProgress {
          received: size,
          total: Some(size)
        }),
        true
      ),
      (reports.first(), reports.last(), increasing)
    );
  }

  fn files_in_top(folder: &std::path::Path) -> Vec<String> {
    let mut names = fs::read_dir(folder)
      .unwrap()
      .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
      .collect::<Vec<_>>();
    names.sort();
    names
  }

  mod helpers {
    use axum::Router;
    use axum::routing::get;
    use std::fs;
    use std::io::{Cursor, Write};
    use std::net::SocketAddr;
    use std::path::Path;
    use std::sync::mpsc;
    use std::thread;
    use zip::ZipWriter;
    use zip::write::SimpleFileOptions;

    pub(super) struct Server {
      address: SocketAddr,
    }

    impl Server {
      pub(super) fn url(&self, name: &str) -> String {
        format!("http://{}/{name}", self.address)
      }
    }

    pub(super) fn serve(bytes: Vec<u8>) -> Server {
      let (send, receive) = mpsc::sync_channel(1);
      thread::spawn(move || {
        let runtime = tokio::runtime::Builder::new_current_thread()
          .enable_all()
          .build()
          .unwrap();
        runtime.block_on(async move {
          let listener = tokio::net::TcpListener::bind("127.0.0.1:0").await.unwrap();
          send.send(listener.local_addr().unwrap()).unwrap();
          let router = Router::new().route("/examples.zip", get(move || async move { bytes }));
          axum::serve(listener, router).await.unwrap();
        });
      });
      Server {
        address: receive.recv().unwrap(),
      }
    }

    pub(super) fn archive(entries: &[(&str, &str)]) -> Vec<u8> {
      let mut writer = ZipWriter::new(Cursor::new(Vec::new()));
      for (name, content) in entries {
        writer.start_file(*name, SimpleFileOptions::default()).unwrap();
        writer.write_all(content.as_bytes()).unwrap();
      }
      writer.finish().unwrap().into_inner()
    }

    pub(super) fn files_in(folder: &Path) -> Vec<(String, String)> {
      let mut files = vec![];
      collect(folder, folder, &mut files);
      files.sort();
      files
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
  }
}

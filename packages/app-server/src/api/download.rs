use crate::api::response::{FileContent, ZipAttachment};
use crate::error::AppError;
use app_commands::runs::files::stream_run_zip;
use axum::body::Body;
use axum::extract::Request;
use axum::http::{HeaderValue, header};
use eyre::Report;
use log::warn;
use percent_encoding::{NON_ALPHANUMERIC, utf8_percent_encode};
use std::path::{Path, PathBuf};
use tokio_util::io::{ReaderStream, SyncIoBridge};
use tower::ServiceExt as _;
use tower_http::services::ServeFile;

const ARCHIVE_BUFFER_SIZE: usize = 64 * 1024;

const CONTENT_TYPES: &[(&str, &str)] = &[
  ("json", "application/json"),
  ("svg", "image/svg+xml"),
  ("png", "image/png"),
  ("zip", "application/zip"),
  ("xz", "application/x-xz"),
  ("gz", "application/gzip"),
  ("zst", "application/zstd"),
  ("nwk", "text/plain; charset=utf-8"),
  ("nexus", "text/plain; charset=utf-8"),
  ("csv", "text/csv; charset=utf-8"),
  ("tsv", "text/tab-separated-values; charset=utf-8"),
  ("fasta", "text/plain; charset=utf-8"),
  ("dot", "text/plain; charset=utf-8"),
  ("jsonl", "application/jsonl"),
];

pub(crate) async fn serve_run_file(path: PathBuf, request: Request) -> Result<FileContent, AppError> {
  let content_type = content_type(&path);
  let disposition = attachment(&file_name(&path))?;
  let mut response = ServeFile::new(&path).oneshot(request).await?.map(Body::new);
  if response.status().is_success() {
    let headers = response.headers_mut();
    headers.insert(header::CONTENT_TYPE, content_type);
    headers.insert(header::CONTENT_DISPOSITION, disposition);
  }
  Ok(FileContent(response))
}

pub(crate) fn stream_run_archive(out_dir: PathBuf, name: &str) -> Result<ZipAttachment, AppError> {
  let (reader, writer) = tokio::io::duplex(ARCHIVE_BUFFER_SIZE);
  let bridge = SyncIoBridge::new(writer);
  let folder = name.to_owned();
  tokio::task::spawn_blocking(move || {
    if let Err(report) = stream_run_zip(&out_dir, &folder, bridge) {
      warn!("The archive of run '{folder}' ended early: {report:#}");
    }
  });
  Ok(ZipAttachment {
    disposition: attachment(&format!("treetime-{name}.zip"))?,
    body: Body::from_stream(ReaderStream::new(reader)),
  })
}

fn content_type(path: &Path) -> HeaderValue {
  let extension = path.extension().and_then(|extension| extension.to_str()).unwrap_or("");
  let content_type = CONTENT_TYPES
    .iter()
    .find(|(known, _)| *known == extension)
    .map_or("application/octet-stream", |(_, content_type)| *content_type);
  HeaderValue::from_static(content_type)
}

fn file_name(path: &Path) -> String {
  path
    .file_name()
    .map_or_else(String::new, |name| name.to_string_lossy().into_owned())
}

fn attachment(name: &str) -> Result<HeaderValue, Report> {
  let fallback: String = name
    .chars()
    .map(|c| {
      if c.is_ascii_graphic() && !matches!(c, '"' | '\\' | '%') || c == ' ' {
        c
      } else {
        '_'
      }
    })
    .collect();
  let encoded = utf8_percent_encode(name, NON_ALPHANUMERIC);
  Ok(HeaderValue::from_str(&format!(
    "attachment; filename=\"{fallback}\"; filename*=UTF-8''{encoded}"
  ))?)
}

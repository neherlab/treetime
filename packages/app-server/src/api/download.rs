use crate::api::response::{FileContent, ZipAttachment};
use crate::error::AppError;
use axum::body::Body;
use axum::extract::Request;
use axum::http::{HeaderValue, header};
use axum::response::Response;
use eyre::Report;
use log::warn;
use percent_encoding::{NON_ALPHANUMERIC, utf8_percent_encode};
use std::path::{Path, PathBuf};

#[cfg(not(target_family = "wasm"))]
use app_backend::runs::files::stream_run_zip;
#[cfg(target_family = "wasm")]
use app_backend::runs::files::{FileChunks, RunZipChunks};
#[cfg(target_family = "wasm")]
use std::convert::Infallible;
#[cfg(target_family = "wasm")]
use std::future::{Ready, ready};
#[cfg(target_family = "wasm")]
use tokio_stream::iter;
#[cfg(not(target_family = "wasm"))]
use tokio_util::io::{ReaderStream, SyncIoBridge};
#[cfg(not(target_family = "wasm"))]
use tower::ServiceExt as _;
#[cfg(not(target_family = "wasm"))]
use tower_http::services::ServeFile;

#[cfg(not(target_family = "wasm"))]
const ARCHIVE_BUFFER_SIZE: usize = 64 * 1024;

const CONTENT_TYPES: &[(&str, &str)] = &[
  ("json", "application/json"),
  ("svg", "image/svg+xml"),
  ("png", "image/png"),
  ("pdf", "application/pdf"),
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
  let mut response = file_response(&path, request).await?;
  if response.status().is_success() {
    let headers = response.headers_mut();
    headers.insert(header::CONTENT_TYPE, content_type);
    headers.insert(header::CONTENT_DISPOSITION, disposition);
  }
  Ok(FileContent(response))
}

pub(crate) fn stream_run_archive(out_dir: &Path, name: &str) -> Result<ZipAttachment, AppError> {
  Ok(ZipAttachment {
    disposition: attachment(&format!("treetime-{name}.zip"))?,
    body: archive_body(out_dir, name),
  })
}

#[cfg(not(target_family = "wasm"))]
async fn file_response(path: &Path, request: Request) -> Result<Response, AppError> {
  Ok(ServeFile::new(path).oneshot(request).await?.map(Body::new))
}

#[cfg(target_family = "wasm")]
fn file_response(path: &Path, _request: Request) -> Ready<Result<Response, AppError>> {
  ready(
    FileChunks::open(path)
      .map(|chunks| Response::new(Body::from_stream(iter(chunks))))
      .map_err(AppError::from),
  )
}

#[cfg(not(target_family = "wasm"))]
#[cfg_attr(
  dylint_lib = "custom",
  expect(
    spawn_handle_dropped,
    reason = "the archive task logs its error and ends the body stream, which the client sees as a cut download"
  )
)]
fn archive_body(out_dir: &Path, folder: &str) -> Body {
  let (reader, writer) = tokio::io::duplex(ARCHIVE_BUFFER_SIZE);
  let bridge = SyncIoBridge::new(writer);
  let out_dir = out_dir.to_owned();
  let folder = folder.to_owned();
  tokio::task::spawn_blocking(move || {
    if let Err(report) = stream_run_zip(&out_dir, &folder, bridge) {
      warn_archive_ended(&folder, &report);
    }
  });
  Body::from_stream(ReaderStream::new(reader))
}

#[cfg(target_family = "wasm")]
fn archive_body(out_dir: &Path, folder: &str) -> Body {
  let chunks = match RunZipChunks::new(out_dir, folder) {
    Ok(chunks) => chunks,
    Err(report) => {
      warn_archive_ended(folder, &report);
      return Body::empty();
    },
  };
  let folder = folder.to_owned();
  Body::from_stream(iter(chunks.map_while(move |chunk| match chunk {
    Ok(bytes) => Some(Ok::<_, Infallible>(bytes)),
    Err(report) => {
      warn_archive_ended(&folder, &report);
      None
    },
  })))
}

fn warn_archive_ended(folder: &str, report: &Report) {
  warn!("The archive of run '{folder}' ended early: {report:#}");
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

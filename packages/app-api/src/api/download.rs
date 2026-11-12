use crate::api::response::{FileContent, ZipAttachment};
use crate::error::AppError;
use axum::body::Body;
use axum::extract::Request;
use axum::http::{HeaderValue, header};
use axum::response::Response;
use eyre::Report;
use log::warn;
use percent_encoding::{NON_ALPHANUMERIC, utf8_percent_encode};
use std::path::Path;

use crate::api::blocking::run_blocking;
use app_backend::runs::downloads::FileDownload;
#[cfg(target_family = "wasm")]
use app_backend::runs::files::FileChunks;
use app_backend::runs::package::RunPackage;
use parking_lot::Mutex;
#[cfg(target_family = "wasm")]
use std::future::{Ready, ready};
use std::iter::repeat;
use std::sync::Arc;
use tokio_stream::StreamExt as _;
use tokio_stream::iter;
#[cfg(not(target_family = "wasm"))]
use tower::ServiceExt as _;
#[cfg(not(target_family = "wasm"))]
use tower_http::services::ServeFile;

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

pub(crate) async fn serve_run_file(
  FileDownload { path, name }: FileDownload,
  request: Request,
) -> Result<FileContent, AppError> {
  let content_type = content_type(&path);
  let disposition = attachment(&name)?;
  let mut response = file_response(&path, request).await?;
  if response.status().is_success() {
    let headers = response.headers_mut();
    headers.insert(header::CONTENT_TYPE, content_type);
    headers.insert(header::CONTENT_DISPOSITION, disposition);
  }
  Ok(FileContent(response))
}

pub(crate) fn stream_run_package(package: RunPackage) -> Result<ZipAttachment, AppError> {
  Ok(ZipAttachment {
    disposition: attachment(package.file_name())?,
    body: package_body(package),
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

fn package_body(package: RunPackage) -> Body {
  let folder = package.folder().to_owned();
  let chunks = Arc::new(Mutex::new(package.chunks()));
  Body::from_stream(
    iter(repeat(()))
      .then(move |()| {
        let chunks = Arc::clone(&chunks);
        run_blocking(move || chunks.lock().next().transpose())
      })
      .map_while(Result::transpose)
      .map(move |chunk| {
        chunk.map_err(|error| {
          let report = error.into_report();
          warn_package_ended(&folder, &report);
          report
        })
      }),
  )
}

fn warn_package_ended(folder: &str, report: &Report) {
  warn!("The run package '{folder}' ended early: {report:#}");
}

fn content_type(path: &Path) -> HeaderValue {
  let extension = path.extension().and_then(|extension| extension.to_str()).unwrap_or("");
  let content_type = CONTENT_TYPES
    .iter()
    .find(|(known, _)| *known == extension)
    .map_or("application/octet-stream", |(_, content_type)| *content_type);
  HeaderValue::from_static(content_type)
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

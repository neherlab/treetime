use crate::error::AppError;
use app_backend::runs::manager::{InputUpload, UploadedInput};
use app_backend::runs::package_import::PackageImport;
use app_backend::runs::record::RunRecord;
use axum::body::Body;
use eyre::Report;
use tokio_stream::StreamExt as _;

#[cfg(target_family = "wasm")]
use std::future::{Ready, ready};

pub(crate) async fn receive_body<S: BodySink>(body: Body, mut sink: S) -> Result<S::Output, AppError> {
  let mut chunks = body.into_data_stream();
  while let Some(chunk) = chunks.next().await {
    let chunk = chunk?;
    sink = run_blocking(move || {
      sink.write(&chunk)?;
      Ok(sink)
    })
    .await?;
  }
  run_blocking(move || sink.finish()).await
}

pub(crate) trait BodySink: Send + 'static {
  type Output: Send + 'static;

  fn write(&mut self, chunk: &[u8]) -> Result<(), Report>;

  fn finish(self) -> Result<Self::Output, Report>;
}

impl BodySink for InputUpload {
  type Output = UploadedInput;

  fn write(&mut self, chunk: &[u8]) -> Result<(), Report> {
    Self::write(self, chunk)
  }

  fn finish(self) -> Result<UploadedInput, Report> {
    Self::finish(self)
  }
}

impl BodySink for PackageImport {
  type Output = RunRecord;

  fn write(&mut self, chunk: &[u8]) -> Result<(), Report> {
    Self::write(self, chunk)
  }

  fn finish(self) -> Result<RunRecord, Report> {
    Self::finish(self)
  }
}

#[cfg(not(target_family = "wasm"))]
pub(crate) async fn run_blocking<T, F>(work: F) -> Result<T, AppError>
where
  T: Send + 'static,
  F: FnOnce() -> Result<T, Report> + Send + 'static,
{
  Ok(tokio::task::spawn_blocking(work).await??)
}

#[cfg(target_family = "wasm")]
pub(crate) fn run_blocking<T, F>(work: F) -> Ready<Result<T, AppError>>
where
  T: Send + 'static,
  F: FnOnce() -> Result<T, Report> + Send + 'static,
{
  ready(work().map_err(AppError::from))
}

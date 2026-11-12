use crate::error::AppError;
use app_backend::runs::manager::{InputUpload, UploadedInput};
use axum::body::Body;
use eyre::Report;
use std::io::Read;
use tokio_stream::StreamExt as _;

#[cfg(target_family = "wasm")]
use axum::body::to_bytes;
#[cfg(target_family = "wasm")]
use std::future::{Ready, ready};
#[cfg(not(target_family = "wasm"))]
use std::io;
#[cfg(not(target_family = "wasm"))]
use tokio_util::io::{StreamReader, SyncIoBridge};

pub(crate) async fn upload_body(body: Body, mut upload: InputUpload) -> Result<UploadedInput, AppError> {
  let mut chunks = body.into_data_stream();
  while let Some(chunk) = chunks.next().await {
    let chunk = chunk?;
    upload = run_blocking(move || {
      upload.write(&chunk)?;
      Ok(upload)
    })
    .await?;
  }
  run_blocking(move || upload.finish()).await
}

#[cfg(not(target_family = "wasm"))]
pub(crate) async fn read_body_blocking<T, F>(body: Body, work: F) -> Result<T, AppError>
where
  T: Send + 'static,
  F: FnOnce(&mut dyn Read) -> Result<T, Report> + Send + 'static,
{
  let stream = body.into_data_stream().map(|chunk| chunk.map_err(io::Error::other));
  let mut reader = SyncIoBridge::new(StreamReader::new(stream));
  run_blocking(move || work(&mut reader)).await
}

#[cfg(target_family = "wasm")]
pub(crate) async fn read_body_blocking<T, F>(body: Body, work: F) -> Result<T, AppError>
where
  T: Send + 'static,
  F: FnOnce(&mut dyn Read) -> Result<T, Report> + Send + 'static,
{
  let bytes = to_bytes(body, usize::MAX).await?;
  run_blocking(move || work(&mut bytes.as_ref())).await
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

use crate::error::AppError;
use axum::body::Body;
use eyre::Report;
use std::io::Read;

#[cfg(target_family = "wasm")]
use axum::body::to_bytes;
#[cfg(target_family = "wasm")]
use std::future::{Ready, ready};
#[cfg(not(target_family = "wasm"))]
use std::io;
#[cfg(not(target_family = "wasm"))]
use tokio_stream::StreamExt as _;
#[cfg(not(target_family = "wasm"))]
use tokio_util::io::{StreamReader, SyncIoBridge};

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

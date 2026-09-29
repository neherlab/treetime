use aide::OperationOutput;
use aide::generate::GenContext;
use aide::openapi::{self, Operation};
use app_commands::runs::record::{RunRecord, RunStatus};
use axum::http::{HeaderMap, HeaderValue, StatusCode, header};
use axum::response::{IntoResponse, Response};
use eyre::Report;
use headers::{ETag, HeaderMapExt, IfNoneMatch};
use itertools::Itertools;
use std::hash::{DefaultHasher, Hash, Hasher};
use treetime_utils::io::json::{JsonPretty, json_write_str};

const REVALIDATE: &str = "private, no-cache";

pub(crate) fn finished_runs_tag(instance: u64, records: &[RunRecord]) -> Result<Option<ETag>, Report> {
  if !records.iter().all(|record| record.status == RunStatus::Ok) {
    return Ok(None);
  }
  let mut hasher = DefaultHasher::new();
  instance.hash(&mut hasher);
  records
    .iter()
    .map(|record| json_write_str(record, JsonPretty(false)))
    .try_collect::<_, Vec<_>, _>()?
    .hash(&mut hasher);
  Ok(Some(format!("W/\"{:016x}\"", hasher.finish()).parse::<ETag>()?))
}

pub(crate) fn is_unchanged(headers: &HeaderMap, tag: Option<&ETag>) -> bool {
  match (headers.typed_get::<IfNoneMatch>(), tag) {
    (Some(if_none_match), Some(tag)) => !if_none_match.precondition_passes(tag),
    _ => false,
  }
}

pub(crate) enum Revalidated<T> {
  Unchanged(ETag),
  Fresh(Option<ETag>, T),
}

impl<T: IntoResponse> IntoResponse for Revalidated<T> {
  fn into_response(self) -> Response {
    let (tag, mut response) = match self {
      Self::Unchanged(tag) => (tag, StatusCode::NOT_MODIFIED.into_response()),
      Self::Fresh(None, body) => return body.into_response(),
      Self::Fresh(Some(tag), body) => (tag, body.into_response()),
    };
    let headers = response.headers_mut();
    headers.typed_insert(tag);
    headers.insert(header::CACHE_CONTROL, HeaderValue::from_static(REVALIDATE));
    response
  }
}

impl<T: OperationOutput> OperationOutput for Revalidated<T> {
  type Inner = T::Inner;

  fn operation_response(ctx: &mut GenContext, operation: &mut Operation) -> Option<openapi::Response> {
    T::operation_response(ctx, operation)
  }

  fn inferred_responses(
    ctx: &mut GenContext,
    operation: &mut Operation,
  ) -> Vec<(Option<openapi::StatusCode>, openapi::Response)> {
    let mut responses = T::inferred_responses(ctx, operation);
    responses.push((
      Some(openapi::StatusCode::Code(304)),
      openapi::Response {
        description: "The answer is unchanged since the response that carried the `If-None-Match` tag".to_owned(),
        ..openapi::Response::default()
      },
    ));
    responses
  }
}

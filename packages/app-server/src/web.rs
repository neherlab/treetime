use crate::error::plain_error;
use crate::routes::{NO_STORE, not_found};
use app_commands::bridge::error::ErrorCode;
use axum::Router;
use axum::body::Body;
use axum::extract::{Request, State};
use axum::handler::Handler as _;
use axum::http::uri::Authority;
use axum::http::{self, HeaderName, HeaderValue, Method, StatusCode, header};
use axum::middleware::{self, Next};
use axum::response::{IntoResponse, Response};
use axum::routing::any;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use tower::util::MapResponseLayer;
use tower::{Layer as _, ServiceExt as _};
use tower_http::CompressionLevel;
use tower_http::compression::CompressionLayer;
use tower_http::compression::predicate::{DefaultPredicate, NotForContentType, Predicate};
use tower_http::services::{ServeDir, ServeFile};
use tower_http::set_header::SetResponseHeaderLayer;
use tower_http::set_header::response::SetMultipleResponseHeadersLayer;

const ASSETS_DIR: &str = "assets";

const INDEX_FILE: &str = "index.html";

const IMMUTABLE: &str = "public, max-age=31536000, immutable";

const REVALIDATE: &str = "no-cache";

const LOOPBACK_HOSTS: &[&str] = &["localhost", "127.0.0.1", "[::1]"];

const LOCALHOST_SUFFIX: &str = ".localhost";

const SECURITY_HEADERS: &[(&str, &str)] = &[
  ("x-content-type-options", "nosniff"),
  ("referrer-policy", "strict-origin-when-cross-origin"),
  ("cross-origin-opener-policy", "same-origin"),
  ("cross-origin-resource-policy", "same-origin"),
  ("content-security-policy", "frame-ancestors 'none'"),
];

#[derive(Clone, Debug, Default)]
pub struct WebOptions {
  pub static_dir: Option<PathBuf>,
  pub allowed_hosts: Vec<String>,
}

pub(crate) fn web_router(api: Router, options: &WebOptions) -> Router {
  let fallback = Router::new()
    .route("/api", any(not_found))
    .route("/api/{*path}", any(not_found))
    .route_layer(SetResponseHeaderLayer::if_not_present(
      header::CACHE_CONTROL,
      HeaderValue::from_static(NO_STORE),
    ));
  let fallback = match &options.static_dir {
    Some(dir) => fallback
      .nest_service(
        &format!("/{ASSETS_DIR}"),
        MapResponseLayer::new(cache_control(IMMUTABLE)).layer(web_files(&dir.join(ASSETS_DIR))),
      )
      .fallback_service(
        MapResponseLayer::new(cache_control(REVALIDATE))
          .layer(web_files(dir).fallback(page_fallback.with_state(index_file(dir)))),
      ),
    None => fallback,
  };
  let router = api.fallback_service(fallback);
  let hosts = Arc::new(AllowedHosts::new(&options.allowed_hosts));
  router
    .layer(
      CompressionLayer::new()
        .quality(CompressionLevel::Fastest)
        .compress_when(compress_predicate()),
    )
    .layer(middleware::from_fn_with_state(hosts, check_host))
    .layer(SetMultipleResponseHeadersLayer::if_not_present(
      SECURITY_HEADERS
        .iter()
        .map(|(name, value)| (HeaderName::from_static(name), HeaderValue::from_static(value)).into())
        .collect(),
    ))
}

fn web_files(dir: &Path) -> ServeDir {
  ServeDir::new(dir).precompressed_br().precompressed_gzip()
}

fn index_file(dir: &Path) -> ServeFile {
  ServeFile::new(dir.join(INDEX_FILE))
    .precompressed_br()
    .precompressed_gzip()
}

async fn page_fallback(State(index): State<ServeFile>, request: Request) -> Response {
  if is_page_navigation(&request) {
    index
      .oneshot(request)
      .await
      .map(|response| response.map(Body::new))
      .into_response()
  } else {
    StatusCode::NOT_FOUND.into_response()
  }
}

fn cache_control<B>(policy: &'static str) -> impl Fn(http::Response<B>) -> http::Response<B> + Clone {
  move |mut response| {
    if response.status().is_success() || response.status() == StatusCode::NOT_MODIFIED {
      response
        .headers_mut()
        .insert(header::CACHE_CONTROL, HeaderValue::from_static(policy));
    }
    response
  }
}

fn is_page_navigation(request: &Request) -> bool {
  matches!(*request.method(), Method::GET | Method::HEAD)
    && request
      .headers()
      .get_all(header::ACCEPT)
      .iter()
      .filter_map(|value| value.to_str().ok())
      .any(|accept| accept.contains("text/html"))
}

fn compress_predicate() -> impl Predicate {
  DefaultPredicate::new()
    .and(NotForContentType::const_new("application/x-xz"))
    .and(NotForContentType::const_new("application/gzip"))
    .and(NotForContentType::const_new("application/zstd"))
    .and(NotForContentType::const_new("font/"))
}

async fn check_host(State(hosts): State<Arc<AllowedHosts>>, request: Request, next: Next) -> Response {
  let host = request
    .headers()
    .get(header::HOST)
    .and_then(|value| value.to_str().ok())
    .and_then(|value| value.parse::<Authority>().ok())
    .map(|authority| authority.host().to_ascii_lowercase());
  match host {
    Some(host) if !hosts.allows(&host) => plain_error(
      ErrorCode::Forbidden,
      format!("the server does not answer requests for host `{host}`"),
    ),
    _ => next.run(request).await,
  }
}

struct AllowedHosts {
  names: Vec<String>,
}

impl AllowedHosts {
  fn new(configured: &[String]) -> Self {
    Self {
      names: configured.iter().map(|name| name.to_ascii_lowercase()).collect(),
    }
  }

  fn allows(&self, host: &str) -> bool {
    LOOPBACK_HOSTS.contains(&host) || host.ends_with(LOCALHOST_SUFFIX) || self.names.iter().any(|name| name == host)
  }
}

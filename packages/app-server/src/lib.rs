#[cfg(test)]
mod __tests__;

mod confine;
mod error;
mod events;
mod openapi;
pub mod routes;
pub mod state;

use crate::state::{AppState, ServerConfig};
use axum::Router;
use eyre::Report;
use std::path::PathBuf;
use std::sync::Arc;
use tower_http::cors::CorsLayer;
use tower_http::services::{ServeDir, ServeFile};

pub fn create_router(config: ServerConfig, static_dir: Option<PathBuf>) -> Result<Router, Report> {
  let api = routes::api_routes(Arc::new(AppState::new(config)?));

  let router = match static_dir {
    Some(static_dir) => {
      let index = static_dir.join("index.html");
      api.fallback_service(ServeDir::new(&static_dir).fallback(ServeFile::new(index)))
    },
    None => api,
  };
  Ok(router.layer(CorsLayer::permissive()))
}

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::init::global::global_init;

  #[ctor(unsafe)]
  fn init() {
    global_init();
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}

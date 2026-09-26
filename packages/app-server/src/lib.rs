#[cfg(test)]
mod __tests__;

mod api;
mod confine;
mod error;
mod events;
mod openapi;
pub mod routes;
pub mod state;

use crate::routes::api_router;
use crate::state::{ServerConfig, server_service};
use axum::Router;
use eyre::Report;
use std::path::PathBuf;
use tower_http::cors::CorsLayer;
use tower_http::services::{ServeDir, ServeFile};

pub fn create_router(config: ServerConfig, static_dir: Option<PathBuf>) -> Result<Router, Report> {
  let (api, _) = api_router(server_service(&config)?, config)?;

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

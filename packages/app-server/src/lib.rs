#[cfg(test)]
mod __tests__;

mod api;
mod app_settings_routes;
mod confine;
mod error;
mod events;
mod openapi;
pub mod routes;
pub mod state;
pub mod web;

use crate::routes::api_router;
use crate::state::{ServerConfig, server_service};
use crate::web::{WebOptions, web_router};
use axum::Router;
use eyre::Report;

pub fn create_router(config: ServerConfig, options: &WebOptions) -> Result<Router, Report> {
  let (api, _) = api_router(server_service(&config)?, config)?;
  Ok(web_router(api, options))
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

#[cfg(test)]
mod __tests__;

mod api;
mod app_settings_routes;
mod confine;
pub mod error;
mod events;
mod openapi;
pub mod port;
pub mod routes;
pub mod state;
#[cfg(not(target_family = "wasm"))]
pub mod web;

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

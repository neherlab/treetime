#[cfg(test)]
mod __tests__;

pub mod command;
pub mod command_config;
pub mod commands;
pub mod config;
pub mod examples_download;
#[cfg(not(target_family = "wasm"))]
mod examples_fetch;
pub mod json_value;
#[cfg(not(target_family = "wasm"))]
pub mod rtt_chart;
#[cfg(not(target_family = "wasm"))]
mod rtt_chart_render;
pub mod version;
pub mod yaml;

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

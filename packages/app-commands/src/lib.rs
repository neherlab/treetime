#[cfg(test)]
mod __tests__;

pub mod atomic_write;
pub mod bridge;
#[cfg(feature = "clap")]
pub mod check_config;
pub mod check_inputs;
pub mod command;
pub mod commands;
pub mod config;
pub mod datasets;
pub mod job;
pub mod json_float;
pub mod results;
pub mod rtt_chart;
mod rtt_chart_render;
pub mod run_checks;
#[cfg(feature = "clap")]
pub mod run_config;
pub mod runs;
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

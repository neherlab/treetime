#[cfg(test)]
mod __tests__;

pub mod check_inputs;
pub mod command;
pub mod commands;
pub mod config;
pub mod job;
pub mod json_float;
pub mod rtt_chart;
mod rtt_chart_render;
pub mod runs;

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

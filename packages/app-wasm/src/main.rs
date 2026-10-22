#[cfg(test)]
mod __tests__;

pub mod backend;
pub mod chunks;
pub mod exchange;
pub mod job;
pub mod jobs;
pub mod start;

fn main() {}

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

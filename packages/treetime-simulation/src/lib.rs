pub mod branch_lengths;
pub mod coalescent;
pub mod evolve;
pub mod gaps;
pub mod substitution_model;
pub mod tc_schedule;

pub const FALSE_FAIL_RATE: f64 = 1e-3;

#[cfg(test)]
mod __tests__;

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

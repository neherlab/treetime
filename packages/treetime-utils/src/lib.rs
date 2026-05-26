pub mod array;
pub mod collections;
pub mod datetime;
pub mod error;
pub mod fmt;
pub mod init;
pub mod interval;
pub mod io;
pub mod iterator;
pub mod sync;
pub mod testing;

#[cfg(test)]
mod tests {
  use crate::init::global::global_init;
  use ctor::ctor;

  #[ctor]
  fn init() {
    global_init();
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}

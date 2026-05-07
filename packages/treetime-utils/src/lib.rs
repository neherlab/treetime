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

#[cfg(any(
  all(target_arch = "aarch64", target_os = "linux", target_env = "gnu"),
  all(target_arch = "aarch64", target_os = "linux", target_env = "musl"),
  all(target_arch = "x86_64", target_os = "linux", target_env = "gnu"),
  all(target_arch = "x86_64", target_os = "linux", target_env = "musl"),
))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

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

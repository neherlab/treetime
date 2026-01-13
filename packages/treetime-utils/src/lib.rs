pub mod assert;
pub mod clap_styles;
pub mod compression;
pub mod console;
pub mod container;
pub mod datetime;
pub mod error;
pub mod file;
pub mod float_fmt;
pub mod fs;
pub mod global_init;
pub mod interval;
pub mod iter;
pub mod iterator;
pub mod manyzip;
pub mod mutex;
pub mod ndarray;
pub mod openblas;
pub mod random;
pub mod serde;
pub mod string;
pub mod vec;

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
  use crate::global_init::global_init;
  use ctor::ctor;

  #[ctor]
  fn init() {
    global_init();
  }
}

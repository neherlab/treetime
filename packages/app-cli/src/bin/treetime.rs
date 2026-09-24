#[cfg(any(
  all(target_arch = "x86_64", target_os = "linux", target_env = "gnu"),
  all(target_arch = "x86_64", target_os = "linux", target_env = "musl"),
  all(target_arch = "aarch64", target_os = "linux", target_env = "gnu"),
  all(target_arch = "aarch64", target_os = "linux", target_env = "musl"),
))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use app_cli::run::run_cli;
use ctor::ctor;
use eyre::Report;
use treetime_utils::init::global::global_init;

#[ctor(unsafe)]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  run_cli()
}

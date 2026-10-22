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
use std::process::ExitCode;
use treetime_utils::init::cpu::require_cpu_features;
use treetime_utils::init::global::{install_panic_hook, stderr_has_colors};

#[ctor(unsafe)]
fn init() {
  require_cpu_features();
  install_panic_hook(stderr_has_colors());
}

fn main() -> ExitCode {
  run_cli()
}

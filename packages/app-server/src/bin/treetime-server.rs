#[cfg(any(
  all(target_arch = "x86_64", target_os = "linux", target_env = "gnu"),
  all(target_arch = "x86_64", target_os = "linux", target_env = "musl"),
  all(target_arch = "aarch64", target_os = "linux", target_env = "gnu"),
  all(target_arch = "aarch64", target_os = "linux", target_env = "musl"),
))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use app_server::create_router;
use app_server::state::ServerConfig;
use clap::Parser;
use ctor::ctor;
use log::{LevelFilter, warn};
use std::io::{self, Write};
use std::path::PathBuf;
use std::thread::available_parallelism;
use treetime_utils::env::env_var_optional;
use treetime_utils::init::global::{global_init, setup_logger};

const HOST_ENV: &str = "HOST";
const PORT_ENV: &str = "PORT";
const STATIC_DIR_ENV: &str = "STATIC_DIR";
const DEFAULT_HOST: &str = "127.0.0.1";
const DEFAULT_PORT: u16 = 3100;

#[ctor]
fn init() {
  global_init();
}

#[derive(Parser, Debug)]
#[command(name = "treetime-server", about = "TreeTime web server")]
struct ServerArgs {
  /// Host address to bind to. Falls back to HOST env var, then 127.0.0.1
  #[arg(long)]
  host: Option<String>,

  /// Port to listen on. Falls back to PORT env var, then 3100
  #[arg(long, short)]
  port: Option<u16>,

  /// Number of processing threads. Defaults to all available CPU threads
  #[arg(long, short = 'j', default_value_t = default_jobs())]
  jobs: usize,

  /// Directory containing input datasets
  #[arg(long)]
  data_dir: PathBuf,

  /// Base directory for output files
  #[arg(long)]
  out_dir: PathBuf,
}

fn default_jobs() -> usize {
  available_parallelism().map_or(1, |n| n.get())
}

#[tokio::main]
async fn main() -> eyre::Result<()> {
  setup_logger(LevelFilter::Warn);
  let args = ServerArgs::parse();

  if args.jobs == 1 {
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .use_current_thread()
      .build_global()?;
  } else {
    rayon::ThreadPoolBuilder::new().num_threads(args.jobs).build_global()?;
  }

  let host = match args.host {
    Some(host) => host,
    None => env_var_optional(HOST_ENV)?.unwrap_or_else(|| DEFAULT_HOST.to_owned()),
  };

  let port = match args.port {
    Some(port) => port,
    None => match env_var_optional(PORT_ENV)? {
      Some(value) => value.parse().unwrap_or_else(|error| {
        warn!("Invalid {PORT_ENV} value '{value}' ({error}), using default {DEFAULT_PORT}");
        DEFAULT_PORT
      }),
      None => DEFAULT_PORT,
    },
  };

  let config = ServerConfig {
    data_dir: args.data_dir,
    out_dir: args.out_dir,
  };

  let static_dir = env_var_optional(STATIC_DIR_ENV)?;

  let addr = format!("{host}:{port}");
  let listener = tokio::net::TcpListener::bind(&addr).await?;
  writeln!(io::stderr().lock(), "TreeTime server listening on http://{addr}")?;
  axum::serve(listener, create_router(config, static_dir))
    .with_graceful_shutdown(shutdown_signal())
    .await?;
  Ok(())
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
#[cfg_attr(
  dylint_lib = "treetime_lints",
  expect(
    debug_remnants,
    reason = "status line on stderr; the shutdown future has no error channel"
  )
)]
async fn shutdown_signal() {
  let ctrl_c = tokio::signal::ctrl_c();
  let mut sigterm = tokio::signal::unix::signal(tokio::signal::unix::SignalKind::terminate())
    .expect("SIGTERM handler registration failed");
  tokio::select! {
    _ = ctrl_c => {},
    _ = sigterm.recv() => {},
  }
  eprintln!("TreeTime server shutting down");
}

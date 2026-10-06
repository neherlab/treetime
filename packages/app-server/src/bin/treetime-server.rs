#[cfg(any(
  all(target_arch = "x86_64", target_os = "linux", target_env = "gnu"),
  all(target_arch = "x86_64", target_os = "linux", target_env = "musl"),
  all(target_arch = "aarch64", target_os = "linux", target_env = "gnu"),
  all(target_arch = "aarch64", target_os = "linux", target_env = "musl"),
))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use app_server::create_router;
use app_server::state::{DEFAULT_MAX_UPLOAD_SIZE, ServerConfig};
use app_server::web::WebOptions;
use clap::Parser;
use ctor::ctor;
use eyre::WrapErr;
use log::{LevelFilter, warn};
use std::io::{self, Write};
use std::path::PathBuf;
use tokio_util::sync::CancellationToken;
use treetime_grid::MaxGridPoints;
use treetime_utils::env::env_var_optional;
use treetime_utils::init::global::{global_init, setup_logger};
use treetime_utils::init::thread_pool::{available_jobs, init_thread_pool};

const HOST_ENV: &str = "HOST";
const PORT_ENV: &str = "PORT";
const STATIC_DIR_ENV: &str = "STATIC_DIR";
const ALLOWED_HOSTS_ENV: &str = "ALLOWED_HOSTS";
const DEFAULT_HOST: &str = "127.0.0.1";
const DEFAULT_PORT: u16 = 3100;

#[ctor(unsafe)]
fn init() {
  global_init();
}

#[tokio::main]
async fn main() -> eyre::Result<()> {
  setup_logger(LevelFilter::Warn);
  let args = ServerArgs::parse();

  init_thread_pool(args.jobs)?;

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

  let allowed_hosts = if args.allowed_host.is_empty() {
    env_var_optional(ALLOWED_HOSTS_ENV)?.map_or_else(Vec::new, |value| {
      value
        .split(',')
        .map(str::trim)
        .filter(|host| !host.is_empty())
        .map(str::to_owned)
        .collect()
    })
  } else {
    args.allowed_host
  };

  let shutdown = CancellationToken::new();
  let config = ServerConfig {
    examples_dir: args.examples_dir,
    runs_dir: args.runs_dir,
    max_upload_size: args.max_upload_size,
    max_grid_points: Some(args.max_grid_points),
    shutdown: shutdown.clone(),
    settings: None,
  };

  let options = WebOptions {
    static_dir: env_var_optional(STATIC_DIR_ENV)?.map(PathBuf::from),
    allowed_hosts,
  };

  let addr = format!("{host}:{port}");
  let listener = tokio::net::TcpListener::bind(&addr)
    .await
    .wrap_err_with(|| format!("When binding the server to {addr}"))?;
  writeln!(io::stderr().lock(), "TreeTime server listening on http://{addr}")
    .wrap_err("When writing the startup line to standard error")?;
  axum::serve(listener, create_router(config, &options)?)
    .with_graceful_shutdown(async move {
      shutdown_signal().await;
      shutdown.cancel();
    })
    .await
    .wrap_err("When serving HTTP requests")?;
  Ok(())
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
  #[arg(long, short = 'j', default_value_t = available_jobs())]
  jobs: usize,

  /// Folder of the example datasets and example configurations that the app lists. Relative input paths of a run
  /// resolve inside it
  #[arg(long)]
  examples_dir: PathBuf,

  /// Directory that holds the runs: one folder per run with its record, events, inputs and outputs
  #[arg(long)]
  runs_dir: PathBuf,

  /// Largest total size, in bytes, of the files uploaded into one run. Defaults to 1 GiB
  #[arg(long, default_value_t = DEFAULT_MAX_UPLOAD_SIZE)]
  max_upload_size: usize,

  /// Largest number of points of one probability grid that a run on this server may use during time inference. A
  /// run that sets no value gets this one; a run may set a lower value, and a higher value is rejected
  #[arg(long, default_value_t = MaxGridPoints::default())]
  max_grid_points: MaxGridPoints,

  /// Host name the server answers besides localhost, `*.localhost`, 127.0.0.1 and [::1]; repeat for several.
  /// Falls back to the comma-separated ALLOWED_HOSTS env var. Requests that name another host get HTTP 403
  #[arg(long)]
  allowed_host: Vec<String>,
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

#[cfg(test)]
mod tests {
  use super::ServerArgs;
  use clap::Parser;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::default(&[],                                1_000_000)]
  #[case::minimum(&["--max-grid-points", "1000"],     1_000)]
  #[case::custom( &["--max-grid-points", "250000"],   250_000)]
  #[trace]
  fn test_treetime_server_max_grid_points(#[case] flag: &[&str], #[case] expected: usize) {
    let args = ServerArgs::try_parse_from(server_args(flag)).unwrap();
    assert_eq!(expected, args.max_grid_points.get());
  }

  #[test]
  fn test_treetime_server_max_grid_points_below_minimum_fails() {
    let error = ServerArgs::try_parse_from(server_args(&["--max-grid-points", "999"])).unwrap_err();
    assert_eq!(
      "error: invalid value '999' for '--max-grid-points <MAX_GRID_POINTS>': the grid point limit must be at least 1000 points, got 999\n\nFor more information, try '--help'.\n",
      error.to_string()
    );
  }

  fn server_args<'a>(extra: &[&'a str]) -> Vec<&'a str> {
    [
      &["treetime-server", "--examples-dir", "examples", "--runs-dir", "runs"],
      extra,
    ]
    .concat()
  }
}

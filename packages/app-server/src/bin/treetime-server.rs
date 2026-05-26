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
use log::LevelFilter;
use std::path::PathBuf;
use std::thread::available_parallelism;
use treetime_utils::init::global::{global_init, setup_logger};

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

  let host = args
    .host
    .or_else(|| std::env::var("HOST").ok())
    .unwrap_or_else(|| "127.0.0.1".to_owned());

  let port: u16 = args
    .port
    .or_else(|| {
      std::env::var("PORT").ok().and_then(|val| {
        val.parse().ok().or_else(|| {
          eprintln!("Warning: invalid PORT value '{val}', using default 3100");
          None
        })
      })
    })
    .unwrap_or(3100);

  let config = ServerConfig {
    data_dir: args.data_dir,
    out_dir: args.out_dir,
  };

  let addr = format!("{host}:{port}");
  let listener = tokio::net::TcpListener::bind(&addr).await?;
  eprintln!("TreeTime server listening on http://{addr}");
  axum::serve(listener, create_router(config))
    .with_graceful_shutdown(shutdown_signal())
    .await?;
  Ok(())
}

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

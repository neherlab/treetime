use app_server::create_router;
use app_server::state::ServerConfig;
use clap::Parser;
use ctor::ctor;
<<<<<<< HEAD
<<<<<<< HEAD
use log::info;
use treetime_utils::init::global::global_init;
=======
use log::LevelFilter;
use std::path::PathBuf;
use std::thread::available_parallelism;
use treetime_utils::init::global::{global_init, setup_logger};
>>>>>>> 6fc31936 (feat(server): add CLI args for jobs, data dir, and output dir)
=======
use log::LevelFilter;
use treetime_utils::init::global::{global_init, setup_logger};
>>>>>>> 3b85da3a (feat(app): add dataset discovery endpoint)

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
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
  let host = std::env::var("HOST").unwrap_or("127.0.0.1".to_owned());
=======
  setup_logger(LevelFilter::Warn);
=======
  setup_logger(LevelFilter::Info);
>>>>>>> 3b85da3a (feat(app): add dataset discovery endpoint)
  let host = std::env::var("HOST").unwrap_or_else(|_| "127.0.0.1".to_owned());
>>>>>>> f5bb3671 (fix(server): default log level to warn (matching CLI default))
=======
  let host = std::env::var("HOST").unwrap_or_else(|_| "127.0.0.1".to_owned());
>>>>>>> b18af097 (fix(app-server): use structured error format and fix unwrap_or_default)
  let port: u16 = match std::env::var("PORT") {
    Ok(val) => val.parse().unwrap_or_else(|_| {
      eprintln!("Warning: invalid PORT value '{val}', using default 3100");
      3100
    }),
    Err(_) => 3100,
=======
  setup_logger(LevelFilter::Warn);
  let args = ServerArgs::parse();

  if args.jobs == 1 {
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .use_current_thread()
      .build_global()?;
  } else {
    rayon::ThreadPoolBuilder::new()
      .num_threads(args.jobs)
      .build_global()?;
  }

  let host = args
    .host
    .or_else(|| std::env::var("HOST").ok())
    .unwrap_or_else(|| "127.0.0.1".to_owned());

  let port: u16 = args.port.or_else(|| {
    std::env::var("PORT").ok().and_then(|val| {
      val.parse().ok().or_else(|| {
        eprintln!("Warning: invalid PORT value '{val}', using default 3100");
        None
      })
    })
  }).unwrap_or(3100);

  let config = ServerConfig {
    data_dir: args.data_dir,
    out_dir: args.out_dir,
>>>>>>> 6fc31936 (feat(server): add CLI args for jobs, data dir, and output dir)
  };

  let addr = format!("{host}:{port}");
  let listener = tokio::net::TcpListener::bind(&addr).await?;
<<<<<<< HEAD
  info!("TreeTime server listening on http://{addr}");
  axum::serve(listener, create_router()).await?;
=======
  eprintln!("TreeTime server listening on http://{addr}");
  axum::serve(listener, create_router(config))
    .with_graceful_shutdown(shutdown_signal())
    .await?;
>>>>>>> a9124f4b (fix: add graceful shutdown to treetime-server on SIGINT)
  Ok(())
}

async fn shutdown_signal() {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
  tokio::signal::ctrl_c().await.ok();
=======
  drop(tokio::signal::ctrl_c().await);
>>>>>>> 5c566b97 (refactor: lint)
=======
  let _ = tokio::signal::ctrl_c().await;
>>>>>>> 7c23b04d (fix(api,napi,server): clean up Cargo deps, fix error contract, fix code style)
=======
  let ctrl_c = tokio::signal::ctrl_c();
  let mut sigterm = tokio::signal::unix::signal(tokio::signal::unix::SignalKind::terminate())
    .expect("SIGTERM handler registration failed");
  tokio::select! {
    _ = ctrl_c => {},
    _ = sigterm.recv() => {},
  }
>>>>>>> 2d60ad7b (feat(server): add SIGTERM handling for graceful shutdown on restart)
  eprintln!("TreeTime server shutting down");
}

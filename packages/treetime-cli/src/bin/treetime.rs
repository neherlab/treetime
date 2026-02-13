use ctor::ctor;
use eyre::Report;
use log::info;
use treetime::commands::ancestral::run::run_ancestral_reconstruction;
use treetime::commands::clock::run::run_clock;
use treetime::commands::homoplasy::run::run_homoplasy;
use treetime::commands::mugration::run::run_mugration;
use treetime::commands::optimize::run::run_optimize;
use treetime::commands::prune::run::run_prune;
use treetime::commands::timetree::run::run_timetree_estimation;
use treetime_cli::cli::treetime_cli::{TreetimeCommands, generate_shell_completions, treetime_parse_cli_args};
use treetime_io::json::{JsonPretty, json_write_str};
use treetime_utils::init::global::global_init;
use treetime_utils::init::openblas::get_openblas_info_str;

#[ctor]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  let args = treetime_parse_cli_args()?;

  info!("# Command line arguments");
  info!("{}", json_write_str(&args, JsonPretty(true))?);

  if args.jobs.jobs == 1 {
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .use_current_thread()
      .build_global()?;
  } else {
    rayon::ThreadPoolBuilder::new()
      .num_threads(args.jobs.jobs)
      .build_global()?;
  }

  match args.command {
    TreetimeCommands::Timetree(timetree_args) => {
      run_timetree_estimation(&timetree_args)?;
    },
    TreetimeCommands::Optimize(optimize_args) => {
      run_optimize(&optimize_args)?;
    },
    TreetimeCommands::Prune(prune_args) => {
      run_prune(&prune_args)?;
    },
    TreetimeCommands::Ancestral(ancestral_args) => {
      run_ancestral_reconstruction(&ancestral_args)?;
    },
    TreetimeCommands::Clock(clock_args) => {
      run_clock(&clock_args)?;
    },
    TreetimeCommands::Homoplasy(homoplasy_args) => {
      run_homoplasy(homoplasy_args)?;
    },
    TreetimeCommands::Mugration(mugration_args) => {
      run_mugration(&mugration_args)?;
    },
    TreetimeCommands::Completions { shell } => {
      generate_shell_completions(&shell)?;
    },
    TreetimeCommands::Arg(arg_args) => {},
    TreetimeCommands::Debug => {
      println!("{}", get_openblas_info_str());
    },
  }

  Ok(())
}

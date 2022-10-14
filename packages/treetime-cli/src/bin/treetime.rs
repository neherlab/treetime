#![allow(clippy::self_named_module_files)]

use ctor::ctor;
use eyre::Report;
use log::info;
use treetime::commands::ancestral::run_ancestral_reconstruction::run_ancestral_reconstruction;
use treetime::commands::clock::run_clock::run_clock;
use treetime::commands::homoplasy::run_homoplasy::run_homoplasy;
use treetime::commands::mugration::run_mugration::run_mugration;
use treetime::commands::timetree::run_timetree_estimation::run_timetree_estimation;
use treetime::utils::global_init::global_init;
use treetime_cli::cli::treetime_cli::{generate_shell_completions, treetime_parse_cli_args, TreetimeCommands};

#[cfg(all(target_os = "linux", target_arch = "x86_64"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[ctor]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  let args = treetime_parse_cli_args()?;

  info!("{:#?}", &args);

  rayon::ThreadPoolBuilder::new().num_threads(args.jobs).build_global()?;

  match args.command {
    TreetimeCommands::Timetree(timetree_args) => {
      run_timetree_estimation(&timetree_args)?;
    }
    TreetimeCommands::Ancestral(ancestral_args) => {
      run_ancestral_reconstruction(&ancestral_args)?;
    }
    TreetimeCommands::Clock(clock_args) => {
      run_clock(&clock_args)?;
    }
    TreetimeCommands::Homoplasy(homoplasy_args) => {
      run_homoplasy(homoplasy_args)?;
    }
    TreetimeCommands::Mugration(mugration_args) => {
      run_mugration(&mugration_args)?;
    }
    TreetimeCommands::Completions { shell } => {
      generate_shell_completions(&shell)?;
    }
    TreetimeCommands::Arg(arg_args) => {}
  }

  Ok(())
}

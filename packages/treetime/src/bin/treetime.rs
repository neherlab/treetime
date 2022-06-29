// use ctor::ctor;
// use eyre::Report;
// use log::info;
// use treetime::ancestral::run_ancestral::run_ancestral;
// use treetime::cli::treetime_cli::{generate_shell_completions, treetime_parse_cli_args, TreetimeCommands};
// use treetime::clock::run_clock::run_clock;
// use treetime::homoplasy::run_homoplasy::run_homoplasy;
// use treetime::mugration::run_mugration::run_mugration;
// use treetime::utils::global_init::global_init;
//
// #[cfg(all(target_family = "linux", target_arch = "x86_64"))]
// #[global_allocator]
// static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
//
// #[ctor]
// fn init() {
//   global_init();
// }
//
// fn main() -> Result<(), Report> {
//   let args = treetime_parse_cli_args()?;
//
//   info!("{:#?}", &args);
//
//   match args.command {
//     TreetimeCommands::Ancestral(ancestral_args) => {
//       run_ancestral(ancestral_args)?;
//     }
//     TreetimeCommands::Clock(clock_args) => {
//       run_clock(clock_args)?;
//     }
//     TreetimeCommands::Homoplasy(homoplasy_args) => {
//       run_homoplasy(homoplasy_args)?;
//     }
//     TreetimeCommands::Mugration(mugration_args) => {
//       run_mugration(mugration_args)?;
//     }
//     TreetimeCommands::Completions { shell } => {
//       generate_shell_completions(&shell)?;
//     }
//   }
//
//   Ok(())
// }

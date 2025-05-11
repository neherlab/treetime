use clap::Parser;
use color_eyre::{Section, SectionExt};
use ctor::ctor;
use eyre::Report;
use treetime::io::compression::remove_compression_ext;
use treetime::make_report;
use treetime::utils::global_init::{global_init, setup_logger};
use treetime_cli::convert::args::{Args, guess_tree_format_from_filename};
use treetime_cli::convert::convert::{ConverterGraph, converter_read_file, converter_write_file};

#[ctor]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

  let args = Args::parse();
  setup_logger(args.verbosity.get_filter_level());

  let input_format = args
    .input_format
    .or_else(|| guess_tree_format_from_filename(remove_compression_ext(&args.input)))
    .ok_or_else(|| {
      make_report!("Input format was not specified and unable to autodetect. Please provide --input-format argument")
    })
    .with_section(|| format!("{:#?}", &args.input).header("Input file:"))?;

  let output_format = args
    .output_format
    .or_else(|| guess_tree_format_from_filename(remove_compression_ext(&args.output)))
    .ok_or_else(|| {
      make_report!("Output format was not specified and unable to autodetect. Please provide --output-format argument")
    })
    .with_section(|| format!("{:#?}", &args.input).header("Output file:"))?;

  let graph: ConverterGraph = converter_read_file(&args, input_format)?;

  converter_write_file(&args, output_format, &graph)?;

  Ok(())
}

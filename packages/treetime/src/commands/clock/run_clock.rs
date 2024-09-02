use crate::cli::rtt_chart::{
  print_clock_regression_chart, write_clock_regression_chart_png, write_clock_regression_chart_svg,
};
use crate::commands::clock::assign_dates::assign_dates;
use crate::commands::clock::clock_args::TreetimeClockArgs;
use crate::commands::clock::clock_graph::ClockGraph;
use crate::commands::clock::clock_regression::{run_clock_regression, ClockOptions};
use crate::commands::clock::reroot::reroot_in_place;
use crate::commands::clock::rtt::{gather_clock_regression_results, write_clock_regression_result_csv};
use crate::io::dates_csv::read_dates;
use crate::io::graphviz::graphviz_write_file;
use crate::io::json::{json_write_file, JsonPretty};
use crate::io::nwk::{nwk_read_file, nwk_write_file, NwkWriteOptions};
use eyre::{Report, WrapErr};

pub fn run_clock(clock_args: &TreetimeClockArgs) -> Result<(), Report> {
  let TreetimeClockArgs {
    aln,
    tree,
    vcf_reference,
    dates,
    name_column,
    date_column,
    sequence_length,
    gtr,
    gtr_params,
    branch_length_mode,
    method_anc,
    clock_filter,
    reroot,
    keep_root,
    prune_short,
    tip_slack,
    covariation,
    allow_negative_rate,
    plot_rtt,
    outdir,
    seed,
  } = clock_args;

  let mut graph: ClockGraph = if let Some(tree) = tree {
    nwk_read_file(tree)
  } else {
    unimplemented!("Tree inference is not implemented")
  }?;

  {
    let dates = read_dates(dates, name_column, date_column).wrap_err("When reading dates")?;
    assign_dates(&graph, &dates)?;
  }

  let options = ClockOptions {
    variance_factor: 0.0,
    variance_offset: 0.0,
  };

  let clock_model = run_clock_regression(&graph, &options)?;
  if !keep_root {
    reroot_in_place(&mut graph, &options)?;
  }

  nwk_write_file(outdir.join("rerooted.nwk"), &graph, &NwkWriteOptions::default())?;
  json_write_file(outdir.join("graph_output.json"), &graph, JsonPretty(true))?;
  graphviz_write_file(outdir.join("graph_output.dot"), &graph)?;

  json_write_file(outdir.join("clock_model.json"), &clock_model, JsonPretty(true))?;

  let results = gather_clock_regression_results(&graph, &clock_model);

  write_clock_regression_result_csv(&results, outdir.join("clock.csv"), b',')?;
  write_clock_regression_chart_svg(&results, &clock_model, outdir.join("clock.svg"))?;
  write_clock_regression_chart_png(&results, &clock_model, outdir.join("clock.png"))?;

  print_clock_regression_chart(&results, &clock_model)?;

  Ok(())
}

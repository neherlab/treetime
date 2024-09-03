use crate::cli::rtt_chart::{
  print_clock_regression_chart, write_clock_regression_chart_png, write_clock_regression_chart_svg,
};
use crate::commands::clock::assign_dates::assign_dates;
use crate::commands::clock::clock_args::TreetimeClockArgs;
use crate::commands::clock::clock_graph::ClockGraph;
use crate::commands::clock::clock_regression::{run_clock_regression, ClockOptions};
use crate::commands::clock::reroot::reroot_in_place;
use crate::commands::clock::rtt::{gather_clock_regression_results, write_clock_regression_result_csv};
use crate::commands::clock::clock_filter::clock_filter_inplace;
use crate::io::dates_csv::read_dates;
use crate::io::graphviz::graphviz_write_file;
use crate::io::json::{json_write_file, JsonPretty};
use crate::io::nwk::{nwk_read_file, nwk_write_file, NwkWriteOptions};
use eyre::{Report, WrapErr};

use super::clock_model::ClockModel;

pub fn get_clock_model(graph: &mut ClockGraph, options: &ClockOptions, keep_root: bool) -> Result<ClockModel, Report> {
  if keep_root {
    run_clock_regression(&graph, &options)
  } else {
    reroot_in_place(graph, &options)?;
    let root = graph.get_exactly_one_root()?;
    let root = root.read_arc().payload().read_arc();
    ClockModel::new(&root.total)
  }
}

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

  let options = if *covariation {
    let seq_len = sequence_length.unwrap_or(0) as f64; // should error if sequence_length is None
    ClockOptions {
    variance_factor: 1.0/seq_len,
    variance_offset: 0.0,
    variance_offset_leaf: 9.0/seq_len/seq_len,
  }
  } else {
    ClockOptions {
      variance_factor: 0.0,
      variance_offset: 0.0,
      variance_offset_leaf: 1.0,
    }
  };

  // if we keep the root, run clock regression and report the model of the current root
  // otherwise, reroot the tree and use clock model of the novel root
  let mut clock_model = get_clock_model(&mut graph, &options, *keep_root)?;

  if *clock_filter>0.0 {
    let new_outliers = clock_filter_inplace(&graph, &clock_model, *clock_filter);
    dbg!(new_outliers);
    if new_outliers > 0 {
      clock_model = get_clock_model(&mut graph, &options, true)?;
    }
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

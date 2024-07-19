use crate::cli::rtt_chart::{draw_rtt_console_chart, write_rtt_svg_chart};
use crate::commands::clock::clock_args::TreetimeClockArgs;
use crate::commands::clock::clock_graph::{create_graph, infer_graph};
use crate::commands::clock::graph_regression::calculate_averages;
use crate::commands::clock::graph_regression_policy::GraphNodeRegressionPolicyReroot;
use crate::commands::clock::run_clock_model::{run_clock_model, RunClockModelParams, RunClockModelResults};
use crate::commands::clock::run_reroot::{run_reroot, RerootParams};
use crate::io::csv::CsvStructFileWriter;
use crate::io::dates_csv::read_dates;
use crate::io::graphviz::graphviz_write_file;
use crate::io::json::{json_write_file, JsonPretty};
use crate::io::nwk::{nwk_write_file, NwkWriteOptions};
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

  let dates = read_dates(dates, name_column, date_column).wrap_err("When reading dates")?;

  let mut graph = match tree {
    None => infer_graph()?,
    Some(tree) => create_graph(tree, &dates)?,
  };

  graphviz_write_file(outdir.join("graph_input.dot"), &graph)?;
  json_write_file(&graph, outdir.join("graph_input.json"), JsonPretty(true))?;

  if *clock_filter > 0.0 {
    unimplemented!("clock_filter")
  }

  calculate_averages::<GraphNodeRegressionPolicyReroot>(&mut graph);

  let slope = None;

  if !keep_root {
    if *covariation {
      unimplemented!("covariation")
    }

    run_reroot(
      &mut graph,
      &RerootParams {
        reroot: *reroot,
        slope,
        force_positive: false,
        keep_node_order: false,
      },
    )?;
  }

  let RunClockModelResults { clock_model, rtt } =
    run_clock_model::<GraphNodeRegressionPolicyReroot>(&mut graph, &RunClockModelParams { slope })?;

  let mut rtt_writer = CsvStructFileWriter::new(outdir.join("rtt.csv"), b',')?;
  rtt.iter().try_for_each(|result| rtt_writer.write(result))?;

  graphviz_write_file(outdir.join("graph_output.dot"), &graph)?;
  json_write_file(&graph, outdir.join("graph_output.json"), JsonPretty(true))?;

  nwk_write_file(outdir.join("rerooted.nwk"), &graph, &NwkWriteOptions::default())?;

  draw_rtt_console_chart(&rtt, &clock_model);
  write_rtt_svg_chart(outdir.join("rtt.svg"), &rtt, &clock_model)?;

  Ok(())
}

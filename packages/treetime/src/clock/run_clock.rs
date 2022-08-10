use crate::clock::clock_args::TreetimeClockArgs;
use crate::clock::clock_graph::{create_graph, infer_graph};
use crate::clock::graph_regression::calculate_averages;
use crate::clock::graph_regression_policy::GraphNodeRegressionPolicyReroot;
use crate::clock::run_clock_model::{run_clock_model, RunClockModelParams};
use crate::clock::run_reroot::{run_reroot, RerootParams};
use crate::io::dates::read_dates;
use crate::io::file::create_file;
use eyre::Report;

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

  let dates = read_dates(dates, name_column, date_column)?;

  let mut graph = match tree {
    None => infer_graph()?,
    Some(tree) => create_graph(tree, &dates)?,
  };

  graph.print_graph(create_file(outdir.join("graph_input.dot"))?)?;

  if *clock_filter > 0.0 {
    unimplemented!("clock_filter")
  }

  calculate_averages::<GraphNodeRegressionPolicyReroot>(&mut graph);

  let slope = 0.0;

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

  let result = run_clock_model::<GraphNodeRegressionPolicyReroot>(&mut graph, &RunClockModelParams { slope })?;

  graph.print_graph(create_file(outdir.join("graph_output.dot"))?)?;

  Ok(())
}

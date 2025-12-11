use clap::Parser;
use ctor::ctor;
use eyre::{Context, Report};
use indexmap::IndexMap;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use treetime::commands::timetree::coalescent::coalescent::compute_coalescent_contributions;
use treetime::commands::timetree::data::date_constraints::load_date_constraints;
use treetime::distribution::distribution::Distribution;
use treetime::graph::node::{GraphNodeKey, Named};
use treetime::io::dates_csv::read_dates;
use treetime::io::nwk::nwk_read_file;
use treetime::o;
use treetime::representation::graph_ancestral::GraphAncestral;
use treetime_convolution::grid::Grid;
use treetime_io::json::json_read_file;
use treetime_utils::float_fmt::FloatFormatExt;
use treetime_utils::global_init::global_init;
use treetime_utils::serde::{array1_from_vec, indexmap_array1_from_map};
use treetime_utils::string::truncate_right_with_ellipsis;

#[ctor]
fn init() {
  global_init();
}

#[derive(Parser, Clone, Serialize, Deserialize)]
#[command(
  name = "coalescent-validation",
  about = "Validates coalescent contribution computation against Python v0 baseline results",
  version
)]
pub struct Args {
  /// Tc values to run (comma-separated: 0.01, 0.1, 1.0, 10.0)
  #[arg(long, value_delimiter = ',', default_value = "0.01,0.1,1.0,10.0")]
  tc_values: Vec<String>,

  /// Output directory base path
  #[arg(long, default_value = "tmp")]
  output_dir: String,

  /// Absolute difference threshold for showing nodes in the detailed table.
  /// If no nodes exceed the threshold, the worst nodes are shown instead.
  #[arg(long, default_value_t = 0.001)]
  abs_diff_threshold: f64,

  /// Maximum number of nodes to show in the detailed table.
  #[arg(long, default_value_t = 25)]
  max_rows: usize,

  /// Optional statistics output file path
  #[arg(long)]
  stats_output: Option<String>,

  /// Enable detailed progress output
  #[arg(long, short = 'v')]
  verbose: bool,
}

#[derive(Debug, Deserialize)]
struct Snapshot {
  #[allow(dead_code)]
  description: String,
  inputs: SnapshotInputs,
  tbp_grid: SnapshotTbpGrid,
  #[serde(deserialize_with = "array1_from_vec")]
  #[allow(dead_code)]
  lineage_counts: Array1<f64>,
  #[serde(deserialize_with = "array1_from_vec")]
  #[allow(dead_code)]
  integral_merger_rate: Array1<f64>,
  #[serde(deserialize_with = "array1_from_vec")]
  #[allow(dead_code)]
  total_merger_rate: Array1<f64>,
  #[serde(deserialize_with = "indexmap_array1_from_map")]
  node_contributions: IndexMap<String, Array1<f64>>,
}

#[derive(Debug, Deserialize)]
struct SnapshotInputs {
  tree_path: String,
  metadata_path: String,
  tc: f64,
  #[allow(dead_code)]
  present_time: f64,
}

#[derive(Debug, Deserialize)]
struct SnapshotTbpGrid {
  start: f64,
  end: f64,
  #[allow(dead_code)]
  padding: f64,
  n_points: usize,
}

#[derive(Clone, Serialize, Deserialize)]
struct TestResult {
  name: String,
  t_grid_tbp: Array1<f64>,
  expected: IndexMap<String, Array1<f64>>,
  actual: IndexMap<String, Array1<f64>>,
}

#[derive(Clone, Serialize, Deserialize)]
struct NodeDiffRow {
  node: String,
  n_points: usize,
  max_abs_diff: f64,
  t_at_max_tbp: f64,
  mae: f64,
  rmse: f64,
  p95_abs_diff: f64,
}

#[derive(Clone, Serialize, Deserialize)]
struct DiffStats {
  n_expected_nodes: usize,
  n_actual_nodes: usize,
  n_compared_nodes: usize,
  n_missing_nodes: usize,
  n_extra_nodes: usize,

  points_per_node: usize,
  total_points_compared: usize,

  global_mae: f64,
  global_rmse: f64,

  node_mae_median: f64,
  node_rmse_median: f64,
  node_max_abs_p95: f64,
  node_max_abs_median: f64,

  worst_node: String,
  worst_t_tbp: f64,
  worst_max_abs_diff: f64,
}

fn main() -> Result<(), Report> {
  let args = Args::parse();
  let mut results = Vec::new();

  println!("# Coalescent Validation Results\n");

  if args.verbose {
    println!("Running tests for Tc values: {}\n", args.tc_values.join(", "));
  }

  for tc_str in &args.tc_values {
    let tc_value: f64 = tc_str
      .parse()
      .map_err(|e| eyre::eyre!("Invalid Tc value '{tc_str}': {e}"))?;
    let snapshot_filename = format!("coalescent_flu_h3n2_20_tc{tc_str}.json");

    match run_coalescent_test(&args, &snapshot_filename, tc_value) {
      Ok(result) => {
        print_comparison_table(&args, &result);
        results.push(result);
      },
      Err(e) => {
        eprintln!("Warning: Test for Tc={tc_value} failed: {e}");
      },
    }
  }

  if results.is_empty() {
    eprintln!("Error: No valid tests were run");
    return Ok(());
  }

  println!("\n## Aggregate Statistics Summary\n");
  print_aggregate_statistics(&results);

  if let Some(stats_output_path) = &args.stats_output {
    write_statistics_to_file(stats_output_path, &results)?;
    if args.verbose {
      println!("Statistics written to: {stats_output_path}");
    }
  }

  Ok(())
}

fn load_snapshot(filename: &str) -> Result<Snapshot, Report> {
  let path = Path::new("test_scripts/snapshots").join(filename);
  json_read_file::<Snapshot, _>(&path)
}

fn load_graph(snapshot: &Snapshot) -> Result<GraphAncestral, Report> {
  let graph = nwk_read_file(&snapshot.inputs.tree_path)?;
  let dates = read_dates(&snapshot.inputs.metadata_path, &Some(o!("name")), &Some(o!("date")))?;
  load_date_constraints(&dates, &graph)?;
  Ok(graph)
}

fn convert_map_keys_to_names<V>(graph: &GraphAncestral, map: &IndexMap<GraphNodeKey, V>) -> Result<IndexMap<String, V>, Report>
where
  V: Clone,
{
  let mut out = IndexMap::with_capacity(map.len());
  for (key, value) in map {
    let node = graph
      .get_node(*key)
      .ok_or_else(|| eyre::eyre!("Node key not found in graph: {key:?}"))?;
    let name = node
      .read_arc()
      .payload()
      .read_arc()
      .name()
      .ok_or_else(|| eyre::eyre!("Node has no name: {key:?}"))?
      .as_ref()
      .to_owned();
    out.insert(name, value.clone());
  }
  Ok(out)
}

fn run_coalescent_test(args: &Args, snapshot_filename: &str, tc_value: f64) -> Result<TestResult, Report> {
  println!("### Test: Coalescent contributions (Tc={tc_value})");

  let snapshot = load_snapshot(snapshot_filename)?;
  let graph = load_graph(&snapshot)?;
  let tc = Distribution::constant(snapshot.inputs.tc);

  let actuals = compute_coalescent_contributions(&graph, &tc)?;
  let actuals = convert_map_keys_to_names(&graph, &actuals).wrap_err("When mapping node keys to node names")?;

  let t_grid = Grid::from_range_n_points(
    snapshot.tbp_grid.start,
    snapshot.tbp_grid.end,
    snapshot.tbp_grid.n_points,
  )?;

  let t_grid_tbp = t_grid.to_array();

  let actuals: IndexMap<String, Array1<f64>> = actuals
    .iter()
    .map(|(name, dist)| {
      let y_resampled = dist
        .eval_many(&t_grid_tbp)
        .wrap_err_with(|| format!("When evaluating node contribution distribution for node '{name}'"))?;
      Ok((name.clone(), y_resampled))
    })
    .collect::<Result<IndexMap<String, Array1<f64>>, Report>>()?;

  if args.verbose {
    println!("  Nodes compared: {}", actuals.len());
    println!("  Grid points: {}", snapshot.tbp_grid.n_points);
  }

  Ok(TestResult {
    name: format!("Tc={tc_value}"),
    t_grid_tbp,
    expected: snapshot.node_contributions,
    actual: actuals,
  })
}

fn validate_comparable(
  t_grid_tbp: &Array1<f64>,
  expected: &IndexMap<String, Array1<f64>>,
  actual: &IndexMap<String, Array1<f64>>,
) -> Result<(Vec<String>, Vec<String>, Vec<String>), Report> {
  let expected_keys: Vec<String> = expected.keys().cloned().collect();
  let actual_keys: Vec<String> = actual.keys().cloned().collect();

  let missing: Vec<String> = expected_keys
    .iter()
    .filter(|k| !actual.contains_key(*k))
    .cloned()
    .collect();
  let extra: Vec<String> = actual_keys
    .iter()
    .filter(|k| !expected.contains_key(*k))
    .cloned()
    .collect();

  let common: Vec<String> = expected_keys
    .iter()
    .filter(|k| actual.contains_key(*k))
    .cloned()
    .collect();

  if t_grid_tbp.is_empty() {
    eyre::bail!("Time grid is empty");
  }

  for node in &common {
    let exp = &expected[node];
    let act = &actual[node];
    if exp.len() != act.len() {
      eyre::bail!(
        "Array length mismatch for node '{node}': expected {} points, actual {} points",
        exp.len(),
        act.len()
      );
    }
    if exp.len() != t_grid_tbp.len() {
      eyre::bail!(
        "Grid length mismatch for node '{node}': grid has {} points, expected has {} points",
        t_grid_tbp.len(),
        exp.len()
      );
    }
  }

  Ok((common, missing, extra))
}

fn median(sorted: &[f64]) -> f64 {
  if sorted.is_empty() {
    return 0.0;
  }
  let mid = sorted.len() / 2;
  if sorted.len() % 2 == 1 {
    sorted[mid]
  } else {
    0.5 * (sorted[mid - 1] + sorted[mid])
  }
}

fn percentile(sorted: &[f64], p: f64) -> f64 {
  if sorted.is_empty() {
    return 0.0;
  }
  let p = p.clamp(0.0, 1.0);
  let idx = ((p * ((sorted.len() - 1) as f64)).ceil() as usize).min(sorted.len() - 1);
  sorted[idx]
}

fn compute_diff_report(result: &TestResult) -> Result<(DiffStats, Vec<NodeDiffRow>), Report> {
  let (common, missing, extra) = validate_comparable(&result.t_grid_tbp, &result.expected, &result.actual)
    .wrap_err("When validating node sets and array lengths")?;

  if common.is_empty() {
    eyre::bail!("No comparable nodes (expected and actual node sets do not overlap)");
  }

  let mut rows = Vec::with_capacity(common.len());

  let mut total_abs = 0.0;
  let mut total_sq = 0.0;
  let mut total_points = 0_usize;

  let mut worst_node = String::new();
  let mut worst_t_tbp = 0.0;
  let mut worst_max_abs_diff = -1.0;

  for node in &common {
    let exp = &result.expected[node];
    let act = &result.actual[node];
    let n = exp.len();

    let mut abs_diffs = Vec::with_capacity(n);
    let mut sum_abs = 0.0;
    let mut sum_sq = 0.0;
    let mut max_abs = -1.0;
    let mut max_idx = 0_usize;

    for (i, (e, a)) in exp.iter().zip(act.iter()).enumerate() {
      let diff = e - a;
      if !diff.is_finite() {
        eyre::bail!("Non-finite difference for node '{node}' at index {i}: expected={e}, actual={a}");
      }
      let abs = diff.abs();
      abs_diffs.push(abs);
      sum_abs += abs;
      sum_sq += diff * diff;
      if abs > max_abs {
        max_abs = abs;
        max_idx = i;
      }
    }

    abs_diffs.sort_by_key(|&x| OrderedFloat(x));
    let mae = sum_abs / (n as f64);
    let rmse = (sum_sq / (n as f64)).sqrt();
    let p95_abs_diff = percentile(&abs_diffs, 0.95);
    let t_at_max_tbp = result.t_grid_tbp[max_idx];

    total_abs += sum_abs;
    total_sq += sum_sq;
    total_points += n;

    if max_abs > worst_max_abs_diff {
      worst_max_abs_diff = max_abs;
      worst_node = node.clone();
      worst_t_tbp = t_at_max_tbp;
    }

    rows.push(NodeDiffRow {
      node: node.clone(),
      n_points: n,
      max_abs_diff: max_abs,
      t_at_max_tbp,
      mae,
      rmse,
      p95_abs_diff,
    });
  }

  rows.sort_by_key(|r| OrderedFloat(-r.max_abs_diff));

  let mut node_maes: Vec<f64> = rows.iter().map(|r| r.mae).collect();
  node_maes.sort_by_key(|&x| OrderedFloat(x));
  let mut node_rmses: Vec<f64> = rows.iter().map(|r| r.rmse).collect();
  node_rmses.sort_by_key(|&x| OrderedFloat(x));
  let mut node_maxes: Vec<f64> = rows.iter().map(|r| r.max_abs_diff).collect();
  node_maxes.sort_by_key(|&x| OrderedFloat(x));

  let points_per_node = result.t_grid_tbp.len();
  let global_mae = total_abs / (total_points as f64);
  let global_rmse = (total_sq / (total_points as f64)).sqrt();

  let stats = DiffStats {
    n_expected_nodes: result.expected.len(),
    n_actual_nodes: result.actual.len(),
    n_compared_nodes: rows.len(),
    n_missing_nodes: missing.len(),
    n_extra_nodes: extra.len(),
    points_per_node,
    total_points_compared: total_points,
    global_mae,
    global_rmse,
    node_mae_median: median(&node_maes),
    node_rmse_median: median(&node_rmses),
    node_max_abs_p95: percentile(&node_maxes, 0.95),
    node_max_abs_median: median(&node_maxes),
    worst_node,
    worst_t_tbp,
    worst_max_abs_diff,
  };

  Ok((stats, rows))
}

fn print_comparison_table(args: &Args, result: &TestResult) {
  println!("### {}\n", result.name);

  let (stats, rows) = match compute_diff_report(result) {
    Ok(v) => v,
    Err(e) => {
      eprintln!("Error: Cannot compute diff report for {}: {e}", result.name);
      println!();
      return;
    },
  };

  println!(
    "Nodes: expected={}, actual={}, compared={} (missing={}, extra={})",
    stats.n_expected_nodes,
    stats.n_actual_nodes,
    stats.n_compared_nodes,
    stats.n_missing_nodes,
    stats.n_extra_nodes
  );
  println!(
    "Grid: points_per_node={}, total_points_compared={}",
    stats.points_per_node, stats.total_points_compared
  );
  println!();

  let mut shown: Vec<&NodeDiffRow> = rows
    .iter()
    .filter(|r| r.max_abs_diff >= args.abs_diff_threshold)
    .collect();

  let n_above = shown.len();
  if shown.is_empty() {
    shown = rows.iter().take(args.max_rows).collect();
  } else if shown.len() > args.max_rows {
    shown.truncate(args.max_rows);
  }

  println!(
    "| {:<50} | {:>12} | {:>12} | {:>12} | {:>12} | {:>12} |",
    "Node",
    "Max |Δ|",
    "t@Max",
    "MAE",
    "RMSE",
    "p95 |Δ|"
  );
  println!(
    "| {:-<50} | {:-<12} | {:-<12} | {:-<12} | {:-<12} | {:-<12} |",
    "", "", "", "", "", ""
  );

  for row in &shown {
    let short_name = truncate_right_with_ellipsis(&row.node, 50);
    let max_str = row.max_abs_diff.to_significant_digits(4);
    let t_str = row.t_at_max_tbp.to_significant_digits(6);
    let mae_str = row.mae.to_significant_digits(4);
    let rmse_str = row.rmse.to_significant_digits(4);
    let p95_str = row.p95_abs_diff.to_significant_digits(4);
    println!("| {short_name:<50} | {max_str:>12} | {t_str:>12} | {mae_str:>12} | {rmse_str:>12} | {p95_str:>12} |");
  }

  if n_above > 0 {
    println!();
    println!(
      "Showing {}/{} nodes with max |Δ| ≥ {:.6} (cap={})",
      shown.len(),
      n_above,
      args.abs_diff_threshold,
      args.max_rows
    );
  } else {
    println!();
    println!(
      "No nodes exceed max |Δ| ≥ {:.6}; showing worst {} nodes",
      args.abs_diff_threshold,
      shown.len()
    );
  }

  println!();
  println!("#### Snapshot Summary");
  println!();
  println!("| {:<32} | {:>14} |", "Metric", "Value");
  println!("| {:-<32} | {:-<14} |", "", "");
  println!(
    "| {:<32} | {:>14} |",
    "Global MAE (all points)",
    stats.global_mae.to_significant_digits(4)
  );
  println!(
    "| {:<32} | {:>14} |",
    "Global RMSE (all points)",
    stats.global_rmse.to_significant_digits(4)
  );
  println!(
    "| {:<32} | {:>14} |",
    "Node median MAE",
    stats.node_mae_median.to_significant_digits(4)
  );
  println!(
    "| {:<32} | {:>14} |",
    "Node median RMSE",
    stats.node_rmse_median.to_significant_digits(4)
  );
  println!(
    "| {:<32} | {:>14} |",
    "Node median max |Δ|",
    stats.node_max_abs_median.to_significant_digits(4)
  );
  println!(
    "| {:<32} | {:>14} |",
    "Node p95 max |Δ|",
    stats.node_max_abs_p95.to_significant_digits(4)
  );
  println!(
    "| {:<32} | {:>14} |",
    "Worst max |Δ|",
    stats.worst_max_abs_diff.to_significant_digits(4)
  );
  println!(
    "| {:<32} | {:>14} |",
    "Worst t@Max (TBP)",
    stats.worst_t_tbp.to_significant_digits(6)
  );
  println!("| {:<32} | {:>14} |", "Worst node", truncate_right_with_ellipsis(&stats.worst_node, 14));
  println!();
}

fn print_aggregate_statistics(results: &[TestResult]) {
  println!(
    "| {:<15} | {:>8} | {:>12} | {:>12} | {:>12} | {:>12} |",
    "Test", "Nodes", "Global MAE", "Global RMSE", "p95 max|Δ|", "Worst max|Δ|"
  );
  println!(
    "| {:-<15} | {:-<8} | {:-<12} | {:-<12} | {:-<12} | {:-<12} |",
    "", "", "", "", "", ""
  );

  let mut total_abs = 0.0;
  let mut total_sq = 0.0;
  let mut total_points = 0_usize;

  let mut overall_worst_max = -1.0;
  let mut overall_worst_test = String::new();
  let mut overall_worst_node = String::new();
  let mut overall_worst_t_tbp = 0.0;

  for result in results {
    let (stats, _) = match compute_diff_report(result) {
      Ok(v) => v,
      Err(e) => {
        eprintln!("Warning: Skipping aggregate stats for {}: {e}", result.name);
        continue;
      },
    };

    println!(
      "| {:<15} | {:>8} | {:>12} | {:>12} | {:>12} | {:>12} |",
      result.name,
      stats.n_compared_nodes,
      stats.global_mae.to_significant_digits(4),
      stats.global_rmse.to_significant_digits(4),
      stats.node_max_abs_p95.to_significant_digits(4),
      stats.worst_max_abs_diff.to_significant_digits(4)
    );

    total_abs += stats.global_mae * (stats.total_points_compared as f64);
    total_sq += (stats.global_rmse * stats.global_rmse) * (stats.total_points_compared as f64);
    total_points += stats.total_points_compared;

    if stats.worst_max_abs_diff > overall_worst_max {
      overall_worst_max = stats.worst_max_abs_diff;
      overall_worst_test = result.name.clone();
      overall_worst_node = stats.worst_node.clone();
      overall_worst_t_tbp = stats.worst_t_tbp;
    }
  }

  let overall_mae = if total_points == 0 {
    0.0
  } else {
    total_abs / (total_points as f64)
  };
  let overall_rmse = if total_points == 0 {
    0.0
  } else {
    (total_sq / (total_points as f64)).sqrt()
  };

  println!(
    "| {:-<15} | {:-<8} | {:-<12} | {:-<12} | {:-<12} | {:-<12} |",
    "", "", "", "", "", ""
  );
  println!(
    "| {:<15} | {:>8} | {:>12} | {:>12} | {:>12} | {:>12} |",
    "Overall",
    "-",
    overall_mae.to_significant_digits(4),
    overall_rmse.to_significant_digits(4),
    "-",
    overall_worst_max.to_significant_digits(4)
  );
  if overall_worst_max >= 0.0 {
    println!(
      "Worst point overall: test={}, node={}, t_tbp={}, max|Δ|={}",
      overall_worst_test,
      overall_worst_node,
      overall_worst_t_tbp.to_significant_digits(6),
      overall_worst_max.to_significant_digits(4)
    );
  }
  println!();
}

fn write_statistics_to_file(path: &str, results: &[TestResult]) -> Result<(), Report> {
  let mut file = File::create(path)?;

  writeln!(
    file,
    "Test,ComparedNodes,ExpectedNodes,ActualNodes,MissingNodes,ExtraNodes,PointsPerNode,TotalPoints,GlobalMAE,GlobalRMSE,NodeMedianMAE,NodeMedianRMSE,NodeMedianMaxAbs,NodeP95MaxAbs,WorstNode,WorstTBP,WorstMaxAbs"
  )?;

  for result in results {
    let (stats, _) = compute_diff_report(result)?;
    writeln!(
      file,
      "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
      result.name,
      stats.n_compared_nodes,
      stats.n_expected_nodes,
      stats.n_actual_nodes,
      stats.n_missing_nodes,
      stats.n_extra_nodes,
      stats.points_per_node,
      stats.total_points_compared,
      stats.global_mae,
      stats.global_rmse,
      stats.node_mae_median,
      stats.node_rmse_median,
      stats.node_max_abs_median,
      stats.node_max_abs_p95,
      stats.worst_node,
      stats.worst_t_tbp,
      stats.worst_max_abs_diff
    )?;
  }

  Ok(())
}

use clap::Parser;
use ctor::ctor;
use eyre::Report;
use indexmap::IndexMap;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
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

  /// Absolute difference threshold for filtering comparison table rows
  #[arg(long, default_value_t = 0.001)]
  abs_diff_threshold: f64,

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
  expected: IndexMap<String, Array1<f64>>,
  actual: IndexMap<String, Array1<f64>>,
}

#[derive(Clone, Serialize, Deserialize)]
struct DiffStats {
  mean_abs_diff: f64,
  median_abs_diff: f64,
  max_abs_diff: f64,
  min_abs_diff: f64,
  std_dev: f64,
  rmse: f64,
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
        print_comparison_table(&result, args.abs_diff_threshold);
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

fn convert_map_keys_to_names<V>(graph: &GraphAncestral, map: &IndexMap<GraphNodeKey, V>) -> IndexMap<String, V>
where
  V: Clone,
{
  map
    .iter()
    .map(|(key, value)| {
      let node = graph.get_node(*key).unwrap();
      let name = node.read_arc().payload().read_arc().name().unwrap().as_ref().to_owned();
      (name, value.clone())
    })
    .collect()
}

fn run_coalescent_test(args: &Args, snapshot_filename: &str, tc_value: f64) -> Result<TestResult, Report> {
  println!("### Test: Coalescent contributions (Tc={tc_value})");

  let snapshot = load_snapshot(snapshot_filename)?;
  let graph = load_graph(&snapshot)?;
  let tc = Distribution::constant(snapshot.inputs.tc);

  let actuals = compute_coalescent_contributions(&graph, &tc)?;
  let actuals = convert_map_keys_to_names(&graph, &actuals);

  let t_grid = Grid::from_range_n_points(
    snapshot.tbp_grid.start,
    snapshot.tbp_grid.end,
    snapshot.tbp_grid.n_points,
  )?;

  let actuals: IndexMap<String, Array1<f64>> = actuals
    .iter()
    .map(|(name, dist)| {
      let y_resampled = dist.eval_many(&t_grid.to_array()).unwrap();
      (name.clone(), y_resampled)
    })
    .collect();

  if args.verbose {
    println!("  Nodes compared: {}", actuals.len());
    println!("  Grid points: {}", snapshot.tbp_grid.n_points);
  }

  Ok(TestResult {
    name: format!("Tc={tc_value}"),
    expected: snapshot.node_contributions,
    actual: actuals,
  })
}

fn calculate_diff_stats(expected: &IndexMap<String, Array1<f64>>, actual: &IndexMap<String, Array1<f64>>) -> DiffStats {
  let diffs: Vec<f64> = expected
    .keys()
    .filter_map(|key| {
      let exp = expected.get(key)?;
      let act = actual.get(key)?;
      let node_diffs: Vec<f64> = exp.iter().zip(act.iter()).map(|(e, a)| (e - a).abs()).collect();
      let node_max_diff = node_diffs.iter().copied().map(OrderedFloat).max().map(|x| x.0)?;
      Some(node_max_diff)
    })
    .collect();

  if diffs.is_empty() {
    return DiffStats {
      mean_abs_diff: 0.0,
      median_abs_diff: 0.0,
      max_abs_diff: 0.0,
      min_abs_diff: 0.0,
      std_dev: 0.0,
      rmse: 0.0,
    };
  }

  let n = diffs.len() as f64;
  let mean_abs_diff = diffs.iter().sum::<f64>() / n;

  let mut sorted_diffs = diffs.clone();
  sorted_diffs.sort_by_key(|&x| OrderedFloat(x));
  let median_abs_diff = sorted_diffs[sorted_diffs.len() / 2];

  let max_abs_diff = diffs.iter().copied().map(OrderedFloat).max().map_or(0.0, |x| x.0);
  let min_abs_diff = diffs.iter().copied().map(OrderedFloat).min().map_or(0.0, |x| x.0);

  let variance = diffs.iter().map(|&d| (d - mean_abs_diff).powi(2)).sum::<f64>() / n;
  let std_dev = variance.sqrt();

  let squared_diffs: f64 = expected
    .keys()
    .filter_map(|key| {
      let exp = expected.get(key)?;
      let act = actual.get(key)?;
      let sum_sq: f64 = exp.iter().zip(act.iter()).map(|(e, a)| (e - a).powi(2)).sum();
      Some(sum_sq)
    })
    .sum();
  let total_points: usize = expected.values().map(|arr| arr.len()).sum();
  let rmse = (squared_diffs / total_points as f64).sqrt();

  DiffStats {
    mean_abs_diff,
    median_abs_diff,
    max_abs_diff,
    min_abs_diff,
    std_dev,
    rmse,
  }
}

fn print_comparison_table(result: &TestResult, abs_diff_threshold: f64) {
  println!("### {}\n", result.name);

  let stats = calculate_diff_stats(&result.expected, &result.actual);

  let mut rows: Vec<_> = result
    .expected
    .keys()
    .filter_map(|key| {
      let expected = result.expected.get(key)?;
      let actual = result.actual.get(key)?;
      let max_diff = expected
        .iter()
        .zip(actual.iter())
        .map(|(e, a)| (e - a).abs())
        .map(OrderedFloat)
        .max()
        .map(|x| x.0)?;
      let mean_diff = expected
        .iter()
        .zip(actual.iter())
        .map(|(e, a)| (e - a).abs())
        .sum::<f64>()
        / expected.len() as f64;
      Some((key.clone(), max_diff, mean_diff))
    })
    .collect();

  rows.sort_by_key(|(_, max_diff, _)| OrderedFloat(-max_diff));

  let above_threshold_count = rows
    .iter()
    .filter(|(_, max_diff, _)| *max_diff >= abs_diff_threshold)
    .count();

  println!("| {:<50} | {:>12} | {:>12} |", "Node", "Max Diff", "Mean Diff");
  println!("| {:-<50} | {:-<12} | {:-<12} |", "", "", "");

  // Show all rows
  for (node, max_diff, mean_diff) in &rows {
    let short_name = truncate_right_with_ellipsis(node, 50);
    let max_diff_str = max_diff.to_significant_digits(4);
    let mean_diff_str = mean_diff.to_significant_digits(4);
    println!("| {short_name:<50} | {max_diff_str:>12} | {mean_diff_str:>12} |");
  }

  if above_threshold_count > 0 {
    println!();
    println!(
      "({} nodes have significant differences ≥{:.6}, {} have small differences <{:.6})",
      above_threshold_count,
      abs_diff_threshold,
      rows.len() - above_threshold_count,
      abs_diff_threshold
    );
  }

  println!();
  println!("#### Statistics");
  println!();
  println!("| {:<30} | {:>12} |", "Metric", "Value");
  println!("| {:-<30} | {:-<12} |", "", "");
  println!(
    "| {:<30} | {:>12} |",
    "Mean Max Abs Diff (per node)",
    stats.mean_abs_diff.to_significant_digits(4)
  );
  println!(
    "| {:<30} | {:>12} |",
    "Median Max Abs Diff (per node)",
    stats.median_abs_diff.to_significant_digits(4)
  );
  println!(
    "| {:<30} | {:>12} |",
    "Max Abs Diff",
    stats.max_abs_diff.to_significant_digits(4)
  );
  println!(
    "| {:<30} | {:>12} |",
    "Min Max Abs Diff",
    stats.min_abs_diff.to_significant_digits(4)
  );
  println!(
    "| {:<30} | {:>12} |",
    "Standard Deviation",
    stats.std_dev.to_significant_digits(4)
  );
  println!(
    "| {:<30} | {:>12} |",
    "RMSE (all points)",
    stats.rmse.to_significant_digits(4)
  );
  println!();
}

fn print_aggregate_statistics(results: &[TestResult]) {
  println!(
    "| {:<15} | {:>18} | {:>20} | {:>14} | {:>15} | {:>11} |",
    "Test", "Mean Max Diff", "Median Max Diff", "Max Diff", "RMSE", "Std Dev"
  );
  println!(
    "| {:-<15} | {:-<18} | {:-<20} | {:-<14} | {:-<15} | {:-<11} |",
    "", "", "", "", "", ""
  );

  for result in results {
    let stats = calculate_diff_stats(&result.expected, &result.actual);
    println!(
      "| {:<15} | {:>18} | {:>20} | {:>14} | {:>15} | {:>11} |",
      result.name,
      stats.mean_abs_diff.to_significant_digits(4),
      stats.median_abs_diff.to_significant_digits(4),
      stats.max_abs_diff.to_significant_digits(4),
      stats.rmse.to_significant_digits(4),
      stats.std_dev.to_significant_digits(4)
    );
  }

  let overall_mean = results
    .iter()
    .map(|r| calculate_diff_stats(&r.expected, &r.actual).mean_abs_diff)
    .sum::<f64>()
    / results.len() as f64;

  let overall_rmse = results
    .iter()
    .map(|r| calculate_diff_stats(&r.expected, &r.actual).rmse)
    .sum::<f64>()
    / results.len() as f64;

  println!(
    "| {:-<15} | {:-<18} | {:-<20} | {:-<14} | {:-<15} | {:-<11} |",
    "", "", "", "", "", ""
  );
  println!(
    "| {:<15} | {:>18} | {:>20} | {:>14} | {:>15} | {:>11} |",
    "Overall",
    overall_mean.to_significant_digits(4),
    "-",
    "-",
    overall_rmse.to_significant_digits(4),
    "-"
  );
  println!();
}

fn write_statistics_to_file(path: &str, results: &[TestResult]) -> Result<(), Report> {
  use std::fs::File;
  use std::io::Write;

  let mut file = File::create(path)?;

  writeln!(
    file,
    "Test,Mean Max Abs Diff,Median Max Abs Diff,Max Abs Diff,Min Max Abs Diff,Std Dev,RMSE"
  )?;

  for result in results {
    let stats = calculate_diff_stats(&result.expected, &result.actual);
    writeln!(
      file,
      "{},{},{},{},{},{},{}",
      result.name,
      stats.mean_abs_diff,
      stats.median_abs_diff,
      stats.max_abs_diff,
      stats.min_abs_diff,
      stats.std_dev,
      stats.rmse
    )?;
  }

  Ok(())
}

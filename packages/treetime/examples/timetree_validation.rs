use clap::Parser;
use ctor::ctor;
use eyre::Report;
use maplit::btreemap;
use ordered_float::OrderedFloat;
use parking_lot::RwLock;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::io::Write;
use std::path::Path;
use std::sync::Arc;
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::marginal::initialize_marginal;
use treetime::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
use treetime::clock::date_constraints::load_date_constraints;
use treetime::clock::find_best_root::params::BranchPointOptimizationParams;
use treetime::gtr::get_gtr::{JC69Params, jc69};
use treetime::representation::partition::fitch::PartitionFitch;
use treetime::representation::partition::marginal_dense::PartitionMarginalDense;
use treetime::representation::partition::timetree::GraphTimetree;
use treetime::representation::partition::traits::PartitionTimetreeAll;
use treetime::representation::payload::timetree::{EdgeTimetree, NodeTimetree};
use treetime::timetree::inference::backward_pass::propagate_distributions_backward;
use treetime::timetree::inference::forward_pass::propagate_distributions_forward;
use treetime::timetree::inference::runner::{BRANCH_GRID_SIZE, run_timetree};
use treetime::timetree::utils::{
  create_poisson_branch_distributions, extract_node_times, initialize_clock_totals_from_time_distributions,
  initialize_node_divergences,
};
use treetime_io::dates_csv::read_dates;
use treetime_io::fasta::read_many_fasta;
use treetime_io::nwk::nwk_read_str;
use treetime_utils::fmt::string::truncate_right_with_ellipsis;
use treetime_utils::init::global::global_init;
use treetime_utils::io::json::{JsonPretty, json_write_file};

#[ctor]
fn init() {
  global_init();
}

const FIXTURES_DIR: &str = concat!(
  env!("CARGO_MANIFEST_DIR"),
  "/src/timetree/inference/__tests__/__fixtures__"
);

#[derive(Parser, Clone, Serialize, Deserialize)]
#[command(
  name = "timetree-validation",
  about = "Validates TreeTime timetree inference against Python v0 baseline results.\n\
           Loads dataset configurations from golden master fixture files and compares\n\
           Rust v1 inference output against Python v0 baseline node times.",
  version
)]
pub struct Args {
  /// Datasets to validate (comma-separated)
  #[arg(long, value_delimiter = ',', default_value = "flu_h3n2_20,ebola_20")]
  datasets: Vec<String>,

  /// Tests to run per dataset (comma-separated: poisson, marginal_sparse, marginal_dense)
  #[arg(
    long,
    value_delimiter = ',',
    default_value = "poisson,marginal_sparse,marginal_dense"
  )]
  tests: Vec<String>,

  /// Output directory base path
  #[arg(long, default_value = "tmp")]
  output_dir: String,

  /// Absolute difference threshold for filtering comparison table rows
  #[arg(long, default_value_t = 0.001)]
  abs_diff_threshold: f64,

  /// Branch grid size for Poisson distribution discretization
  #[arg(long, default_value_t = BRANCH_GRID_SIZE)]
  branch_grid_size: usize,

  /// Optional statistics output file path
  #[arg(long)]
  stats_output: Option<String>,

  /// Enable detailed progress output
  #[arg(long, short = 'v')]
  verbose: bool,
}

#[derive(Debug, Deserialize)]
struct FixtureOutputs {
  rerooted_tree_nwk: String,
  clock_rate: f64,
  sequence_length: usize,
  poisson: BTreeMap<String, f64>,
  marginal_dense: BTreeMap<String, f64>,
}

#[derive(Debug, Deserialize)]
struct FixtureInput {
  #[allow(dead_code)]
  tree_path: String,
  aln_path: String,
  metadata_path: String,
  name_column: Option<String>,
}

struct DatasetConfig {
  name: String,
  rerooted_tree_nwk: String,
  clock_rate: f64,
  sequence_length: usize,
  metadata_path: String,
  aln_path: String,
  name_column: Option<String>,
  expected_poisson: BTreeMap<String, f64>,
  expected_marginal: BTreeMap<String, f64>,
}

#[derive(Clone, Serialize, Deserialize)]
struct TestResult {
  name: String,
  expected: BTreeMap<String, f64>,
  actual: BTreeMap<String, f64>,
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

  println!("# TreeTime Validation Results\n");

  let configs = load_dataset_configs(&args.datasets)?;

  if args.verbose {
    println!(
      "Datasets: {}\nTests: {}\n",
      args.datasets.join(", "),
      args.tests.join(", ")
    );
  }

  let mut results = Vec::new();

  for config in &configs {
    for test_name in &args.tests {
      let test_lower = test_name.to_lowercase();

      if let Some(reason) = skip_reason(&config.name, &test_lower) {
        println!("### {}/{test_lower}: SKIPPED ({reason})\n", config.name);
        continue;
      }

      let result = match test_lower.as_str() {
        "poisson" => run_poisson_test(config, &args)?,
        "marginal_sparse" => run_marginal_sparse_test(config, &args)?,
        "marginal_dense" => run_marginal_dense_test(config, &args)?,
        _ => {
          eprintln!("Warning: Unknown test '{test_name}', skipping");
          continue;
        },
      };

      print_comparison_table(&result, args.abs_diff_threshold);
      results.push(result);
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

fn load_dataset_configs(datasets: &[String]) -> Result<Vec<DatasetConfig>, Report> {
  let fixtures_dir = Path::new(FIXTURES_DIR);

  let outputs_json = fs::read_to_string(fixtures_dir.join("gm_runner_outputs.json"))?;
  let outputs: BTreeMap<String, FixtureOutputs> = serde_json::from_str(&outputs_json)?;

  let inputs_json = fs::read_to_string(fixtures_dir.join("gm_runner_inputs.json"))?;
  let inputs: BTreeMap<String, FixtureInput> = serde_json::from_str(&inputs_json)?;

  let mut configs = Vec::with_capacity(datasets.len());
  for name in datasets {
    let output = outputs
      .get(name)
      .ok_or_else(|| eyre::eyre!("dataset {name:?} not found in gm_runner_outputs.json"))?;
    let input = inputs
      .get(name)
      .ok_or_else(|| eyre::eyre!("dataset {name:?} not found in gm_runner_inputs.json"))?;

    configs.push(DatasetConfig {
      name: name.clone(),
      rerooted_tree_nwk: output.rerooted_tree_nwk.clone(),
      clock_rate: output.clock_rate,
      sequence_length: output.sequence_length,
      metadata_path: input.metadata_path.clone(),
      aln_path: input.aln_path.clone(),
      name_column: input.name_column.clone(),
      expected_poisson: output.poisson.clone(),
      expected_marginal: output.marginal_dense.clone(),
    });
  }

  Ok(configs)
}

fn skip_reason(dataset: &str, test: &str) -> Option<&'static str> {
  match (dataset, test) {
    ("rsv_a_20" | "tb_20", _) => Some("crashes in distribution functions (zero-length branches)"),
    _ => None,
  }
}

fn dump_graph(graph: &GraphTimetree, output_dir: &str, filename: &str) -> Result<(), Report> {
  let output_path = Path::new(output_dir);
  json_write_file(output_path.join(filename), graph, JsonPretty(true))?;
  Ok(())
}

fn run_poisson_test(config: &DatasetConfig, args: &Args) -> Result<TestResult, Report> {
  let output_dir_str = format!("{}/test_timetree_{}_poisson", args.output_dir, config.name);

  println!("### {}/poisson (branch length input mode)", config.name);

  let graph: GraphTimetree = nwk_read_str(&config.rerooted_tree_nwk)?;
  let dates = read_dates(&config.metadata_path, &config.name_column, &None)?;
  load_date_constraints(&dates, &graph)?;

  create_poisson_branch_distributions(&graph, config.clock_rate, config.sequence_length, args.branch_grid_size)?;
  dump_graph(
    &graph,
    &output_dir_str,
    "001_after_create_poisson_branch_distributions.json",
  )?;

  propagate_distributions_backward(&graph, None)?;
  dump_graph(
    &graph,
    &output_dir_str,
    "002_after_propagate_distributions_backward.json",
  )?;

  propagate_distributions_forward(&graph)?;
  dump_graph(
    &graph,
    &output_dir_str,
    "003_after_propagate_distributions_forward.json",
  )?;

  let actual = extract_node_times(&graph);

  println!("Completed\n");

  Ok(TestResult {
    name: format!("{}/poisson", config.name),
    expected: config.expected_poisson.clone(),
    actual,
  })
}

fn run_marginal_sparse_test(config: &DatasetConfig, args: &Args) -> Result<TestResult, Report> {
  let output_dir_str = format!("{}/test_timetree_{}_marginal_sparse", args.output_dir, config.name);

  println!("### {}/marginal_sparse", config.name);

  let mut graph: GraphTimetree = nwk_read_str(&config.rerooted_tree_nwk)?;
  let dates = read_dates(&config.metadata_path, &config.name_column, &None)?;
  load_date_constraints(&dates, &graph)?;

  let alphabet = Alphabet::default();
  let aln = read_many_fasta(&[config.aln_path.as_str()], &alphabet)?;

  let fitch = PartitionFitch::compress(&graph, 0, alphabet, &aln)?;
  let sparse_partition = Arc::new(RwLock::new(
    fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?,
  ));
  dump_graph(&graph, &output_dir_str, "001_after_compress_sequences.json")?;

  let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>> = vec![sparse_partition];

  initialize_marginal(&graph, &partitions, &aln)?;
  dump_graph(&graph, &output_dir_str, "002_after_run_marginal.json")?;

  initialize_node_divergences(&graph);

  let clock_model = estimate_clock_model_with_reroot(
    &mut graph,
    &ClockParams::default(),
    Some(config.clock_rate),
    true,
    &BranchPointOptimizationParams::default(),
    None,
  )?;
  dump_graph(&graph, &output_dir_str, "003_after_clock_model.json")?;

  initialize_clock_totals_from_time_distributions(&graph)?;
  dump_graph(&graph, &output_dir_str, "004_after_initialize_node_times.json")?;

  run_timetree(&mut graph, &partitions, &clock_model, None, false)?;
  dump_graph(&graph, &output_dir_str, "005_after_run_timetree.json")?;

  let actual = extract_node_times(&graph);

  println!("Completed\n");

  Ok(TestResult {
    name: format!("{}/marginal_sparse", config.name),
    expected: config.expected_marginal.clone(),
    actual,
  })
}

fn run_marginal_dense_test(config: &DatasetConfig, args: &Args) -> Result<TestResult, Report> {
  let output_dir_str = format!("{}/test_timetree_{}_marginal_dense", args.output_dir, config.name);

  println!("### {}/marginal_dense", config.name);

  let mut graph: GraphTimetree = nwk_read_str(&config.rerooted_tree_nwk)?;
  let dates = read_dates(&config.metadata_path, &config.name_column, &None)?;
  load_date_constraints(&dates, &graph)?;

  let alphabet = Alphabet::default();
  let aln = read_many_fasta(&[config.aln_path.as_str()], &alphabet)?;

  let dense_partition = Arc::new(RwLock::new(PartitionMarginalDense {
    index: 0,
    gtr: jc69(JC69Params::default())?,
    alphabet,
    length: config.sequence_length,
    nodes: btreemap! {},
    edges: btreemap! {},
  }));

  let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>> = vec![dense_partition];

  initialize_marginal(&graph, &partitions, &aln)?;
  dump_graph(&graph, &output_dir_str, "001_after_run_marginal.json")?;

  initialize_node_divergences(&graph);

  let clock_model = estimate_clock_model_with_reroot(
    &mut graph,
    &ClockParams::default(),
    Some(config.clock_rate),
    true,
    &BranchPointOptimizationParams::default(),
    None,
  )?;
  dump_graph(&graph, &output_dir_str, "003_after_clock_model.json")?;

  initialize_clock_totals_from_time_distributions(&graph)?;
  dump_graph(&graph, &output_dir_str, "004_after_initialize_node_times.json")?;

  run_timetree(&mut graph, &partitions, &clock_model, None, false)?;
  dump_graph(&graph, &output_dir_str, "005_after_run_timetree.json")?;

  let actual = extract_node_times(&graph);

  println!("Completed\n");

  Ok(TestResult {
    name: format!("{}/marginal_dense", config.name),
    expected: config.expected_marginal.clone(),
    actual,
  })
}

fn calculate_diff_stats(expected: &BTreeMap<String, f64>, actual: &BTreeMap<String, f64>) -> DiffStats {
  let diffs: Vec<f64> = expected
    .keys()
    .filter_map(|key| {
      let exp = expected.get(key)?;
      let act = actual.get(key)?;
      Some((exp - act).abs())
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
      Some((exp - act).powi(2))
    })
    .sum();
  let rmse = (squared_diffs / n).sqrt();

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
    .map(|key| {
      let expected = result.expected.get(key).copied().unwrap_or(0.0);
      let actual = result.actual.get(key).copied().unwrap_or(0.0);
      let diff = (expected - actual).abs();
      (key.clone(), expected, actual, diff)
    })
    .collect();

  rows.sort_by_key(|(_, _, _, diff)| OrderedFloat(-diff));

  let filtered_rows: Vec<_> = rows
    .iter()
    .filter(|(_, _, _, diff)| *diff >= abs_diff_threshold)
    .collect();
  let hidden_count = rows.len() - filtered_rows.len();

  println!(
    "| {:<50} | {:>12} | {:>12} | {:>12} |",
    "Node", "Expected", "Actual", "Abs Diff"
  );
  println!("| {:-<50} | {:-<12} | {:-<12} | {:-<12} |", "", "", "", "");

  for (node, expected, actual, diff) in &filtered_rows {
    let short_name = truncate_right_with_ellipsis(node, 50);
    println!("| {short_name:<50} | {expected:>12.6} | {actual:>12.6} | {diff:>12.6} |");
  }

  if hidden_count > 0 {
    println!();
    println!("({hidden_count} nodes with abs diff < {abs_diff_threshold:.6} are hidden)");
  }

  println!();
  println!("#### Statistics");
  println!();
  println!("| {:<30} | {:>12} |", "Metric", "Value");
  println!("| {:-<30} | {:-<12} |", "", "");
  println!("| {:<30} | {:>12.6} |", "Mean Absolute Difference", stats.mean_abs_diff);
  println!(
    "| {:<30} | {:>12.6} |",
    "Median Absolute Difference", stats.median_abs_diff
  );
  println!("| {:<30} | {:>12.6} |", "Max Absolute Difference", stats.max_abs_diff);
  println!("| {:<30} | {:>12.6} |", "Min Absolute Difference", stats.min_abs_diff);
  println!("| {:<30} | {:>12.6} |", "Standard Deviation", stats.std_dev);
  println!("| {:<30} | {:>12.6} |", "RMSE", stats.rmse);
  println!();
}

fn print_aggregate_statistics(results: &[TestResult]) {
  println!(
    "| {:<35} | {:>15} | {:>17} | {:>14} | {:>10} | {:>11} |",
    "Test", "Mean Abs Diff", "Median Abs Diff", "Max Abs Diff", "RMSE", "Std Dev"
  );
  println!(
    "| {:-<35} | {:-<15} | {:-<17} | {:-<14} | {:-<10} | {:-<11} |",
    "", "", "", "", "", ""
  );

  for result in results {
    let stats = calculate_diff_stats(&result.expected, &result.actual);
    println!(
      "| {:<35} | {:>15.6} | {:>17.6} | {:>14.6} | {:>10.6} | {:>11.6} |",
      result.name, stats.mean_abs_diff, stats.median_abs_diff, stats.max_abs_diff, stats.rmse, stats.std_dev
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
    "| {:-<35} | {:-<15} | {:-<17} | {:-<14} | {:-<10} | {:-<11} |",
    "", "", "", "", "", ""
  );
  println!(
    "| {:<35} | {:>15.6} | {:>17} | {:>14} | {:>10.6} | {:>11} |",
    "Overall", overall_mean, "-", "-", overall_rmse, "-"
  );
  println!();
}

fn write_statistics_to_file(path: &str, results: &[TestResult]) -> Result<(), Report> {
  let mut file = fs::File::create(path)?;

  writeln!(
    file,
    "Test,Mean Abs Diff,Median Abs Diff,Max Abs Diff,Min Abs Diff,Std Dev,RMSE"
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

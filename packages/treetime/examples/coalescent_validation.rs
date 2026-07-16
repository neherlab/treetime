use clap::Parser;
use clap::builder::{PossibleValuesParser, TypedValueParser};
use ctor::ctor;
use eyre::{Context, Report};
use indexmap::IndexMap;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::fs::{DirEntry, read_dir};
use std::io::Write;
use std::path::Path;
use std::sync::LazyLock;
use treetime::clock::date_constraints::load_date_constraints;
use treetime::coalescent::coalescent::compute_coalescent_model;
use treetime::o;
use treetime::partition::timetree::GraphTimetree;
use treetime_distribution::Distribution;
use treetime_graph::node::Named;
use treetime_grid::grid::Grid;
use treetime_io::dates_csv::read_dates;
use treetime_io::nwk::nwk_read_file;
use treetime_utils::array::serde::{array1_from_vec, indexmap_array1_from_map};
use treetime_utils::fmt::string::truncate_right_with_ellipsis;
use treetime_utils::init::clap_styles::styles;
use treetime_utils::init::global::global_init;
use treetime_utils::io::json::json_read_file;

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
#[clap(styles = styles())]
pub struct Args {
  /// Only run snapshots with this virus path. If omitted, runs all discovered snapshots.
  #[arg(long, value_parser = PossibleValuesParser::new(coalescent_possible_virus_paths()))]
  virus_path: Option<String>,

  /// Only run snapshots with these Tc values (comma-separated, no spaces). If omitted, runs all discovered snapshots.
  #[arg(
    long,
    value_delimiter = ',',
    value_parser = PossibleValuesParser::new(coalescent_possible_tc_values())
      .map(|s| s.parse::<f64>().unwrap())
  )]
  tc_values: Option<Vec<f64>>,

  /// List discovered snapshots (virus paths and Tc values) and exit.
  #[arg(long)]
  list_snapshots: bool,

  /// Output directory base path
  #[arg(long, default_value = "tmp")]
  output_dir: String,

  /// Pass/fail threshold for max absolute difference.
  /// Tests with worst error below this threshold pass.
  #[arg(long, default_value_t = 1e-5)]
  pass_threshold: f64,

  /// Maximum number of nodes to show in the per-node table (verbose mode only). Use -1 for unlimited.
  #[arg(long, default_value_t = 10, allow_hyphen_values = true)]
  max_rows: isize,

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
  virus_path: String,
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
  dataset: String,
  tc: f64,
  t_grid_tbp: Array1<f64>,
  expected: IndexMap<String, Array1<f64>>,
  actual: IndexMap<String, Array1<f64>>,
}

impl TestResult {
  fn display_name(&self) -> String {
    format!("{} Tc={}", self.dataset, self.tc)
  }
}

/// Per-node error statistics comparing expected vs actual contribution values.
#[derive(Clone, Serialize, Deserialize)]
struct NodeDiffRow {
  node: String,
  n_points: usize,
  /// Maximum Absolute Error: max|expected[i] - actual[i]| over all grid points
  max_abs_err: f64,
  /// Time (in time-before-present units) where maximum absolute error occurred
  t_at_max: f64,
  /// Mean Absolute Error: (1/n) * Σ|expected[i] - actual[i]|
  mae: f64,
  /// Root Mean Square Error: sqrt((1/n) * Σ(expected[i] - actual[i])^2)
  rmse: f64,
  /// 95th percentile of absolute errors
  p95_abs_err: f64,
}

/// Aggregate error statistics across all compared nodes.
#[derive(Clone, Serialize, Deserialize)]
struct DiffStats {
  n_expected_nodes: usize,
  n_actual_nodes: usize,
  n_compared_nodes: usize,
  n_missing_nodes: usize,
  n_extra_nodes: usize,

  points_per_node: usize,
  total_points_compared: usize,

  /// Global Mean Absolute Error across all points from all nodes
  global_mae: f64,
  /// Global Root Mean Square Error across all points from all nodes
  global_rmse: f64,

  /// Median of per-node MAE values
  node_mae_median: f64,
  /// Median of per-node RMSE values
  node_rmse_median: f64,
  /// 95th percentile of per-node Maximum Absolute Error values
  node_max_abs_err_p95: f64,
  /// Median of per-node Maximum Absolute Error values
  node_max_abs_err_median: f64,

  /// Node with the largest Maximum Absolute Error
  worst_node: String,
  /// Time (in time-before-present units) where worst error occurred
  worst_t: f64,
  /// Largest Maximum Absolute Error across all nodes
  worst_max_abs_err: f64,
}

fn main() -> Result<(), Report> {
  let args = Args::parse();
  let mut results = Vec::new();
  let mut failed_tests = Vec::new();

  let snapshots_dir = Path::new(SNAPSHOTS_DIR);
  let discovered = discover_coalescent_snapshots(snapshots_dir)
    .wrap_err_with(|| format!("When discovering snapshots in {}", snapshots_dir.display()))?;

  if args.list_snapshots {
    print_discovered_snapshots(&discovered);
    return Ok(());
  }

  let selected = select_snapshots_owned(&args, discovered)?;

  for discovered in selected {
    match run_coalescent_test(&discovered.file_name, discovered.snapshot) {
      Ok(result) => {
        results.push(result);
      },
      Err(e) => {
        eprintln!("Warning: Test for snapshot '{}' failed: {e}", discovered.file_name);
      },
    }
  }

  if results.is_empty() {
    eprintln!("Error: No valid tests were run");
    return Ok(());
  }

  let all_passed = compute_pass_status(&args, &results, &mut failed_tests);

  let has_failures = !all_passed;
  let show_details = args.verbose || has_failures;
  let only_failures = !args.verbose;

  if show_details {
    let heading = if args.verbose {
      "## Per-Test Details"
    } else {
      "## Failures"
    };
    println!("{heading}\n");

    for result in &results {
      print_test_details(&args, result, only_failures);
    }
  }

  print_summary_table(&args, &results);

  println!();
  if all_passed {
    println!("✅ All tests passed (Max Err < {})", format_sci(args.pass_threshold));
  } else {
    println!(
      "\x1b[31m❌ {} test(s) exceeded Max Err threshold from --pass-threshold={}:\x1b[0m",
      failed_tests.len(),
      format_sci(args.pass_threshold)
    );
    for test_name in &failed_tests {
      println!("\x1b[31m  - {test_name}\x1b[0m");
    }
  }

  if let Some(stats_output_path) = &args.stats_output {
    write_statistics_to_file(stats_output_path, &results)?;
    if args.verbose {
      println!("\nStatistics written to: {stats_output_path}");
    }
  }

  if !all_passed {
    std::process::exit(1);
  }

  Ok(())
}

fn load_graph(snapshot: &Snapshot) -> Result<GraphTimetree, Report> {
  let fixtures_dir = Path::new(SNAPSHOTS_DIR);
  let tree_path = fixtures_dir.join(&snapshot.inputs.tree_path);
  let metadata_path = fixtures_dir.join(&snapshot.inputs.metadata_path);
  let graph = nwk_read_file(&tree_path)?;
  let dates = read_dates(&metadata_path, &[], &Some(o!("name")), &Some(o!("date")))?;
  load_date_constraints(&dates, &graph)?;
  Ok(graph)
}

fn format_sci(v: f64) -> String {
  if v == 0.0 {
    o!("0")
  } else if v.abs() >= 0.01 && v.abs() < 1000.0 {
    format!("{v:.4}")
  } else {
    format!("{v:.2e}")
  }
}

fn run_coalescent_test(_snapshot_filename: &str, snapshot: Snapshot) -> Result<TestResult, Report> {
  let graph = load_graph(&snapshot)?;
  let tc_value = snapshot.inputs.tc;
  let tc = Distribution::constant(tc_value);

  let model = compute_coalescent_model(&graph, &tc)?;

  let t_grid = Grid::from_range_n_points(
    snapshot.tbp_grid.start,
    snapshot.tbp_grid.end,
    snapshot.tbp_grid.n_points,
  )?;

  let t_grid_tbp = t_grid.to_array();

  let calendar_grid = t_grid_tbp.mapv(|time| snapshot.inputs.present_time - time);
  let actuals: IndexMap<String, Array1<f64>> = graph
    .get_nodes()
    .iter()
    .filter_map(|node| {
      let node = node.read_arc();
      let name = node.payload().read_arc().name()?.as_ref().to_owned();
      Some((name, node.outbound().len()))
    })
    .map(|(name, n_children)| {
      let values = calendar_grid
        .iter()
        .map(|&time| {
          if n_children == 0 {
            Ok(model.leaf_cost(time))
          } else {
            model.internal_cost(time, n_children)
          }
        })
        .collect::<Result<Vec<_>, Report>>()?;
      Ok((name, Array1::from_vec(values)))
    })
    .collect::<Result<IndexMap<String, Array1<f64>>, Report>>()?;

  Ok(TestResult {
    dataset: snapshot.virus_path,
    tc: tc_value,
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

fn print_node_table(args: &Args, rows: &[NodeDiffRow]) {
  let max_rows_limit = if args.max_rows < 0 {
    None
  } else {
    Some(args.max_rows as usize)
  };

  let shown: Vec<&NodeDiffRow> = if let Some(limit) = max_rows_limit {
    rows.iter().take(limit).collect()
  } else {
    rows.iter().collect()
  };

  println!(
    "| {:<3} | {:<40} | {:>10} | {:>10} | {:>10} |",
    "", "Node", "Max Err", "MAE", "RMSE"
  );
  println!("| {:-<3} | {:-<40} | {:-<10}:| {:-<10}:| {:-<10}:|", "", "", "", "", "");

  for row in &shown {
    let short_name = truncate_right_with_ellipsis(&row.node, 40);
    let status = if row.max_abs_err < args.pass_threshold {
      "✅"
    } else {
      "❌"
    };
    println!(
      "| {:<2} | {:<40} | {:>10} | {:>10} | {:>10} |",
      status,
      short_name,
      format_sci(row.max_abs_err),
      format_sci(row.mae),
      format_sci(row.rmse)
    );
  }

  if let Some(limit) = max_rows_limit {
    if rows.len() > limit {
      println!(
        "({} more nodes omitted due to --max-rows={})",
        rows.len() - limit,
        limit
      );
    }
  }
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
  let mut worst_t = 0.0;
  let mut worst_max_abs_err = -1.0;

  for node in &common {
    let exp = &result.expected[node];
    let act = &result.actual[node];
    let n = exp.len();

    let mut abs_errors = Vec::with_capacity(n);
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
      abs_errors.push(abs);
      sum_abs += abs;
      sum_sq += diff * diff;
      if abs > max_abs {
        max_abs = abs;
        max_idx = i;
      }
    }

    abs_errors.sort_by_key(|&x| OrderedFloat(x));
    let mae = sum_abs / (n as f64);
    let rmse = (sum_sq / (n as f64)).sqrt();
    let p95_abs_err = percentile(&abs_errors, 0.95);
    let t_at_max = result.t_grid_tbp[max_idx];

    total_abs += sum_abs;
    total_sq += sum_sq;
    total_points += n;

    if max_abs > worst_max_abs_err {
      worst_max_abs_err = max_abs;
      worst_node = node.clone();
      worst_t = t_at_max;
    }

    rows.push(NodeDiffRow {
      node: node.clone(),
      n_points: n,
      max_abs_err: max_abs,
      t_at_max,
      mae,
      rmse,
      p95_abs_err,
    });
  }

  rows.sort_by_key(|r| OrderedFloat(-r.max_abs_err));

  let mut node_maes: Vec<f64> = rows.iter().map(|r| r.mae).collect();
  node_maes.sort_by_key(|&x| OrderedFloat(x));
  let mut node_rmses: Vec<f64> = rows.iter().map(|r| r.rmse).collect();
  node_rmses.sort_by_key(|&x| OrderedFloat(x));
  let mut node_maxes: Vec<f64> = rows.iter().map(|r| r.max_abs_err).collect();
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
    node_max_abs_err_p95: percentile(&node_maxes, 0.95),
    node_max_abs_err_median: median(&node_maxes),
    worst_node,
    worst_t,
    worst_max_abs_err,
  };

  Ok((stats, rows))
}

fn compute_pass_status(args: &Args, results: &[TestResult], failed_tests: &mut Vec<String>) -> bool {
  let mut all_passed = true;

  for result in results {
    let Ok((stats, _)) = compute_diff_report(result) else {
      continue;
    };

    let passed = stats.worst_max_abs_err < args.pass_threshold;
    if !passed {
      all_passed = false;
      failed_tests.push(result.display_name());
    }
  }

  all_passed
}

fn print_summary_table(args: &Args, results: &[TestResult]) {
  println!("## Summary\n");
  println!(
    "| {:<2} | {:<15} | {:>6} | {:>6} | {:>12} | {:>12} | {:>12} |",
    "", "Dataset", "Tc", "Nodes", "Global MAE", "Global RMSE", "Max Err"
  );
  println!(
    "| {:-<2} | {:-<15} | {:-<6}:| {:-<6}:| {:-<12}:| {:-<12}:| {:-<12}:|",
    "", "", "", "", "", "", ""
  );

  for result in results {
    let (stats, _) = match compute_diff_report(result) {
      Ok(v) => v,
      Err(e) => {
        eprintln!("Warning: Cannot compute stats for {}: {e}", result.display_name());
        continue;
      },
    };

    let passed = stats.worst_max_abs_err < args.pass_threshold;
    let status = if passed { "✅" } else { "❌" };

    println!(
      "| {:^} | {:<15} | {:>6} | {:>6} | {:>12} | {:>12} | {:>12} |",
      status,
      result.dataset,
      result.tc,
      stats.n_compared_nodes,
      format_sci(stats.global_mae),
      format_sci(stats.global_rmse),
      format_sci(stats.worst_max_abs_err)
    );
  }
}

fn print_test_details(args: &Args, result: &TestResult, only_failures: bool) {
  let (stats, rows) = match compute_diff_report(result) {
    Ok(v) => v,
    Err(e) => {
      eprintln!("Error: Cannot compute diff report for {}: {e}", result.display_name());
      return;
    },
  };

  let test_passed = stats.worst_max_abs_err < args.pass_threshold;

  if only_failures && test_passed {
    return;
  }

  let failed_rows: Vec<NodeDiffRow> = if only_failures {
    rows
      .into_iter()
      .filter(|r| r.max_abs_err >= args.pass_threshold)
      .collect()
  } else {
    rows
  };

  if failed_rows.is_empty() {
    return;
  }

  let failed_count = failed_rows
    .iter()
    .filter(|r| r.max_abs_err >= args.pass_threshold)
    .count();

  let display_name = result.display_name();
  if failed_count > 0 {
    println!(
      "### {display_name} ({failed_count} node(s) above threshold from --pass-threshold={})\n",
      format_sci(args.pass_threshold)
    );
  } else {
    println!("### {display_name}\n");
  }

  if !only_failures {
    println!(
      "Nodes: {} compared ({} expected, {} actual)",
      stats.n_compared_nodes, stats.n_expected_nodes, stats.n_actual_nodes
    );
    if stats.n_missing_nodes > 0 || stats.n_extra_nodes > 0 {
      println!("  Missing: {}, Extra: {}", stats.n_missing_nodes, stats.n_extra_nodes);
    }
    println!("Grid: {} points/node\n", stats.points_per_node);
  }

  print_node_table(args, &failed_rows);
  println!();
}

fn write_statistics_to_file(path: &str, results: &[TestResult]) -> Result<(), Report> {
  let mut file = File::create(path)?;

  writeln!(
    file,
    "Test,ComparedNodes,ExpectedNodes,ActualNodes,MissingNodes,ExtraNodes,PointsPerNode,TotalPoints,GlobalMAE,GlobalRMSE,NodeMedianMAE,NodeMedianRMSE,NodeMedianMaxAbsErr,NodeP95MaxAbsErr,WorstNode,WorstT,WorstMaxAbsErr"
  )?;

  for result in results {
    let (stats, _) = compute_diff_report(result)?;
    writeln!(
      file,
      "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
      result.display_name(),
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
      stats.node_max_abs_err_median,
      stats.node_max_abs_err_p95,
      stats.worst_node,
      stats.worst_t,
      stats.worst_max_abs_err
    )?;
  }

  Ok(())
}

const SNAPSHOTS_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/src/coalescent/__tests__/__fixtures__");

static COALESCENT_POSSIBLE_VIRUS_PATHS: LazyLock<Vec<&'static str>> = LazyLock::new(|| {
  let snapshots_dir = Path::new(SNAPSHOTS_DIR);
  let Ok(discovered) = discover_coalescent_snapshots(snapshots_dir) else {
    return Vec::new();
  };

  let mut set: BTreeSet<String> = BTreeSet::new();
  for d in discovered {
    set.insert(d.snapshot.virus_path);
  }

  set
    .into_iter()
    .map(|s| {
      let leaked: &'static str = Box::leak(s.into_boxed_str());
      leaked
    })
    .collect()
});

static COALESCENT_POSSIBLE_TC_VALUES: LazyLock<Vec<&'static str>> = LazyLock::new(|| {
  let snapshots_dir = Path::new(SNAPSHOTS_DIR);
  let Ok(discovered) = discover_coalescent_snapshots(snapshots_dir) else {
    return Vec::new();
  };

  let mut set: BTreeSet<OrderedFloat<f64>> = BTreeSet::new();
  for d in discovered {
    set.insert(OrderedFloat(d.snapshot.inputs.tc));
  }

  set
    .into_iter()
    .map(|v| v.into_inner().to_string())
    .map(|s| {
      let leaked: &'static str = Box::leak(s.into_boxed_str());
      leaked
    })
    .collect()
});

struct DiscoveredSnapshot {
  file_name: String,
  snapshot: Snapshot,
}

fn discover_coalescent_snapshots(snapshots_dir: &Path) -> Result<Vec<DiscoveredSnapshot>, Report> {
  let mut entries: Vec<DirEntry> = read_dir(snapshots_dir)
    .wrap_err_with(|| format!("Cannot read directory {}", snapshots_dir.display()))?
    .collect::<Result<Vec<_>, _>>()?;

  entries.sort_by_key(|e| e.file_name());

  let mut out = Vec::new();
  for entry in entries {
    let file_name = entry.file_name().to_string_lossy().to_string();
    let path = entry.path();
    let is_json = path
      .extension()
      .and_then(|s| s.to_str())
      .is_some_and(|s| s.eq_ignore_ascii_case("json"));

    if !file_name.starts_with("gm_coalescent_") || !is_json {
      continue;
    }
    let snapshot =
      json_read_file::<Snapshot, _>(&path).wrap_err_with(|| format!("When reading snapshot {}", path.display()))?;

    out.push(DiscoveredSnapshot { file_name, snapshot });
  }

  if out.is_empty() {
    eyre::bail!(
      "No coalescent snapshots found in {} (expected files like gm_coalescent_*.json)",
      snapshots_dir.display()
    );
  }

  out.sort_by(|a, b| {
    a.snapshot
      .virus_path
      .cmp(&b.snapshot.virus_path)
      .then_with(|| OrderedFloat(a.snapshot.inputs.tc).cmp(&OrderedFloat(b.snapshot.inputs.tc)))
      .then_with(|| a.file_name.cmp(&b.file_name))
  });

  Ok(out)
}

fn coalescent_possible_virus_paths() -> Vec<&'static str> {
  COALESCENT_POSSIBLE_VIRUS_PATHS.clone()
}

fn coalescent_possible_tc_values() -> Vec<&'static str> {
  COALESCENT_POSSIBLE_TC_VALUES.clone()
}

fn print_discovered_snapshots(discovered: &[DiscoveredSnapshot]) {
  let mut by_virus: BTreeMap<String, BTreeSet<OrderedFloat<f64>>> = BTreeMap::new();
  for d in discovered {
    by_virus
      .entry(d.snapshot.virus_path.clone())
      .or_default()
      .insert(OrderedFloat(d.snapshot.inputs.tc));
  }

  println!("## Discovered coalescent snapshots\n");
  for (virus_path, tc_values) in by_virus {
    let tc_list = tc_values
      .into_iter()
      .map(|v| format_sci(v.into_inner()))
      .collect::<Vec<_>>()
      .join(", ");
    println!("- {virus_path}: Tc=[{tc_list}]");
  }
}

fn select_snapshots_owned(args: &Args, discovered: Vec<DiscoveredSnapshot>) -> Result<Vec<DiscoveredSnapshot>, Report> {
  let (after_virus, available_tc) = {
    let mut out = Vec::new();
    let mut tc_set: BTreeSet<OrderedFloat<f64>> = BTreeSet::new();
    for d in discovered {
      let matches = args.virus_path.as_ref().is_none_or(|v| d.snapshot.virus_path == *v);
      if matches {
        tc_set.insert(OrderedFloat(d.snapshot.inputs.tc));
        out.push(d);
      }
    }
    let tc_list = tc_set.into_iter().map(|v| v.into_inner()).collect::<Vec<_>>();
    (out, tc_list)
  };

  if after_virus.is_empty() {
    if let Some(virus_path) = &args.virus_path {
      eyre::bail!("No snapshots found for --virus-path={virus_path}");
    }
    eyre::bail!("No snapshots found");
  }

  if let Some(tc_values) = &args.tc_values {
    let mut missing = Vec::new();
    for wanted in tc_values {
      let found = available_tc.iter().any(|tc| (*tc - *wanted).abs() < 1e-12);
      if !found {
        missing.push(*wanted);
      }
    }

    if !missing.is_empty() {
      let available = available_tc
        .iter()
        .map(|v| format_sci(*v))
        .collect::<Vec<_>>()
        .join(", ");
      let missing = missing.iter().map(|v| format_sci(*v)).collect::<Vec<_>>().join(", ");
      eyre::bail!("Unknown Tc value(s): {missing}. Available Tc values: [{available}]");
    }
  }

  let selected: Vec<DiscoveredSnapshot> = after_virus
    .into_iter()
    .filter(|d| {
      args
        .tc_values
        .as_ref()
        .is_none_or(|wanted| wanted.iter().any(|v| (d.snapshot.inputs.tc - *v).abs() < 1e-12))
    })
    .collect();

  if selected.is_empty() {
    eyre::bail!("No snapshots selected after applying filters");
  }

  Ok(selected)
}

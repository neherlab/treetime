//! Runtime-parametrized benchmark for the branch-length grid and its long-branch tail extension.
//!
//! Drives the marginal-dense timetree inference over a set of `BranchGridConfig` cases and, for
//! each, dumps the inferred node dates and measures the wall-clock cost. Grid size, base grid
//! extent, and tail parameters are read from the sweep file at runtime, so the whole study runs
//! without recompiling.
//!
//!     ./dev/docker/run ./dev/dev Er branch_grid_tail_benchmark -- \
//!       --cases <sweep.json> --dump-dir <dir> [--out-csv <path>]
//!
//! Two regimes are supported through the sweep file:
//!
//! - Reroot-to-best on v0's rerooted tree, at v0's fitted clock rate (`reroot=true`,
//!   `tree_source="rerooted"`, `compare_v0=true`). Node dates are comparable to the v0 golden
//!   master `gm_runner_outputs.json`, so accuracy against v0 is reported.
//! - Keep-root on the original input tree at a fixed clock rate (`reroot=false`,
//!   `tree_source="original"`, `clock_rate=<x>`, `compare_v0=false`). There is no v0 oracle for
//!   this configuration, so only the dumped node dates and timing are produced; decay on/off/real
//!   effects are computed from the dumps in post-processing.
//!
//! Timing isolates the grid-dependent `run_timetree` call: the grid-independent setup (alignment
//! load, marginal init, clock model) is rebuilt each repeat but excluded from the timer. This is a
//! single inference pass with no ML branch-length optimization, so keep-root node dates isolate the
//! tail effect but do not reproduce the full CLI pipeline's values.

use clap::Parser;
use ctor::ctor;
use eyre::{Report, eyre};
use ordered_float::OrderedFloat;
use parking_lot::RwLock;
use serde::Deserialize;
use std::collections::BTreeMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::marginal::initialize_marginal;
use treetime::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
use treetime::clock::date_constraints::load_date_constraints;
use treetime::clock::find_best_root::params::BranchPointOptimizationParams;
use treetime::gtr::get_gtr::{JC69Params, jc69};
use treetime::partition::marginal_dense::PartitionMarginalDense;
use treetime::partition::timetree::{GraphTimetree, PartitionTimetree, PartitionTimetreeAllVec};
use treetime::timetree::inference::branch_length_likelihood::BranchGridConfig;
use treetime::timetree::inference::runner::run_timetree;
use treetime::timetree::utils::{
  extract_node_times, initialize_clock_totals_from_time_distributions, initialize_node_divergences,
};
use treetime_io::csv::default_name_candidates;
use treetime_io::dates_csv::read_dates;
use treetime_io::fasta::read_many_fasta;
use treetime_io::nwk::nwk_read_str;
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

#[derive(Parser)]
#[command(
  name = "branch_grid_tail_benchmark",
  about = "Sweep branch-length grid/extent/tail parameters, dumping node dates and timing.",
  version
)]
struct Args {
  /// Sweep specification file (dataset, regime, and the list of parameter cases).
  #[arg(long, value_hint = clap::ValueHint::FilePath)]
  cases: PathBuf,

  /// Directory for per-case node-date JSON dumps (`<label>.json`).
  #[arg(long, value_hint = clap::ValueHint::DirPath)]
  dump_dir: Option<PathBuf>,

  /// Optional path for the per-case CSV summary (v0-comparable sweeps only).
  #[arg(long, value_hint = clap::ValueHint::FilePath)]
  out_csv: Option<PathBuf>,

  /// Number of worst-offending nodes to tabulate per case (v0-comparable sweeps only).
  #[arg(long, default_value_t = 5)]
  worst_rows: usize,
}

fn main() -> Result<(), Report> {
  rayon::ThreadPoolBuilder::new()
    .num_threads(1)
    .build_global()
    .expect("rayon global thread pool initialization failed");

  let args = Args::parse();
  let spec: SweepSpec = load_json(&args.cases)?;

  if spec.mode != "marginal_dense" {
    return Err(eyre!(
      "mode {:?} is not supported; only \"marginal_dense\" has a v0 oracle and tail extension",
      spec.mode
    ));
  }

  let golden = load_golden()?;
  let dataset_paths = load_dataset_paths()?;
  let case = golden
    .get(&spec.dataset)
    .ok_or_else(|| eyre!("dataset {:?} not found in gm_runner_outputs.json", spec.dataset))?;
  let paths = dataset_paths
    .get(&spec.dataset)
    .ok_or_else(|| eyre!("dataset {:?} not found in gm_runner_inputs.json", spec.dataset))?;

  let tree_nwk = load_tree(&spec, case, paths)?;
  let clock_rate = spec.clock_rate.unwrap_or(case.clock_rate);

  println!("# Branch grid / extent / tail benchmark\n");
  println!("- dataset:     {}", spec.dataset);
  println!("- regime:      reroot={}, tree={}", spec.reroot, spec.tree_source);
  println!("- clock rate:  {clock_rate:.6} subs/site/year");
  println!("- repeats:     {}, threads: 1 (pinned)\n", spec.repeats);

  if let Some(dir) = &args.dump_dir {
    fs::create_dir_all(dir)?;
  }

  let mut results: Vec<CaseResult> = Vec::with_capacity(spec.cases.len());
  for case_spec in &spec.cases {
    let config = case_spec.to_config();
    let mut timings_s: Vec<f64> = Vec::with_capacity(spec.repeats.max(1));
    let mut actual: BTreeMap<String, f64> = BTreeMap::new();

    // Untimed warmup, then `repeats` timed runs. Setup is rebuilt each run because `run_timetree`
    // mutates the graph in place; only the grid-dependent inference is timed.
    for repeat in 0..=spec.repeats {
      let run = run_marginal_dense_timed(
        &tree_nwk,
        paths,
        case.sequence_length,
        clock_rate,
        spec.reroot,
        spec.no_indels,
        &config,
      )?;
      if repeat > 0 {
        timings_s.push(run.run_timetree_seconds);
      }
      actual = run.node_times;
    }

    if let Some(dir) = &args.dump_dir {
      json_write_file(dir.join(format!("{}.json", case_spec.label)), &actual, JsonPretty(true))?;
    }

    let metrics = spec
      .compare_v0
      .then(|| compute_metrics(&case.marginal_dense, &actual, spec.tolerance));
    let time_min = timings_s.iter().copied().map(OrderedFloat).min().map_or(0.0, |x| x.0);
    let time_median = median(&mut timings_s);

    let result = CaseResult {
      label: case_spec.label.clone(),
      config,
      metrics,
      time_min_s: time_min,
      time_median_s: time_median,
    };
    print_case(&result, args.worst_rows);
    results.push(result);
  }

  if spec.compare_v0 {
    print_summary_table(&results, spec.tolerance);
    if let Some(path) = &args.out_csv {
      write_csv(path, &spec, &results)?;
      println!("\nRaw v0-comparison data written to {}", path.display());
    }
  }
  if let Some(dir) = &args.dump_dir {
    println!("Per-case node dates written to {}/", dir.display());
  }

  Ok(())
}

struct TimedRun {
  node_times: BTreeMap<String, f64>,
  run_timetree_seconds: f64,
}

fn run_marginal_dense_timed(
  tree_nwk: &str,
  paths: &DatasetPaths,
  sequence_length: usize,
  clock_rate: f64,
  reroot: bool,
  no_indels: bool,
  config: &BranchGridConfig,
) -> Result<TimedRun, Report> {
  let mut graph: GraphTimetree = nwk_read_str(tree_nwk)?;

  let metadata_path = project_root().join(&paths.metadata_path);
  let dates = read_dates(&metadata_path, &default_name_candidates(), &paths.name_column, &None)?;
  load_date_constraints(&dates, &graph)?;

  let alphabet = Alphabet::default();
  let aln_path = project_root().join(&paths.aln_path);
  let aln = read_many_fasta(&[&aln_path], &alphabet)?;

  let dense_partition = Arc::new(RwLock::new(PartitionTimetree::Dense(PartitionMarginalDense::new(
    0,
    jc69(JC69Params::default())?,
    alphabet,
    sequence_length,
  ))));
  let partitions: PartitionTimetreeAllVec = vec![dense_partition];

  initialize_marginal(&graph, &partitions, &aln)?.value();
  initialize_node_divergences(&graph)?;

  let clock_model = estimate_clock_model_with_reroot(
    &mut graph,
    &ClockParams::default(),
    Some(clock_rate),
    reroot,
    &BranchPointOptimizationParams::default(),
    None,
  )?;

  initialize_clock_totals_from_time_distributions(&graph)?;

  // Time only the grid-dependent inference: this builds branch-length distributions at the
  // configured resolution/extent and convolves them, which is what the sweep is measuring.
  let start = Instant::now();
  run_timetree(&mut graph, &partitions, &clock_model, None, no_indels, config)?;
  let run_timetree_seconds = start.elapsed().as_secs_f64();

  Ok(TimedRun {
    node_times: extract_node_times(&graph),
    run_timetree_seconds,
  })
}

fn load_tree(spec: &SweepSpec, case: &GoldenOutputs, paths: &DatasetPaths) -> Result<String, Report> {
  match spec.tree_source.as_str() {
    "rerooted" => Ok(case.rerooted_tree_nwk.clone()),
    "original" => {
      let path = project_root().join(&paths.tree_path);
      fs::read_to_string(&path).map_err(|err| eyre!("reading original tree {path:?}: {err}"))
    },
    other => Err(eyre!("tree_source {other:?} is not one of \"rerooted\", \"original\"")),
  }
}

/// Per-node absolute-difference statistics of v1 node dates against the v0 baseline.
fn compute_metrics(expected: &BTreeMap<String, f64>, actual: &BTreeMap<String, f64>, tolerance: f64) -> Metrics {
  let mut diffs: Vec<NodeDiff> = Vec::with_capacity(expected.len());
  let mut missing = 0_usize;
  for (node, &exp) in expected {
    match actual.get(node) {
      Some(&act) => diffs.push(NodeDiff {
        node: node.clone(),
        expected: exp,
        actual: act,
        abs_diff: (exp - act).abs(),
      }),
      None => missing += 1,
    }
  }
  diffs.sort_by_key(|d| OrderedFloat(-d.abs_diff));

  let n = diffs.len();
  let mut abs: Vec<f64> = diffs.iter().map(|d| d.abs_diff).collect();
  let max = diffs.first().map_or(0.0, |d| d.abs_diff);
  let worst_node = diffs.first().map(|d| d.node.clone());
  let mean = if n == 0 {
    0.0
  } else {
    abs.iter().sum::<f64>() / n as f64
  };
  let rms = if n == 0 {
    0.0
  } else {
    (abs.iter().map(|d| d * d).sum::<f64>() / n as f64).sqrt()
  };
  let median = median(&mut abs);
  let p95 = percentile(&abs, 95.0);
  let n_over_tol = diffs.iter().filter(|d| d.abs_diff > tolerance).count();

  Metrics {
    n_nodes: n,
    missing,
    max,
    mean,
    median,
    p95,
    rms,
    worst_node,
    n_over_tol,
    worst_rows: diffs,
  }
}

fn print_case(result: &CaseResult, worst_rows: usize) {
  let c = &result.config;
  println!(
    "## {} (n={}, extent={}x/cap{}, E={}, F={:.0e})\n",
    result.label, c.grid_size, c.peak_extent_multiple, c.max_branch_length, c.tail_max_grid_growth, c.tail_rel_floor
  );
  println!(
    "- run_timetree:   {:.3} s (min) / {:.3} s (median)",
    result.time_min_s, result.time_median_s
  );
  let Some(m) = &result.metrics else {
    println!();
    return;
  };
  println!("- nodes compared: {} (missing {})", m.n_nodes, m.missing);
  println!(
    "- max abs diff:   {:.6} years  (node {})",
    m.max,
    m.worst_node.as_deref().unwrap_or("-")
  );
  println!("- rms / mean:     {:.6} / {:.6} years", m.rms, m.mean);
  println!("- median / p95:   {:.6} / {:.6} years", m.median, m.p95);
  println!("- nodes > tol:    {}\n", m.n_over_tol);

  let shown: Vec<&NodeDiff> = m.worst_rows.iter().take(worst_rows).collect();
  if !shown.is_empty() {
    println!("| {:<40} | {:>13} | {:>13} | {:>12} |", "Node", "v0", "v1", "Abs diff");
    println!("| {:-<40} | {:-<13} | {:-<13} | {:-<12} |", "", "", "", "");
    for d in shown {
      println!(
        "| {:<40} | {:>13.6} | {:>13.6} | {:>12.6} |",
        d.node, d.expected, d.actual, d.abs_diff
      );
    }
    println!();
  }
}

fn print_summary_table(results: &[CaseResult], tolerance: f64) {
  println!("## Summary vs v0 (tolerance {tolerance:.1e})\n");
  println!(
    "| {:<24} | {:>5} | {:>7} | {:>3} | {:>7} | {:>10} | {:>10} | {:>8} | {:>9} |",
    "Case", "n", "extent", "E", "F", "max (y)", "rms (y)", ">tol", "time (s)"
  );
  println!(
    "| {:-<24} | {:-<5} | {:-<7} | {:-<3} | {:-<7} | {:-<10} | {:-<10} | {:-<8} | {:-<9} |",
    "", "", "", "", "", "", "", "", ""
  );
  for r in results {
    let m = r.metrics.as_ref();
    println!(
      "| {:<24} | {:>5} | {:>7} | {:>3} | {:>7.0e} | {:>10.4} | {:>10.4} | {:>8} | {:>9.3} |",
      r.label,
      r.config.grid_size,
      r.config.peak_extent_multiple,
      r.config.tail_max_grid_growth,
      r.config.tail_rel_floor,
      m.map_or(0.0, |m| m.max),
      m.map_or(0.0, |m| m.rms),
      m.map_or(0, |m| m.n_over_tol),
      r.time_median_s
    );
  }
  println!();
}

fn write_csv(path: &Path, spec: &SweepSpec, results: &[CaseResult]) -> Result<(), Report> {
  if let Some(parent) = path.parent() {
    fs::create_dir_all(parent)?;
  }
  let mut file = fs::File::create(path)?;
  writeln!(
    file,
    "dataset,mode,label,grid_size,peak_extent_multiple,max_branch_length,tail_max_grid_growth,\
     tail_rel_floor,repeats,tolerance,nodes,missing,max_abs_diff,mean_abs_diff,median_abs_diff,\
     p95_abs_diff,rms_abs_diff,worst_node,nodes_over_tol,time_min_s,time_median_s"
  )?;
  for r in results {
    let m = r.metrics.as_ref();
    writeln!(
      file,
      "{},{},{},{},{},{},{},{:e},{},{:e},{},{},{},{},{},{},{},{},{},{},{}",
      spec.dataset,
      spec.mode,
      r.label,
      r.config.grid_size,
      r.config.peak_extent_multiple,
      r.config.max_branch_length,
      r.config.tail_max_grid_growth,
      r.config.tail_rel_floor,
      spec.repeats,
      spec.tolerance,
      m.map_or(0, |m| m.n_nodes),
      m.map_or(0, |m| m.missing),
      m.map_or(0.0, |m| m.max),
      m.map_or(0.0, |m| m.mean),
      m.map_or(0.0, |m| m.median),
      m.map_or(0.0, |m| m.p95),
      m.map_or(0.0, |m| m.rms),
      m.and_then(|m| m.worst_node.as_deref()).unwrap_or(""),
      m.map_or(0, |m| m.n_over_tol),
      r.time_min_s,
      r.time_median_s
    )?;
  }
  Ok(())
}

fn median(values: &mut [f64]) -> f64 {
  if values.is_empty() {
    return 0.0;
  }
  values.sort_by_key(|&x| OrderedFloat(x));
  let mid = values.len() / 2;
  if values.len().is_multiple_of(2) {
    f64::midpoint(values[mid - 1], values[mid])
  } else {
    values[mid]
  }
}

fn percentile(abs_diffs: &[f64], pct: f64) -> f64 {
  if abs_diffs.is_empty() {
    return 0.0;
  }
  let mut v: Vec<f64> = abs_diffs.to_vec();
  v.sort_by_key(|&x| OrderedFloat(x));
  let rank = (pct / 100.0 * (v.len() - 1) as f64).round() as usize;
  v[rank.min(v.len() - 1)]
}

fn load_json<T: for<'de> Deserialize<'de>>(path: &Path) -> Result<T, Report> {
  let content = fs::read_to_string(path).map_err(|err| eyre!("reading {path:?}: {err}"))?;
  serde_json::from_str(&content).map_err(|err| eyre!("parsing {path:?}: {err}"))
}

fn load_golden() -> Result<BTreeMap<String, GoldenOutputs>, Report> {
  load_json(&Path::new(FIXTURES_DIR).join("gm_runner_outputs.json"))
}

fn load_dataset_paths() -> Result<BTreeMap<String, DatasetPaths>, Report> {
  load_json(&Path::new(FIXTURES_DIR).join("gm_runner_inputs.json"))
}

fn project_root() -> PathBuf {
  PathBuf::from(env!("CARGO_MANIFEST_DIR"))
    .parent()
    .and_then(Path::parent)
    .map(Path::to_path_buf)
    .expect("failed to resolve project root from CARGO_MANIFEST_DIR")
}

fn default_true() -> bool {
  true
}

fn default_tree_source() -> String {
  "rerooted".to_owned()
}

fn default_tolerance() -> f64 {
  1e-6
}

fn default_peak_extent_multiple() -> f64 {
  5.0
}

fn default_max_branch_length() -> f64 {
  5.0
}

#[derive(Debug, Deserialize)]
struct SweepSpec {
  dataset: String,
  mode: String,
  repeats: usize,
  #[serde(default = "default_true")]
  reroot: bool,
  #[serde(default = "default_tree_source")]
  tree_source: String,
  #[serde(default)]
  clock_rate: Option<f64>,
  #[serde(default = "default_true")]
  compare_v0: bool,
  #[serde(default)]
  no_indels: bool,
  #[serde(default = "default_tolerance")]
  tolerance: f64,
  cases: Vec<CaseSpec>,
}

#[derive(Debug, Deserialize)]
struct CaseSpec {
  label: String,
  grid_size: usize,
  #[serde(default = "default_peak_extent_multiple")]
  peak_extent_multiple: f64,
  #[serde(default = "default_max_branch_length")]
  max_branch_length: f64,
  tail_max_grid_growth: usize,
  tail_rel_floor: f64,
}

impl CaseSpec {
  fn to_config(&self) -> BranchGridConfig {
    BranchGridConfig {
      grid_size: self.grid_size,
      peak_extent_multiple: self.peak_extent_multiple,
      max_branch_length: self.max_branch_length,
      tail_max_grid_growth: self.tail_max_grid_growth,
      tail_rel_floor: self.tail_rel_floor,
    }
  }
}

#[derive(Debug, Deserialize)]
struct GoldenOutputs {
  rerooted_tree_nwk: String,
  clock_rate: f64,
  sequence_length: usize,
  marginal_dense: BTreeMap<String, f64>,
}

#[derive(Debug, Deserialize)]
struct DatasetPaths {
  tree_path: String,
  aln_path: String,
  metadata_path: String,
  #[serde(default)]
  name_column: Option<String>,
}

struct CaseResult {
  label: String,
  config: BranchGridConfig,
  metrics: Option<Metrics>,
  time_min_s: f64,
  time_median_s: f64,
}

struct Metrics {
  n_nodes: usize,
  missing: usize,
  max: f64,
  mean: f64,
  median: f64,
  p95: f64,
  rms: f64,
  worst_node: Option<String>,
  n_over_tol: usize,
  worst_rows: Vec<NodeDiff>,
}

struct NodeDiff {
  node: String,
  expected: f64,
  actual: f64,
  abs_diff: f64,
}

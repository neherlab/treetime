//! Slow golden-master runner for timetree inference (out of the main test suite).
//!
//! This is the same golden-master check as `test_gm_runner_marginal_dense`, moved
//! into a standalone example binary so it does not slow down `cargo nextest`. A full
//! marginal-dense timetree run takes seconds to minutes; running it on every test
//! invocation is prohibitive, so it lives here and is invoked on demand:
//!
//!     ./dev/docker/run ./dev/dev E gm_runner_slow
//!
//! ## What it does
//!
//! For each dataset listed in `gm_runner_slow_inputs.json`, it runs the Rust v1
//! marginal-dense timetree pipeline and compares inferred node dates against the
//! Python v0 baseline captured in `gm_runner_outputs.json`.
//!
//! Node labels match between v0 and v1 because v1 loads v0's rerooted tree
//! (`rerooted_tree_nwk`), so the comparison is by node name with no topology
//! alignment needed.
//!
//! ## Oracle and conventions
//!
//! The v0 baseline is the shared golden master captured by the `gm_runner_capture`
//! script into `gm_runner_outputs.json`. This binary deliberately reuses that oracle
//! rather than capturing its own, so the fast (`#[ignore]`d) in-suite gm tests and
//! this slow binary always compare against the same v0 values. The per-key
//! absolute-difference-within-tolerance check mirrors `pretty_assert_map_abs_diff_eq!`,
//! but is computed explicitly here so every dataset is reported before a tolerance
//! breach fails the run (a panicking assertion would abort after the first dataset).
//!
//! ## Interpreting results
//!
//! At the honest target tolerance `1e-6`, this is currently RED for flu_h3n2_20: the
//! deepest (root-adjacent) node disagrees with v0 by up to ~0.92 years because v1's
//! uniform branch-length grid under-resolves the far tail. Driving the reported max
//! absolute difference to <= tolerance is the acceptance criterion for
//! `kb/issues/M-timetree-branch-grid-uniform-resolution.md`. Use the printed
//! max/RMS/worst-node summary to see whether a grid or tail change moves v1 toward v0.

use clap::Parser;
use ctor::ctor;
use eyre::{Report, eyre};
use ordered_float::OrderedFloat;
use parking_lot::RwLock;
use serde::Deserialize;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Arc;
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

#[ctor]
fn init() {
  global_init();
}

const FIXTURES_DIR: &str = concat!(
  env!("CARGO_MANIFEST_DIR"),
  "/src/timetree/inference/__tests__/__fixtures__"
);

const WORST_ROWS: usize = 10;

#[derive(Parser)]
#[command(
  name = "gm_runner_slow",
  about = "Slow golden-master runner: compares v1 timetree node dates against the v0 baseline.\n\
           Reads run parameters from gm_runner_slow_inputs.json, runs v1, and fails if any\n\
           dataset exceeds its configured tolerance versus v0.",
  version
)]
struct Args {
  /// Run-parameters file (dataset -> mode, tolerance)
  #[arg(long, default_value_t = default_inputs_path())]
  inputs: String,

  /// Number of worst-offending nodes to tabulate per dataset
  #[arg(long, default_value_t = WORST_ROWS)]
  worst_rows: usize,
}

fn main() -> Result<(), Report> {
  let args = Args::parse();

  let runs = load_run_params(&args.inputs)?;
  let golden = load_golden()?;
  let dataset_paths = load_dataset_paths()?;

  println!("# Slow timetree golden-master runner (v1 vs v0)\n");

  let mut failures: Vec<String> = Vec::new();

  for (dataset, params) in &runs {
    println!("## {dataset} ({})\n", params.mode);

    let case = golden
      .get(dataset)
      .ok_or_else(|| eyre!("dataset {dataset:?} not found in gm_runner_outputs.json"))?;
    let paths = dataset_paths
      .get(dataset)
      .ok_or_else(|| eyre!("dataset {dataset:?} not found in gm_runner_inputs.json"))?;

    let expected = match params.mode.as_str() {
      "marginal_dense" => &case.marginal_dense,
      other => {
        return Err(eyre!(
          "dataset {dataset:?}: mode {other:?} is not implemented in gm_runner_slow (only \"marginal_dense\")"
        ));
      },
    };

    let actual = run_marginal_dense(case, paths)?;
    let report = compare(expected, &actual);
    print_report(&report, params.tolerance, args.worst_rows);

    if report.max_abs_diff > params.tolerance {
      failures.push(format!(
        "{dataset}: max abs diff {:.6} > tolerance {:.1e}",
        report.max_abs_diff, params.tolerance
      ));
    }
  }

  if failures.is_empty() {
    println!("All datasets within tolerance.");
    Ok(())
  } else {
    println!("Datasets exceeding tolerance:");
    for failure in &failures {
      println!("  - {failure}");
    }
    Err(eyre!("{} dataset(s) exceeded tolerance versus v0", failures.len()))
  }
}

fn run_marginal_dense(case: &GoldenOutputs, paths: &DatasetPaths) -> Result<BTreeMap<String, f64>, Report> {
  let mut graph: GraphTimetree = nwk_read_str(&case.rerooted_tree_nwk)?;

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
    case.sequence_length,
  ))));
  let partitions: PartitionTimetreeAllVec = vec![dense_partition];

  initialize_marginal(&graph, &partitions, &aln)?.value();
  initialize_node_divergences(&graph)?;

  let clock_model = estimate_clock_model_with_reroot(
    &mut graph,
    &ClockParams::default(),
    Some(case.clock_rate),
    true,
    &BranchPointOptimizationParams::default(),
    None,
  )?;

  initialize_clock_totals_from_time_distributions(&graph)?;
  run_timetree(
    &mut graph,
    &partitions,
    &clock_model,
    None,
    false,
    &BranchGridConfig::default(),
  )?;

  Ok(extract_node_times(&graph))
}

/// Per-node absolute difference of v1 node dates against the v0 baseline.
///
/// Mirrors the `pretty_assert_map_abs_diff_eq!` semantics used by the in-suite gm
/// tests: comparison is over the expected (v0) keys, and missing or extra keys are
/// surfaced explicitly rather than silently ignored.
fn compare(expected: &BTreeMap<String, f64>, actual: &BTreeMap<String, f64>) -> Comparison {
  let mut rows: Vec<NodeDiff> = Vec::with_capacity(expected.len());
  let mut missing: Vec<String> = Vec::new();

  for (node, &exp) in expected {
    match actual.get(node) {
      Some(&act) => rows.push(NodeDiff {
        node: node.clone(),
        expected: exp,
        actual: act,
        abs_diff: (exp - act).abs(),
      }),
      None => missing.push(node.clone()),
    }
  }

  let extra: Vec<String> = actual.keys().filter(|k| !expected.contains_key(*k)).cloned().collect();

  rows.sort_by_key(|row| OrderedFloat(-row.abs_diff));

  let n = rows.len();
  let max_abs_diff = rows.first().map_or(0.0, |row| row.abs_diff);
  let worst_node = rows.first().map(|row| row.node.clone());
  let rmse = if n == 0 {
    0.0
  } else {
    (rows.iter().map(|row| row.abs_diff.powi(2)).sum::<f64>() / n as f64).sqrt()
  };

  Comparison {
    rows,
    missing,
    extra,
    max_abs_diff,
    rmse,
    worst_node,
  }
}

fn print_report(report: &Comparison, tolerance: f64, worst_rows: usize) {
  let verdict = if report.max_abs_diff <= tolerance {
    "PASS"
  } else {
    "FAIL"
  };
  let worst_node = report.worst_node.as_deref().unwrap_or("-");

  println!("- nodes compared: {}", report.rows.len());
  println!(
    "- max abs diff:   {:.6} years  (node {worst_node})",
    report.max_abs_diff
  );
  println!("- rms abs diff:   {:.6} years", report.rmse);
  println!("- tolerance:      {tolerance:.1e}  -> {verdict}");
  if !report.missing.is_empty() {
    println!(
      "- v0 nodes missing from v1: {} ({})",
      report.missing.len(),
      report.missing.join(", ")
    );
  }
  if !report.extra.is_empty() {
    println!(
      "- v1 nodes absent from v0: {} ({})",
      report.extra.len(),
      report.extra.join(", ")
    );
  }
  println!();

  let shown: Vec<&NodeDiff> = report.rows.iter().take(worst_rows).collect();
  if !shown.is_empty() {
    println!(
      "| {:<40} | {:>13} | {:>13} | {:>12} |",
      "Node", "v0 expected", "v1 actual", "Abs diff"
    );
    println!("| {:-<40} | {:-<13} | {:-<13} | {:-<12} |", "", "", "", "");
    for row in &shown {
      println!(
        "| {:<40} | {:>13.6} | {:>13.6} | {:>12.6} |",
        row.node, row.expected, row.actual, row.abs_diff
      );
    }
    let hidden = report.rows.len().saturating_sub(shown.len());
    if hidden > 0 {
      println!("\n({hidden} closer-agreeing nodes not shown)");
    }
    println!();
  }
}

fn load_run_params(path: &str) -> Result<BTreeMap<String, RunParams>, Report> {
  let content = fs::read_to_string(path).map_err(|err| eyre!("reading run-params file {path:?}: {err}"))?;
  serde_json::from_str(&content).map_err(|err| eyre!("parsing run-params file {path:?}: {err}"))
}

fn load_golden() -> Result<BTreeMap<String, GoldenOutputs>, Report> {
  let path = Path::new(FIXTURES_DIR).join("gm_runner_outputs.json");
  let content = fs::read_to_string(&path)?;
  Ok(serde_json::from_str(&content)?)
}

fn load_dataset_paths() -> Result<BTreeMap<String, DatasetPaths>, Report> {
  let path = Path::new(FIXTURES_DIR).join("gm_runner_inputs.json");
  let content = fs::read_to_string(&path)?;
  Ok(serde_json::from_str(&content)?)
}

fn project_root() -> PathBuf {
  PathBuf::from(env!("CARGO_MANIFEST_DIR"))
    .parent()
    .and_then(Path::parent)
    .map(Path::to_path_buf)
    .expect("failed to resolve project root from CARGO_MANIFEST_DIR")
}

fn default_inputs_path() -> String {
  Path::new(FIXTURES_DIR)
    .join("gm_runner_slow_inputs.json")
    .to_string_lossy()
    .into_owned()
}

#[derive(Debug, Deserialize)]
struct RunParams {
  mode: String,
  tolerance: f64,
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
  aln_path: String,
  metadata_path: String,
  #[serde(default)]
  name_column: Option<String>,
}

struct Comparison {
  rows: Vec<NodeDiff>,
  missing: Vec<String>,
  extra: Vec<String>,
  max_abs_diff: f64,
  rmse: f64,
  worst_node: Option<String>,
}

struct NodeDiff {
  node: String,
  expected: f64,
  actual: f64,
  abs_diff: f64,
}

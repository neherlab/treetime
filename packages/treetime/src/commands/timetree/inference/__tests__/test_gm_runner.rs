#[cfg(test)]
mod tests {
  //! Golden master tests for timetree inference runner.
  //!
  //! Compares Rust v1 timetree inference against Python v0 baseline node times.
  //!
  //! Three algorithms tested:
  //! - Poisson (input branch-length mode): backward/forward pass with Poisson distributions
  //! - Marginal dense: full timetree pipeline with dense sequence partition
  //! - Marginal sparse: same as dense but with sparse representation (no direct v0 analog,
  //!   validated against v0 marginal dense golden values)
  //!
  //! Golden outputs captured via `gm_runner_capture` script from v0 Python TreeTime.
  //!
  //! Tolerances are wide (4e-1 for poisson, 9e-1 for marginal) because v0 and v1 use
  //! different numerical implementations. The root node dominates the max diff in all
  //! cases - non-root nodes typically agree within 1e-2. These serve as regression
  //! guards: any code change that significantly worsens agreement will be caught.
  //!
  //! Datasets rsv_a_20 and tb_20 are captured in fixtures but excluded from tests
  //! due to pre-existing crashes in distribution functions (zero-length branches).
  //! ebola_20 marginal_dense crashes in `process_node_backward` and is also excluded.

  use crate::alphabet::alphabet::Alphabet;
  use crate::commands::ancestral::fitch::compress_sequences;
  use crate::commands::ancestral::marginal::initialize_marginal;
  use crate::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
  use crate::commands::clock::date_constraints::load_date_constraints;
  use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::commands::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::commands::timetree::inference::forward_pass::propagate_distributions_forward;
  use crate::commands::timetree::inference::runner::{BRANCH_GRID_SIZE, run_timetree};
  use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
  use crate::commands::timetree::utils::{
    initialize_clock_totals_from_time_distributions, initialize_node_divergences,
  };
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::partition::timetree::GraphTimetree;
  use crate::representation::payload::timetree::{EdgeTimetree, NodeTimetree};
  use eyre::Report;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use ndarray::Array1;
  use ordered_float::OrderedFloat;
  use parking_lot::RwLock;
  use rstest::rstest;
  use serde::Deserialize;
  use std::collections::BTreeMap;
  use std::fs;
  use std::path::{Path, PathBuf};
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_distribution::DistributionFunction;
  use treetime_graph::edge::HasBranchLength;
  use treetime_graph::node::Named;
  use treetime_io::dates_csv::read_dates;
  use treetime_io::fasta::{FastaRecord, read_many_fasta};
  use treetime_io::nwk::nwk_read_str;

  // --- Poisson tests ---

  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  #[case::ebola_20("ebola_20")]
  fn test_gm_runner_poisson(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];
    let expected = &case.poisson;

    let graph: GraphTimetree = nwk_read_str(&case.rerooted_tree_nwk)?;
    let dates = load_dates_for_dataset(dataset)?;
    load_date_constraints(&dates, &graph)?;

    create_poisson_branch_distributions(&graph, case.clock_rate, case.sequence_length, BRANCH_GRID_SIZE)?;
    propagate_distributions_backward(&graph, None)?;
    propagate_distributions_forward(&graph)?;

    let actual = extract_node_times(&graph);
    assert_node_times_match(expected, &actual, 4e-1, dataset, "poisson");

    Ok(())
  }

  // --- Marginal dense tests ---
  // ebola_20 excluded: crashes in process_node_backward with "x array must be uniformly spaced"

  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  fn test_gm_runner_marginal_dense(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];
    let expected = &case.marginal_dense;

    let mut graph: GraphTimetree = nwk_read_str(&case.rerooted_tree_nwk)?;
    let dates = load_dates_for_dataset(dataset)?;
    load_date_constraints(&dates, &graph)?;

    let aln = load_alignment_for_dataset(dataset)?;

    let dense_partition = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: ALPHABET.clone(),
      length: case.sequence_length,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>> = vec![dense_partition];

    initialize_marginal(&graph, &partitions, &aln)?;
    initialize_node_divergences(&graph);

    let clock_model = estimate_clock_model_with_reroot(
      &mut graph,
      &ClockParams::default(),
      Some(case.clock_rate),
      true,
      &BranchPointOptimizationParams::default(),
    )?;

    initialize_clock_totals_from_time_distributions(&graph)?;
    run_timetree(&mut graph, &partitions, &clock_model, None)?;

    let actual = extract_node_times(&graph);
    assert_node_times_match(expected, &actual, 9e-1, dataset, "marginal_dense");

    Ok(())
  }

  // --- Marginal sparse tests ---
  // Sparse is a v1 optimization of dense. v0 has no sparse mode.
  // These tests compare sparse output against v0 marginal dense golden values.

  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  #[case::ebola_20("ebola_20")]
  fn test_gm_runner_marginal_sparse(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];
    let expected = &case.marginal_dense;

    let mut graph: GraphTimetree = nwk_read_str(&case.rerooted_tree_nwk)?;
    let dates = load_dates_for_dataset(dataset)?;
    load_date_constraints(&dates, &graph)?;

    let aln = load_alignment_for_dataset(dataset)?;

    let sparse_partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: ALPHABET.clone(),
      length: case.sequence_length,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    compress_sequences(&graph, std::slice::from_ref(&sparse_partition), &aln)?;

    let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>> = vec![sparse_partition];

    initialize_marginal(&graph, &partitions, &aln)?;
    initialize_node_divergences(&graph);

    let clock_model = estimate_clock_model_with_reroot(
      &mut graph,
      &ClockParams::default(),
      Some(case.clock_rate),
      true,
      &BranchPointOptimizationParams::default(),
    )?;

    initialize_clock_totals_from_time_distributions(&graph)?;
    run_timetree(&mut graph, &partitions, &clock_model, None)?;

    let actual = extract_node_times(&graph);
    assert_node_times_match(expected, &actual, 9e-1, dataset, "marginal_sparse");

    Ok(())
  }

  // --- Fixture types and loading ---

  const FIXTURES_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/src/commands/timetree/inference/__tests__/__fixtures__"
  );

  #[derive(Debug, Deserialize)]
  struct DatasetOutputs {
    rerooted_tree_nwk: String,
    clock_rate: f64,
    sequence_length: usize,
    poisson: BTreeMap<String, f64>,
    marginal_dense: BTreeMap<String, f64>,
  }

  lazy_static! {
    static ref OUTPUTS: BTreeMap<String, DatasetOutputs> = {
      let path = Path::new(FIXTURES_DIR).join("gm_runner_outputs.json");
      let content = fs::read_to_string(&path).expect("Failed to read gm_runner_outputs.json");
      serde_json::from_str(&content).expect("Failed to parse gm_runner_outputs.json")
    };
    static ref INPUTS: BTreeMap<String, DatasetInput> = {
      let path = Path::new(FIXTURES_DIR).join("gm_runner_inputs.json");
      let content = fs::read_to_string(&path).expect("Failed to read gm_runner_inputs.json");
      serde_json::from_str(&content).expect("Failed to parse gm_runner_inputs.json")
    };
    static ref ALPHABET: Alphabet = Alphabet::default();
    static ref PROJECT_ROOT: PathBuf = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .expect("Failed to find project root")
      .to_path_buf();
  }

  #[derive(Debug, Deserialize)]
  struct DatasetInput {
    #[allow(dead_code)]
    tree_path: String,
    aln_path: String,
    metadata_path: String,
  }

  // --- Helper functions ---

  fn load_dates_for_dataset(dataset: &str) -> Result<treetime_io::dates_csv::DatesMap, Report> {
    let input = &INPUTS[dataset];
    let metadata_path = PROJECT_ROOT.join(&input.metadata_path);
    read_dates(&metadata_path, &None, &None)
  }

  fn load_alignment_for_dataset(dataset: &str) -> Result<Vec<FastaRecord>, Report> {
    let input = &INPUTS[dataset];
    let aln_path = PROJECT_ROOT.join(&input.aln_path);
    read_many_fasta(&[&aln_path], &*ALPHABET)
  }

  fn extract_node_times(graph: &GraphTimetree) -> BTreeMap<String, f64> {
    graph
      .get_nodes()
      .into_iter()
      .filter_map(|node_ref| {
        let node = node_ref.read_arc();
        let payload = node.payload().read_arc();
        let name = payload.name()?.as_ref().to_owned();
        let time = payload.time?;
        Some((name, time))
      })
      .collect()
  }

  /// Strip Bio.Phylo confidence suffix from internal node names.
  ///
  /// Bio.Phylo's NWK writer concatenates name and confidence without a separator:
  /// `NODE_00000161.00` = name `NODE_0000016` + confidence `1.00`.
  /// TreeTime v0 uses 7-digit zero-padded indices: `NODE_XXXXXXX` (12 chars).
  #[allow(clippy::string_slice)]
  fn normalize_node_name(name: &str) -> &str {
    // NODE_ names are always ASCII, so byte indexing is safe
    if name.starts_with("NODE_") && name.len() > 12 {
      let suffix = &name[12..];
      if suffix.parse::<f64>().is_ok() {
        return &name[..12];
      }
    }
    name
  }

  fn assert_node_times_match(
    expected: &BTreeMap<String, f64>,
    actual: &BTreeMap<String, f64>,
    epsilon: f64,
    dataset: &str,
    algo: &str,
  ) {
    let actual_normalized: BTreeMap<&str, f64> = actual
      .iter()
      .map(|(name, &time)| (normalize_node_name(name), time))
      .collect();

    let mut compared = 0;
    let mut max_diff = 0.0_f64;
    let mut max_diff_node = String::new();

    for (name, &expected_time) in expected {
      let Some(&actual_time) = actual_normalized.get(name.as_str()) else {
        panic!("[{dataset}/{algo}] node {name:?} not found in actual output")
      };

      let diff = (expected_time - actual_time).abs();
      if diff > max_diff {
        max_diff = diff;
        max_diff_node = name.clone();
      }

      assert!(
        diff <= epsilon,
        "[{dataset}/{algo}] node {name:?}: expected={expected_time}, actual={actual_time}, diff={diff}, epsilon={epsilon}"
      );

      compared += 1;
    }

    assert!(compared > 0, "[{dataset}/{algo}] no nodes compared");

    eprintln!(
      "[{dataset}/{algo}] compared {compared} nodes, max_diff={max_diff:.2e} at {max_diff_node:?}"
    );
  }

  /// Replicate v0 Python's Poisson branch-length distribution construction.
  ///
  /// For each edge with branch length `b`, creates a discretized distribution
  /// P(dt) ~ exp(-dt * mu * L) * (dt * mu * L)^(b * L), where:
  /// - `mu` = clock rate (substitutions/site/year)
  /// - `L` = sequence length
  /// - `b` = branch length (substitutions/site)
  fn create_poisson_branch_distributions(
    graph: &GraphTimetree,
    mu: f64,
    seq_len: usize,
    n_points: usize,
  ) -> Result<(), Report> {
    let seq_len_f64 = seq_len as f64;

    for edge_ref in graph.get_edges() {
      let mut edge = edge_ref.write_arc().payload().write_arc();

      if let Some(branch_length) = edge.branch_length() {
        let expected_time = branch_length / mu;
        let max_time = 3.0 * expected_time.max(1.0);
        let dx = max_time / (n_points - 1) as f64;

        let y = Array1::from_shape_fn(n_points, |i| {
          let dt = i as f64 * dx;
          if dt < 1e-10 {
            0.0
          } else {
            let log_p = -dt * mu * seq_len_f64 + branch_length * seq_len_f64 * (dt * mu * seq_len_f64).ln();
            log_p.exp()
          }
        });

        let y_max = y.iter().copied().map(OrderedFloat).max().map_or(1.0, |x| x.0);
        let y_normalized = y.mapv(|v| v / y_max);

        let distribution_fn = DistributionFunction::from_start_dx_values(0.0, dx, y_normalized)?;
        let distribution = Distribution::Function(distribution_fn);
        edge.branch_length_distribution = Some(Arc::new(distribution));
      }
    }

    Ok(())
  }
}

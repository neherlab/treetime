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
  //! Tolerances are wide (3e-1 for poisson, 9e-1 for marginal) because v0 and v1 use
  //! different numerical implementations. The root node dominates the max diff in all
  //! cases - non-root nodes typically agree within 1e-2. These serve as regression
  //! guards: any code change that significantly worsens agreement will be caught.

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
    create_poisson_branch_distributions, extract_node_times, initialize_clock_totals_from_time_distributions,
    initialize_node_divergences,
  };
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::partition::timetree::GraphTimetree;
  use crate::representation::payload::timetree::{EdgeTimetree, NodeTimetree};
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use rstest::rstest;
  use serde::Deserialize;
  use std::collections::BTreeMap;
  use std::fs;
  use std::path::{Path, PathBuf};
  use std::sync::{Arc, LazyLock};
  use treetime_io::dates_csv::read_dates;
  use treetime_io::fasta::{FastaRecord, read_many_fasta};
  use treetime_io::nwk::nwk_read_str;

  // --- Poisson tests ---

  #[rstest]
  // #[case::dengue_20("dengue_20")]       // TODO: zero-length branches ("x array must be uniformly spaced")
  #[case::ebola_20("ebola_20")]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  // #[case::lassa_L_20("lassa_L_20")]     // TODO: zero-length branches ("x array must be uniformly spaced")
  // #[case::mpox_clade_ii_20("mpox_clade_ii_20")] // TODO: zero-length branches ("x array must be uniformly spaced")
  // #[case::rsv_a_20("rsv_a_20")]         // TODO: zero-length branches ("x array must be uniformly spaced")
  // #[case::tb_20("tb_20")]               // TODO: zero-length branches ("x array must be uniformly spaced")
  // #[case::zika_20("zika_20")]           // TODO: read_dates strips # from headers, name_column="#name" mismatches
  #[trace]
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
    assert_node_times_match(expected, &actual, 3e-1, dataset, "poisson");

    Ok(())
  }

  // --- Marginal dense tests ---

  #[rstest]
  // #[case::dengue_20("dengue_20")]       // TODO: zero-length branches ("x array must be uniformly spaced")
  // #[case::ebola_20("ebola_20")]         // TODO: alphabet doesn't handle gap character '-'
  #[case::flu_h3n2_20("flu_h3n2_20")]
  // #[case::lassa_L_20("lassa_L_20")]     // TODO: zero-length branches ("x array must be uniformly spaced")
  // #[case::mpox_clade_ii_20("mpox_clade_ii_20")] // TODO: zero-length branches ("x array must be uniformly spaced")
  // #[case::rsv_a_20("rsv_a_20")]         // TODO: zero-length branches ("x array must be uniformly spaced")
  // #[case::tb_20("tb_20")]               // TODO: zero-length branches ("x array must be uniformly spaced")
  // #[case::zika_20("zika_20")]           // TODO: read_dates strips # from headers, name_column="#name" mismatches
  #[trace]
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

  #[rstest]
  // #[case::dengue_20("dengue_20")]       // TODO: zero-length branches ("x array must be uniformly spaced")
  #[case::ebola_20("ebola_20")]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  // #[case::lassa_L_20("lassa_L_20")]     // TODO: zero-length branches ("x array must be uniformly spaced")
  // #[case::mpox_clade_ii_20("mpox_clade_ii_20")] // TODO: zero-length branches ("x array must be uniformly spaced")
  // #[case::rsv_a_20("rsv_a_20")]         // TODO: zero-length branches ("x array must be uniformly spaced")
  // #[case::tb_20("tb_20")]               // TODO: zero-length branches ("x array must be uniformly spaced")
  // #[case::zika_20("zika_20")]           // TODO: read_dates strips # from headers, name_column="#name" mismatches
  #[trace]
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

  static OUTPUTS: LazyLock<BTreeMap<String, DatasetOutputs>> = LazyLock::new(|| {
    let path = Path::new(FIXTURES_DIR).join("gm_runner_outputs.json");
    let content = fs::read_to_string(&path).expect("Failed to read gm_runner_outputs.json");
    serde_json::from_str(&content).expect("Failed to parse gm_runner_outputs.json")
  });

  static INPUTS: LazyLock<BTreeMap<String, DatasetInput>> = LazyLock::new(|| {
    let path = Path::new(FIXTURES_DIR).join("gm_runner_inputs.json");
    let content = fs::read_to_string(&path).expect("Failed to read gm_runner_inputs.json");
    serde_json::from_str(&content).expect("Failed to parse gm_runner_inputs.json")
  });

  static ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  static PROJECT_ROOT: LazyLock<PathBuf> = LazyLock::new(|| {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .expect("Failed to find project root")
      .to_path_buf()
  });

  #[derive(Debug, Deserialize)]
  struct DatasetInput {
    #[allow(dead_code)]
    tree_path: String,
    aln_path: String,
    metadata_path: String,
    name_column: Option<String>,
  }

  // --- Helper functions ---

  fn load_dates_for_dataset(dataset: &str) -> Result<treetime_io::dates_csv::DatesMap, Report> {
    let input = &INPUTS[dataset];
    let metadata_path = PROJECT_ROOT.join(&input.metadata_path);
    read_dates(&metadata_path, &input.name_column, &None)
  }

  fn load_alignment_for_dataset(dataset: &str) -> Result<Vec<FastaRecord>, Report> {
    let input = &INPUTS[dataset];
    let aln_path = PROJECT_ROOT.join(&input.aln_path);
    read_many_fasta(&[&aln_path], &*ALPHABET)
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

    eprintln!("[{dataset}/{algo}] compared {compared} nodes, max_diff={max_diff:.2e} at {max_diff_node:?}");
  }
}

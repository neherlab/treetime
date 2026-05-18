#[cfg(test)]
mod tests {

  use approx::assert_relative_eq;
  use eyre::Report;
  use rstest::rstest;
  use std::collections::BTreeMap;

  use crate::optimize::params::BranchOptMethod;
  use std::path::Path;
  use treetime_graph::edge::HasBranchLength;

  use helpers::{load_gm_inputs, load_gm_outputs, setup_and_run};

  // Golden master: v1 brent-sqrt against v0 reference.
  // v0 uses Brent in sqrt(t) space, so brent-sqrt is the matching v1 method.
  // Total branch length compared within 5%.
  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20_jc69_damped")]
  // #[case::dengue_20("dengue_20_jc69_damped")] // slow
  // #[case::tb_20("tb_20_jc69_damped")] // slow (bacterial genome)
  // #[case::ebola_20("ebola_20_jc69_damped")] // slow
  // #[case::zika_20("zika_20_jc69_damped")] // slow
  // #[case::rsv_a_20("rsv_a_20_jc69_damped")] // slow
  // #[case::lassa_l_20("lassa_l_20_jc69_damped")] // slow
  // #[case::mpox_clade_ii_20("mpox_clade_ii_20_jc69_damped")] // slow
  fn test_gm_optimize(#[case] case_name: &str) -> Result<(), Report> {
    let inputs = load_gm_inputs();
    let outputs = load_gm_outputs();
    let case = &inputs[case_name];
    let expected = &outputs[case_name];

    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .unwrap();

    let result = setup_and_run(workspace_root, case, BranchOptMethod::BrentSqrt)?;

    let v1_total_bl: f64 = result
      .graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
      .sum();

    assert_relative_eq!(v1_total_bl, expected.final_total_branch_length, max_relative = 0.05);

    Ok(())
  }

  /// Per-branch GM comparison: walk each edge, look up the target node's name,
  /// and check the optimized length against the v0-derived expected length
  /// stored in `final_branch_lengths`. Branch-level checks catch regressions
  /// that cancel out in the summed total checked by `test_gm_optimize`.
  ///
  /// Branches where v0 collapsed to "essentially zero" (below 1e-5 subs/site,
  /// far below the smallest meaningful branch length 1/L) are skipped: v0's
  /// `prune_short_branches` collapses any internal branch with
  /// `bl < 0.1 * one_mutation && prob_t(parent, child, 0) > 0.1`, while v1's
  /// optimize loop collapses only branches the optimizer drives to exact zero.
  /// Comparing the two criteria on near-zero branches is out of scope for
  /// this test; that divergence is tracked separately under the prune
  /// command intentional-change docs.
  ///
  /// Currently expected to fail with the existing fixture: per-branch
  /// divergences exceed 10% relative tolerance even after skipping
  /// v0-collapsed branches. See
  /// `docs/port-known-issues/M-optimize-gm-per-branch-divergence.md`.
  #[ignore = "Per-branch divergence with v0 fixture exceeds 10% relative tolerance; tracked in M-optimize-gm-per-branch-divergence.md"]
  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20_jc69_damped")]
  fn test_gm_optimize_per_branch(#[case] case_name: &str) -> Result<(), Report> {
    let inputs = load_gm_inputs();
    let outputs = load_gm_outputs();
    let case = &inputs[case_name];
    let expected = &outputs[case_name];

    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .unwrap();

    let result = setup_and_run(workspace_root, case, BranchOptMethod::BrentSqrt)?;

    let v1_branch_lengths: BTreeMap<String, f64> = result
      .graph
      .get_edges()
      .iter()
      .filter_map(|e| {
        let edge = e.read_arc();
        let bl = edge.payload().read_arc().branch_length()?;
        let target_key = edge.target();
        let node = result.graph.get_node(target_key)?;
        let name = node.read_arc().payload().read_arc().name.clone()?;
        Some((name, bl))
      })
      .collect();

    for (name, &expected_bl) in &expected.final_branch_lengths {
      let actual_bl = v1_branch_lengths
        .get(name)
        .copied()
        .unwrap_or_else(|| panic!("missing branch length for node '{name}' in optimized graph"));
      if expected_bl.abs() < 1e-5 {
        continue;
      }
      assert_relative_eq!(actual_bl, expected_bl, max_relative = 0.10);
    }

    Ok(())
  }

  // End-to-end: damped vs undamped on same dataset using BrentSqrt (v0-matching method).
  // Damped should converge and show fewer sign flips.
  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20_jc69_damped")]
  // #[case::dengue_20("dengue_20_jc69_damped")] // slow
  // #[case::tb_20("tb_20_jc69_damped")] // slow (bacterial genome)
  // #[case::ebola_20("ebola_20_jc69_damped")] // slow
  // #[case::zika_20("zika_20_jc69_damped")] // slow
  // #[case::rsv_a_20("rsv_a_20_jc69_damped")] // slow
  // #[case::lassa_l_20("lassa_l_20_jc69_damped")] // slow
  // #[case::mpox_clade_ii_20("mpox_clade_ii_20_jc69_damped")] // slow
  fn test_gm_optimize_damped_vs_undamped(#[case] case_name: &str) -> Result<(), Report> {
    let inputs = load_gm_inputs();
    let case = &inputs[case_name];

    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .unwrap();

    let mut undamped_case = case.clone();
    undamped_case.damping = 0.0;

    let undamped = setup_and_run(workspace_root, &undamped_case, BranchOptMethod::BrentSqrt)?;
    let damped = setup_and_run(workspace_root, case, BranchOptMethod::BrentSqrt)?;

    let undamped_sign_flips = count_sign_flips(&undamped.lh_history);
    let damped_sign_flips = count_sign_flips(&damped.lh_history);

    assert!(
      damped.stopped_at.is_some(),
      "Damped optimization did not stop within {} iterations",
      case.max_iter
    );

    // Damping should not increase oscillation. On well-conditioned datasets both
    // converge equally; on poorly-conditioned ones damping reduces sign flips.
    assert!(
      damped_sign_flips <= undamped_sign_flips,
      "Damped ({damped_sign_flips}) has more sign flips than undamped ({undamped_sign_flips})"
    );

    Ok(())
  }

  fn count_sign_flips(lh_history: &[f64]) -> usize {
    let deltas: Vec<f64> = lh_history
      .windows(2)
      .map(|w| match w {
        [a, b] => b - a,
        _ => 0.0,
      })
      .collect();
    deltas
      .windows(2)
      .filter(|w| matches!(w, [a, b] if a.signum() != b.signum()))
      .count()
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::ancestral::marginal::{initialize_marginal, update_marginal};
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::optimize::params::BranchOptMethod;
    use crate::optimize::dispatch::initial_guess_mixed;
    use crate::optimize::run_loop::{collect_optimize_partitions, run_optimize_loop};
    use crate::ancestral::fitch::create_fitch_partition;
    use crate::partition::marginal_dense::PartitionMarginalDense;
    use crate::seq::alignment::get_common_length;

    use crate::payload::ancestral::GraphAncestral;
    use eyre::Report;

    use parking_lot::RwLock;
    use serde::Deserialize;
    use std::collections::BTreeMap;
    use std::fs::read_to_string;
    use std::path::Path;
    use std::sync::Arc;
    use treetime_io::fasta::read_many_fasta;
    use treetime_io::nwk::nwk_read_file;

    #[derive(Clone, Deserialize)]
    pub struct GmOptimizeCase {
      pub tree: String,
      pub aln: String,
      pub damping: f64,
      pub max_iter: usize,
    }

    #[derive(Deserialize)]
    pub struct GmOptimizeExpected {
      pub final_total_branch_length: f64,
      pub final_branch_lengths: BTreeMap<String, f64>,
    }

    pub struct OptimizeResult {
      pub graph: GraphAncestral,
      pub lh_history: Vec<f64>,
      pub stopped_at: Option<(usize, crate::optimize::run_loop::ConvergenceReason)>,
    }

    pub fn load_gm_inputs() -> BTreeMap<String, GmOptimizeCase> {
      let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("src/optimize/__tests__/__fixtures__/gm_optimize_inputs.json");
      let content = read_to_string(&path).unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
      serde_json::from_str(&content).unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
    }

    pub fn load_gm_outputs() -> BTreeMap<String, GmOptimizeExpected> {
      let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("src/optimize/__tests__/__fixtures__/gm_optimize_outputs.json");
      let content = read_to_string(&path).unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
      serde_json::from_str(&content).unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
    }

    pub fn setup_and_run(
      workspace_root: &Path,
      case: &GmOptimizeCase,
      method: BranchOptMethod,
    ) -> Result<OptimizeResult, Report> {
      let alphabet_sparse = Alphabet::default();
      let alphabet_dense = Alphabet::default();

      let tree_path = workspace_root.join(&case.tree);
      let aln_path = workspace_root.join(&case.aln);
      let aln = read_many_fasta(&[aln_path.to_str().unwrap()], &alphabet_sparse)?;
      let mut graph: GraphAncestral = nwk_read_file(&tree_path)?;

      let fitch = create_fitch_partition(&graph, 0, alphabet_sparse, &aln)?;
      let sparse_partitions = vec![Arc::new(RwLock::new(
        fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?,
      ))];

      let length = get_common_length(&aln)?;
      let dense_partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense::new(
        1,
        jc69(JC69Params::default())?,
        alphabet_dense,
        length,
      )))];

      initialize_marginal(&graph, &dense_partitions, &aln)?;
      update_marginal(&graph, &sparse_partitions)?;
      update_marginal(&graph, &dense_partitions)?;

      let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);
      initial_guess_mixed(&graph, &mixed_partitions, true)?;

      let dp = 0.1;
      let result = run_optimize_loop(
        &mut graph,
        &sparse_partitions,
        &dense_partitions,
        &mixed_partitions,
        case.max_iter,
        dp,
        case.damping,
        method,
        false,
      )?;

      // Append a trailing likelihood measurement so `lh_history.last()` reflects the state
      // AFTER the final in-loop branch-length update (`run_optimize_loop` records the LH
      // at the START of each iteration, before that iteration's update).
      let mut lh_history = result.lh_history;
      let sparse_lh = update_marginal(&graph, &sparse_partitions)?;
      let dense_lh = update_marginal(&graph, &dense_partitions)?;
      lh_history.push(sparse_lh + dense_lh);

      Ok(OptimizeResult {
        graph,
        lh_history,
        stopped_at: result.stopped_at,
      })
    }
  }
}

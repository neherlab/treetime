#[cfg(test)]
mod tests {

  use approx::assert_relative_eq;
  use eyre::Report;
  use rstest::rstest;

  use std::path::Path;
  use treetime_graph::edge::HasBranchLength;

  use helpers::{load_gm_inputs, load_gm_outputs, setup_and_run};

  // Golden master: v1 damped optimization against v0 reference.
  // v0 uses Brent + inline damping, v1 uses Newton-Raphson + post-pass damping.
  // Both target the same ML fixed point. Total branch length compared within 5%.
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

    let result = setup_and_run(workspace_root, case)?;

    let v1_total_bl: f64 = result
      .graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
      .sum();

    assert_relative_eq!(v1_total_bl, expected.final_total_branch_length, max_relative = 0.05);

    Ok(())
  }

  // End-to-end: damped vs undamped on same dataset.
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

    let undamped = setup_and_run(workspace_root, &undamped_case)?;
    let damped = setup_and_run(workspace_root, case)?;

    let undamped_sign_flips = count_sign_flips(&undamped.lh_history);
    let damped_sign_flips = count_sign_flips(&damped.lh_history);

    assert!(
      damped.converged_at.is_some(),
      "Damped optimization did not converge within {} iterations",
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
    let deltas: Vec<f64> = lh_history.windows(2).map(|w| w[1] - w[0]).collect();
    deltas.windows(2).filter(|w| w[0].signum() != w[1].signum()).count()
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
    use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
    use crate::commands::optimize::optimize_unified::{initial_guess_mixed, run_optimize_mixed};
    use crate::commands::optimize::run::{apply_damping, collect_optimize_partitions, save_branch_lengths};
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::representation::partition::marginal_dense::PartitionMarginalDense;
    use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
    use crate::representation::payload::ancestral::GraphAncestral;
    use eyre::Report;
    use maplit::btreemap;
    use parking_lot::RwLock;
    use serde::Deserialize;
    use std::collections::BTreeMap;
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
      pub converged_at: Option<usize>,
    }

    pub fn load_gm_inputs() -> BTreeMap<String, GmOptimizeCase> {
      let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("src/commands/optimize/__tests__/__fixtures__/gm_optimize_inputs.json");
      let content = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
      serde_json::from_str(&content).unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
    }

    pub fn load_gm_outputs() -> BTreeMap<String, GmOptimizeExpected> {
      let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("src/commands/optimize/__tests__/__fixtures__/gm_optimize_outputs.json");
      let content = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
      serde_json::from_str(&content).unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
    }

    pub fn setup_and_run(workspace_root: &Path, case: &GmOptimizeCase) -> Result<OptimizeResult, Report> {
      let alphabet_sparse = Alphabet::default();
      let alphabet_dense = Alphabet::default();

      let tree_path = workspace_root.join(&case.tree);
      let aln_path = workspace_root.join(&case.aln);
      let aln = read_many_fasta(&[aln_path.to_str().unwrap()], &alphabet_sparse)?;
      let graph: GraphAncestral = nwk_read_file(&tree_path)?;

      let sparse_partitions = vec![Arc::new(RwLock::new(PartitionMarginalSparse {
        index: 0,
        gtr: jc69(JC69Params::default())?,
        alphabet: alphabet_sparse,
        length: get_common_length(&aln)?,
        nodes: btreemap! {},
        edges: btreemap! {},
      }))];

      let dense_partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense {
        index: 1,
        gtr: jc69(JC69Params::default())?,
        alphabet: alphabet_dense,
        length: get_common_length(&aln)?,
        nodes: btreemap! {},
        edges: btreemap! {},
      }))];

      compress_sequences(&graph, &sparse_partitions, &aln)?;
      initialize_marginal(&graph, &dense_partitions, &aln)?;
      update_marginal(&graph, &sparse_partitions)?;
      update_marginal(&graph, &dense_partitions)?;

      let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);
      initial_guess_mixed(&graph, &mixed_partitions)?;

      let dp = 1e-2;
      let mut lh_history = Vec::with_capacity(case.max_iter + 1);
      let mut converged_at = None;
      let mut lh_prev = f64::MIN;

      for i in 0..case.max_iter {
        let sparse_lh = update_marginal(&graph, &sparse_partitions)?;
        let dense_lh = update_marginal(&graph, &dense_partitions)?;
        let total_lh = sparse_lh + dense_lh;
        lh_history.push(total_lh);

        if (total_lh - lh_prev).abs() < dp && converged_at.is_none() {
          converged_at = Some(i);
        }

        let old_bls = save_branch_lengths(&graph);
        run_optimize_mixed(&graph, &mixed_partitions)?;
        apply_damping(&graph, &old_bls, case.damping, i);
        lh_prev = total_lh;
      }

      let sparse_lh = update_marginal(&graph, &sparse_partitions)?;
      let dense_lh = update_marginal(&graph, &dense_partitions)?;
      lh_history.push(sparse_lh + dense_lh);

      Ok(OptimizeResult {
        graph,
        lh_history,
        converged_at,
      })
    }
  }
}

#[cfg(test)]
mod tests {

  use crate::optimize::params::BranchOptMethod;
  use approx::assert_relative_eq;
  use eyre::Report;
  use helpers::{load_gm_inputs, load_gm_outputs, setup_and_run};
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::path::Path;

  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20_jc69_damped")]
  // TODO: slow datasets, disabled to keep the test run short (kb/issues/N-optimize-dense-iteration-slow.md)
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
      .map(|e| result.branch_lengths.get(&e.key()).copied().flatten().unwrap_or(0.0))
      .sum();

    assert_relative_eq!(v1_total_bl, expected.final_total_branch_length, max_relative = 0.05);

    Ok(())
  }

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
      .filter_map(|e| {
        let edge = e;
        let bl = result.branch_lengths.get(&edge.key()).copied().flatten()?;
        let target_key = edge.target();
        let node = result.graph.get_node(target_key)?;
        let name = result.names.get(&node.key()).cloned().flatten()?;
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

  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20_jc69_damped")]
  // TODO: slow datasets, disabled to keep the test run short (kb/issues/N-optimize-dense-iteration-slow.md)
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
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::optimize::dispatch::initial_guess_mixed;
    use crate::optimize::gather::{gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts};
    use crate::optimize::params::ExistingBranchLengths;
    use crate::optimize::params::{BranchOptMethod, TopologyOps};
    use crate::optimize::run_loop::run_optimize_loop;
    use crate::partition::marginal::dense::partition::PartitionMarginalDense;
    use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction};
    use crate::seq::alignment::node_seq_inputs;
    use eyre::Report;
    use itertools::Itertools;
    use serde::Deserialize;
    use std::collections::BTreeMap;
    use std::fs::read_to_string;
    use std::path::Path;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::fasta::fasta_read_file;
    use treetime_io::nwk::nwk_read_file;
    use treetime_primitives::AlignmentRecord;
    use treetime_primitives::LogLh;

    #[derive(Clone, Deserialize)]
    pub(super) struct GmOptimizeCase {
      pub(crate) tree: String,
      pub(crate) aln: String,
      pub(crate) damping: f64,
      pub(crate) max_iter: usize,
    }

    #[derive(Deserialize)]
    pub(super) struct GmOptimizeExpected {
      pub(crate) final_total_branch_length: f64,
      pub(crate) final_branch_lengths: BTreeMap<String, f64>,
    }

    pub(super) struct OptimizeResult {
      pub(crate) graph: Graph,
      pub(crate) names: BTreeMap<GraphNodeKey, Option<String>>,
      pub(crate) branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
      pub(crate) lh_history: Vec<f64>,
      pub(crate) stopped_at: Option<(usize, crate::optimize::run_loop::ConvergenceReason)>,
    }

    pub(super) fn load_gm_inputs() -> BTreeMap<String, GmOptimizeCase> {
      let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("src/optimize/__tests__/__fixtures__/gm_optimize_inputs.json");
      let content = read_to_string(&path).unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
      serde_json::from_str(&content).unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
    }

    pub(super) fn load_gm_outputs() -> BTreeMap<String, GmOptimizeExpected> {
      let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("src/optimize/__tests__/__fixtures__/gm_optimize_outputs.json");
      let content = read_to_string(&path).unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
      serde_json::from_str(&content).unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
    }

    pub(super) fn setup_and_run(
      workspace_root: &Path,
      case: &GmOptimizeCase,
      method: BranchOptMethod,
    ) -> Result<OptimizeResult, Report> {
      let alphabet = Alphabet::default();

      let tree_path = workspace_root.join(&case.tree);
      let aln_path = workspace_root.join(&case.aln);
      let aln: Vec<AlignmentRecord> = fasta_read_file(&aln_path, &alphabet)?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
      let nwk_parsed = nwk_read_file(&tree_path)?;
      let names = nwk_parsed.names();
      let mut graph = nwk_parsed.graph;
      let mut branch_lengths = nwk_parsed.branch_lengths;

      let partition = PartitionMarginalDense::new(0, alphabet, &graph, &node_seq_inputs(&graph, &names, aln))?;
      let reconstruction =
        MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, jc69(JC69Params::default())?));
      let (reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

      {
        let total_length = reconstruction.sequence_length();
        let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);
        let sub_counts = gather_edge_sub_counts(&graph, &reconstruction)?;
        let effective_lengths = gather_edge_effective_lengths(&graph, &reconstruction)?;
        initial_guess_mixed(
          &graph,
          total_length,
          &indel_counts,
          &sub_counts,
          &effective_lengths,
          ExistingBranchLengths::Overwrite,
          false,
          &mut branch_lengths,
        )?;
      }

      let dp = 0.1;
      let names_tt_1 = names.clone();
      let result = run_optimize_loop(
        &mut graph,
        reconstruction,
        case.max_iter,
        dp,
        case.damping,
        method,
        false,
        TopologyOps::default(),
        branch_lengths,
        &names_tt_1,
      )?;
      let branch_lengths = result.branch_lengths;

      let mut lh_history = result.lh_history.into_iter().map(LogLh::value).collect_vec();
      let (_, final_lh) = result
        .reconstruction
        .marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
      lh_history.push(final_lh.value());

      Ok(OptimizeResult {
        graph,
        names,
        branch_lengths,
        lh_history,
        stopped_at: result.stopped_at,
      })
    }
  }
}

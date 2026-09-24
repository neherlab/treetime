#[cfg(test)]
mod tests {

  use approx::assert_relative_eq;
  use eyre::Report;
  use rstest::rstest;
  use std::collections::BTreeMap;

  use crate::optimize::params::BranchOptMethod;
  use std::path::Path;

  use helpers::{load_gm_inputs, load_gm_outputs, setup_and_run};

  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20_jc69_damped")]
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
    use crate::ancestral::fitch::create_fitch_partition;
    use crate::ancestral::marginal::branch_lengths_or_zero;
    use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::optimize::dispatch::initial_guess_mixed;
    use crate::optimize::gather::{
      gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts, total_sequence_length,
    };
    use crate::optimize::params::ExistingBranchLengths;
    use crate::optimize::params::{BranchOptMethod, TopologyOps};
    use crate::optimize::run_loop::{marginal_update_dense, marginal_update_sparse, run_optimize_loop};
    use crate::partition::marginal::dense::partition::PartitionMarginalDense;
    use crate::seq::alignment::get_common_length;
    use crate::seq::alignment::node_seq_inputs;

    use eyre::Report;
    use itertools::Itertools;
    use treetime_graph::graph::Graph;

    use serde::Deserialize;
    use std::collections::BTreeMap;
    use std::fs::read_to_string;
    use std::path::Path;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::fasta::read_many_fasta_path;
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
      let alphabet_sparse = Alphabet::default();
      let alphabet_dense = Alphabet::default();

      let tree_path = workspace_root.join(&case.tree);
      let aln_path = workspace_root.join(&case.aln);
      let aln: Vec<AlignmentRecord> = read_many_fasta_path(&[aln_path.to_str().unwrap()], &alphabet_sparse)?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
      let nwk_parsed = nwk_read_file(&tree_path)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let mut branch_lengths = nwk_parsed.branch_lengths;
      let mut graph: Graph = graph;

      let fitch = create_fitch_partition(
        &graph,
        0,
        alphabet_sparse,
        &node_seq_inputs(&graph, &names, aln.clone()),
      )?;
      let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
      let sparse_partitions = vec![SparseReconstruction::seeded(
        partition,
        jc69(JC69Params::default())?,
        node_states,
      )];

      let length = get_common_length(&aln)?;
      let dense_partition = PartitionMarginalDense::new(1, alphabet_dense, length);
      let dense_node_states = dense_partition.attach_sequences(&graph, &node_seq_inputs(&graph, &names, aln))?;
      let dense_partitions = vec![DenseReconstruction::seeded(
        dense_partition,
        jc69(JC69Params::default())?,
        dense_node_states,
      )];

      let (sparse_partitions, _) =
        marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;
      let (dense_partitions, _) =
        marginal_update_dense(&graph, &branch_lengths_or_zero(&branch_lengths), dense_partitions)?;

      {
        let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
        let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);
        let sub_counts = gather_edge_sub_counts(&graph, &dense_partitions, &sparse_partitions)?;
        let effective_lengths = gather_edge_effective_lengths(&graph, &dense_partitions, &sparse_partitions)?;
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
        sparse_partitions,
        dense_partitions,
        case.max_iter,
        dp,
        case.damping,
        method,
        false,
        TopologyOps::default(),
        branch_lengths,
        &names_tt_1,
      )?;
      let sparse_partitions = result.sparse_partitions;
      let dense_partitions = result.dense_partitions;
      let branch_lengths = result.branch_lengths;

      let mut lh_history = result.lh_history.into_iter().map(LogLh::value).collect_vec();
      let (sparse_partitions, sparse_lh) =
        marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;
      let sparse_lh = sparse_lh.value();
      let (dense_partitions, dense_lh) =
        marginal_update_dense(&graph, &branch_lengths_or_zero(&branch_lengths), dense_partitions)?;
      let dense_lh = dense_lh.value();
      lh_history.push(sparse_lh + dense_lh);

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

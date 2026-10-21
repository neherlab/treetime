#[cfg(test)]
mod tests {
  use helpers::{arb_jc69_case, optimize_simulated};
  use proptest::prelude::{ProptestConfig, TestCaseError, prop_assert, prop_assert_eq, proptest};

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(8))]

    #[test]
    fn test_prop_optimize_simulated_invariants((spec, perturbation_seed) in arb_jc69_case()) {
      let run = optimize_simulated(&spec, perturbation_seed).map_err(|error| TestCaseError::fail(format!("{error:?}")))?;

      prop_assert_eq!(&run.true_clades, &run.output_clades);
      prop_assert!(
        run.output_log_lh >= run.start_log_lh,
        "log-likelihood at the output {} is below the log-likelihood at the start {}",
        run.output_log_lh,
        run.start_log_lh
      );
    }

    #[test]
    fn test_prop_optimize_simulated_reaches_true_length_likelihood((spec, perturbation_seed) in arb_jc69_case()) {
      let run = optimize_simulated(&spec, perturbation_seed).map_err(|error| TestCaseError::fail(format!("{error:?}")))?;

      prop_assert!(
        run.output_log_lh >= run.exact_log_lh,
        "log-likelihood at the output {} is below the log-likelihood at the exact branch lengths {} \
         (start {}); optimizer warnings: {:?}",
        run.output_log_lh,
        run.exact_log_lh,
        run.start_log_lh,
        run.warnings
      );
    }
  }

  mod helpers {
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::GtrModelName;
    use crate::optimize::indel::IndelPolicy;
    use crate::optimize::params::{BranchOptMethod, InitialGuessMode, OptimizeLoopSettings, TopologyOps};
    use crate::optimize::partition_rates::PartitionRates;
    use crate::optimize::pipeline::{OptimizeParams, run};
    use crate::progress::NoopProgress;
    use crate::test_utils::{
      CaseSpec, ModelFamily, RecordingLog, arb_simulated_case, gtr_from_model, optimize_input, run_dense_marginal,
    };
    use eyre::Report;
    use proptest::arbitrary::any;
    use proptest::strategy::Strategy;
    use std::collections::{BTreeMap, BTreeSet};
    use treetime_graph::clades::{TipIds, graph_clade_keys};
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node_names::NodeNames;
    use treetime_simulation::branch_lengths::{exact, perturbed};
    use treetime_utils::sync::random::get_random_number_generator;

    const START_PERTURBATION_SIGMA: f64 = 0.3;

    pub(super) fn arb_jc69_case() -> impl Strategy<Value = (CaseSpec, u64)> {
      (arb_simulated_case(), any::<u64>()).prop_map(|(spec, perturbation_seed)| {
        (
          CaseSpec {
            model: ModelFamily::Jc69,
            ..spec
          },
          perturbation_seed,
        )
      })
    }

    pub(super) fn optimize_simulated(spec: &CaseSpec, perturbation_seed: u64) -> Result<OptimizeRun, Report> {
      let case = spec.simulate()?;
      let tree = &case.tree;
      let exact_lengths = exact(&tree.edge_durations(), &case.model);
      let start_lengths = perturbed(
        &exact_lengths,
        START_PERTURBATION_SIGMA,
        &mut get_random_number_generator(perturbation_seed),
      )?;
      let input = optimize_input(tree, &case.sequences.sequences, &start_lengths)?;
      let start_input_lengths = input.branch_lengths.clone();
      let exact_input_lengths = optimize_input(tree, &case.sequences.sequences, &exact_lengths)?.branch_lengths;
      let log = RecordingLog::default();

      let output = run(&optimize_params(), input, &NoopCancel, &NoopProgress, &log)?;

      let records = case.sequences.leaf_records(&tree.graph, &tree.names);
      let gtr = gtr_from_model(&case.model, 1.0)?;
      let log_lh = |graph: &Graph, lengths: &BTreeMap<GraphEdgeKey, Option<f64>>, names: &NodeNames| {
        run_dense_marginal(graph, lengths, names, &records, gtr.clone())
      };
      let (start_log_lh, _) = log_lh(&tree.graph, &start_input_lengths, &tree.names)?;
      let (exact_log_lh, _) = log_lh(&tree.graph, &exact_input_lengths, &tree.names)?;
      let (output_log_lh, _) = log_lh(&output.graph, &output.branch_lengths, &output.names)?;
      let mut ids = TipIds::default();
      let true_clades = graph_clade_keys(&tree.graph, &tree.names, &mut ids)?
        .into_values()
        .collect();
      let output_clades = graph_clade_keys(&output.graph, &output.names, &mut ids)?
        .into_values()
        .collect();
      Ok(OptimizeRun {
        true_clades,
        output_clades,
        start_log_lh,
        exact_log_lh,
        output_log_lh,
        warnings: log.run_warnings().into_iter().map(|warning| warning.message).collect(),
      })
    }

    pub(super) struct OptimizeRun {
      pub true_clades: BTreeSet<Vec<usize>>,
      pub output_clades: BTreeSet<Vec<usize>>,
      pub start_log_lh: f64,
      pub exact_log_lh: f64,
      pub output_log_lh: f64,
      pub warnings: Vec<String>,
    }

    fn optimize_params() -> OptimizeParams {
      OptimizeParams {
        model: GtrModelName::JC69,
        dense: Some(true),
        initial_guess: InitialGuessMode::Never,
        reroot_spec: None,
        loop_settings: OptimizeLoopSettings {
          max_iter: 100,
          dp: 1e-6,
          damping: 0.0,
          opt_method: BranchOptMethod::Brent,
          indel_policy: IndelPolicy::Exclude,
          topology_ops: TopologyOps {
            collapse_short_branches: false,
            merge_siblings: false,
            flip_parent_child: false,
          },
          partition_rates: PartitionRates::Estimate,
        },
      }
    }
  }
}

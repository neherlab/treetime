#[cfg(test)]
mod tests {
  use crate::test_utils::{decrease_check, recovery_check, require_passed, run_replicates};
  use eyre::Report;
  use helpers::{
    MARGIN, SITE_COMPARISON_BASE_SEED, SITE_COMPARISON_REPLICATES, TOTAL_LENGTH_BASE_SEED, TOTAL_LENGTH_REPLICATES,
    branch_length_rmse_log_ratio, total_length_log_ratio,
  };
  use treetime_simulation::FALSE_FAIL_RATE;

  #[test]
  fn test_sim_optimize_total_branch_length_recovers_the_exact_total() -> Result<(), Report> {
    let errors = run_replicates(TOTAL_LENGTH_BASE_SEED, TOTAL_LENGTH_REPLICATES, total_length_log_ratio)?;

    require_passed(&[recovery_check(
      "log ratio of the optimized and the exact total branch length",
      &errors,
      MARGIN,
      FALSE_FAIL_RATE,
    )?])
  }

  #[test]
  fn test_sim_optimize_branch_length_error_falls_with_more_sites() -> Result<(), Report> {
    let log_ratios = run_replicates(
      SITE_COMPARISON_BASE_SEED,
      SITE_COMPARISON_REPLICATES,
      branch_length_rmse_log_ratio,
    )?;

    require_passed(&[decrease_check(
      "log ratio of the branch length RMSE with all and with the first sites",
      &log_ratios,
      FALSE_FAIL_RATE,
    )?])
  }

  mod helpers {
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::{GtrModel, JC69Params};
    use crate::optimize::indel::IndelPolicy;
    use crate::optimize::params::{BranchOptMethod, InitialGuessMode, OptimizeLoopSettings, TopologyOps};
    use crate::optimize::partition_rates::PartitionRates;
    use crate::optimize::pipeline::{OptimizeOutput, OptimizeParams, run};
    use crate::progress::NoopProgress;
    use crate::test_utils::{CaseSpec, ModelFamily, SimulatedCase, first_columns, optimize_input};
    use eyre::Report;
    use itertools::Itertools;
    use rand::RngExt;
    use std::collections::BTreeMap;
    use treetime_graph::clades::{TipIds, graph_clade_keys};
    use treetime_graph::edge::Edge;
    use treetime_graph::key_map::EdgeMap;
    use treetime_simulation::branch_lengths::{exact, perturbed};
    use treetime_simulation::coalescent::SampleBatch;
    use treetime_simulation::tc_schedule::TcSchedule;
    use treetime_utils::cast::ToF64Exact;
    use treetime_utils::make_report;
    use treetime_utils::sync::random::get_random_number_generator;

    pub(super) const MARGIN: f64 = 0.03;
    pub(super) const TOTAL_LENGTH_BASE_SEED: u64 = 1;
    pub(super) const TOTAL_LENGTH_REPLICATES: usize = 1430;
    pub(super) const SITE_COMPARISON_BASE_SEED: u64 = 1;
    pub(super) const SITE_COMPARISON_REPLICATES: usize = 9;
    const BATCH_COUNT: usize = 8;
    const TIPS_PER_BATCH: usize = 5;
    const SAMPLING_WINDOW: f64 = 2.0;
    const CONSTANT_TC: f64 = 5.0;
    const MU: f64 = 0.01;
    const TOTAL_LENGTH_SITES: usize = 1000;
    const FIRST_SITES: usize = 500;
    const ALL_SITES: usize = 4000;
    const START_PERTURBATION_SIGMA: f64 = 0.3;

    pub(super) fn total_length_log_ratio(seed: u64) -> Result<f64, Report> {
      let replicate = Replicate::simulate(seed, TOTAL_LENGTH_SITES)?;

      let output = replicate.optimize(TOTAL_LENGTH_SITES)?;

      let optimized_total = output
        .branch_lengths
        .iter()
        .map(|(edge, length)| length.ok_or_else(|| make_report!("Edge {edge} has no optimized branch length")))
        .sum::<Result<f64, Report>>()?;
      let exact_total: f64 = replicate.exact.values().sum();
      Ok((optimized_total / exact_total).ln())
    }

    pub(super) fn branch_length_rmse_log_ratio(seed: u64) -> Result<f64, Report> {
      let replicate = Replicate::simulate(seed, ALL_SITES)?;

      let first_sites_rmse = replicate.branch_length_rmse(&replicate.optimize(FIRST_SITES)?)?;
      let all_sites_rmse = replicate.branch_length_rmse(&replicate.optimize(ALL_SITES)?)?;

      Ok((all_sites_rmse / first_sites_rmse).ln())
    }

    struct Replicate {
      case: SimulatedCase,
      exact: EdgeMap<f64>,
      start: EdgeMap<f64>,
    }

    impl Replicate {
      fn simulate(seed: u64, length: usize) -> Result<Self, Report> {
        let mut rng = get_random_number_generator(seed);
        let spec = CaseSpec {
          seed: rng.random(),
          sampling: forty_tip_sampling(),
          schedule: TcSchedule::constant(CONSTANT_TC)?,
          model: ModelFamily::Jc69,
          mu: MU,
          length,
        };
        let case = spec.simulate()?;
        let exact = exact(&case.tree.edge_durations(), &case.model)?;
        let start = perturbed(&exact, START_PERTURBATION_SIGMA, &mut rng)?;
        Ok(Self { case, exact, start })
      }

      fn optimize(&self, columns: usize) -> Result<OptimizeOutput, Report> {
        let sequences = first_columns(&self.case.sequences.sequences, columns)?;
        let input = optimize_input(
          &self.case.tree,
          &sequences,
          &self.start,
          &GtrModel::JC69(JC69Params::default()),
        )?;
        Ok(run(
          &optimize_params(),
          input,
          &NoopCancel,
          &NoopProgress,
          &NoopProgress,
        )?)
      }

      fn branch_length_rmse(&self, output: &OptimizeOutput) -> Result<f64, Report> {
        let tree = &self.case.tree;
        let mut ids = TipIds::default();
        let true_clades = graph_clade_keys(&tree.graph, &tree.names, &mut ids)?;
        let output_clades = graph_clade_keys(&output.graph, &output.names, &mut ids)?;
        let optimized = output
          .graph
          .get_edges()
          .map(|edge| {
            let length = output.branch_lengths[&edge.key()]
              .ok_or_else(|| make_report!("Edge {} has no optimized branch length", edge.key()))?;
            Ok((output_clades[&edge.target()].as_slice(), length))
          })
          .collect::<Result<BTreeMap<_, _>, Report>>()?;
        let error = |edges: &[&Edge]| {
          edges
            .iter()
            .map(|edge| {
              let length = optimized
                .get(true_clades[&edge.target()].as_slice())
                .ok_or_else(|| make_report!("The output has no edge above node {}", &tree.names[edge.target()]))?;
              Ok(length - self.exact[&edge.key()])
            })
            .sum::<Result<f64, Report>>()
        };
        let root = tree.graph.get_exactly_one_root()?.key();
        let (root_edges, other_edges): (Vec<_>, Vec<_>) =
          tree.graph.get_edges().partition(|edge| edge.source() == root);
        let errors = other_edges
          .iter()
          .map(|&edge| error(&[edge]))
          .chain([error(&root_edges)])
          .try_collect::<_, Vec<_>, _>()?;
        let squares: f64 = errors.iter().map(|error| error * error).sum();
        Ok((squares / errors.len().to_f64_exact()).sqrt())
      }
    }

    fn forty_tip_sampling() -> Vec<SampleBatch> {
      (0..BATCH_COUNT)
        .map(|batch| SampleBatch {
          tau: SAMPLING_WINDOW * batch.to_f64_exact() / (BATCH_COUNT - 1).to_f64_exact(),
          count: TIPS_PER_BATCH,
        })
        .collect()
    }

    fn optimize_params() -> OptimizeParams {
      OptimizeParams {
        dense: None,
        initial_guess: InitialGuessMode::Never,
        reroot_spec: None,
        loop_settings: OptimizeLoopSettings {
          max_iter: 10,
          dp: 0.1,
          damping: 0.0,
          opt_method: BranchOptMethod::Brent,
          indel_policy: IndelPolicy::Include,
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

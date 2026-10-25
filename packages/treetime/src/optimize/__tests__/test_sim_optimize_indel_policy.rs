#[cfg(test)]
mod tests {
  use eyre::Report;
  use helpers::{Representation, check_poisson_relation};

  #[test]
  fn test_sim_optimize_indel_policy_dense_include_follows_the_poisson_relation_to_exclude() -> Result<(), Report> {
    check_poisson_relation(Representation::Dense)
  }

  #[test]
  fn test_sim_optimize_indel_policy_sparse_include_follows_the_poisson_relation_to_exclude() -> Result<(), Report> {
    check_poisson_relation(Representation::Sparse)
  }

  mod helpers {
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::{GtrModel, JC69Params};
    use crate::optimize::indel::IndelPolicy;
    use crate::optimize::params::{BranchOptMethod, InitialGuessMode, OptimizeLoopSettings, TopologyOps};
    use crate::optimize::partition_rates::PartitionRates;
    use crate::optimize::pipeline::{OptimizeOutput, OptimizeParams, run};
    use crate::progress::NoopProgress;
    use crate::test_utils::{CaseSpec, ModelFamily, optimize_input, recovery_check, require_passed, run_replicates};
    use eyre::Report;
    use itertools::Itertools;
    use ndarray::Array1;
    use rand::RngExt;
    use std::collections::BTreeMap;
    use treetime_graph::edge::{Edge, GraphEdgeKey};
    use treetime_graph::graph::Graph;
    use treetime_simulation::FALSE_FAIL_RATE;
    use treetime_simulation::branch_lengths::{exact, perturbed};
    use treetime_simulation::coalescent::SampleBatch;
    use treetime_simulation::gaps::{GAP, apply_gap_events, disjoint_gap_events};
    use treetime_simulation::tc_schedule::TcSchedule;
    use treetime_utils::cast::ToF64Exact;
    use treetime_utils::collections::container::get_exactly_one;
    use treetime_utils::sync::random::get_random_number_generator;
    use treetime_utils::{make_error, make_report, pretty_assert_ulps_eq};

    const MARGIN: f64 = 0.01;
    const BASE_SEED: u64 = 1;
    const REPLICATES: usize = 8;
    const BATCH_COUNT: usize = 8;
    const TIPS_PER_BATCH: usize = 5;
    const SAMPLING_WINDOW: f64 = 2.0;
    const CONSTANT_TC: f64 = 5.0;
    const MU: f64 = 3e-3;
    const SITES: usize = 1000;
    const GAP_EVENTS: usize = 10;
    const GAP_LENGTH: usize = 3;
    const START_PERTURBATION_SIGMA: f64 = 0.3;
    const MIN_SUBSTITUTIONS: f64 = 0.5;

    pub(super) fn check_poisson_relation(representation: Representation) -> Result<(), Report> {
      let replicates = run_replicates(BASE_SEED, REPLICATES, |seed| paired_runs(seed, representation))?;

      let expected_rates: Array1<f64> = replicates.iter().map(|pair| pair.expected_indel_rate).collect();
      let rates: Array1<f64> = replicates.iter().map(|pair| pair.indel_rate).collect();
      pretty_assert_ulps_eq!(expected_rates, rates, max_ulps = 4);
      let ratios = replicates
        .iter()
        .flat_map(|pair| pair.ratios.iter().copied())
        .collect_vec();
      require_passed(&[
        recovery_check(
          "log ratio of the Include length and its Poisson prediction (L t + k) / (L + rate), branches with gap events",
          &ratios.iter().filter_map(|ratio| ratio.with_gap()).collect_vec(),
          MARGIN,
          FALSE_FAIL_RATE / 2.0,
        )?,
        recovery_check(
          "log ratio of the Include length and its Poisson prediction L t / (L + rate), branches without gap events",
          &ratios.iter().filter_map(|ratio| ratio.without_gap()).collect_vec(),
          MARGIN,
          FALSE_FAIL_RATE / 2.0,
        )?,
      ])
    }

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Representation {
      Dense,
      Sparse,
    }

    struct PairedRuns {
      expected_indel_rate: f64,
      indel_rate: f64,
      ratios: Vec<EdgeRatio>,
    }

    #[derive(Clone, Copy, Debug)]
    enum EdgeRatio {
      WithGap(f64),
      WithoutGap(f64),
    }

    impl EdgeRatio {
      fn with_gap(self) -> Option<f64> {
        match self {
          Self::WithGap(ratio) => Some(ratio),
          Self::WithoutGap(_) => None,
        }
      }

      fn without_gap(self) -> Option<f64> {
        match self {
          Self::WithGap(_) => None,
          Self::WithoutGap(ratio) => Some(ratio),
        }
      }
    }

    fn paired_runs(seed: u64, representation: Representation) -> Result<PairedRuns, Report> {
      let mut rng = get_random_number_generator(seed);
      let case = CaseSpec {
        seed: rng.random(),
        sampling: forty_tip_sampling(),
        schedule: TcSchedule::constant(CONSTANT_TC)?,
        model: ModelFamily::Jc69,
        mu: MU,
        length: SITES,
      }
      .simulate()?;
      let tree = &case.tree;
      let durations = tree.edge_durations();
      let events = disjoint_gap_events(&durations, GAP_EVENTS, GAP_LENGTH, SITES, &mut rng)?;
      let sequences = apply_gap_events(&tree.graph, &case.sequences.sequences, &events)?;
      let start = perturbed(&exact(&durations, &case.model)?, START_PERTURBATION_SIGMA, &mut rng)?;
      let optimize = |indel_policy| -> Result<OptimizeOutput, Report> {
        let input = optimize_input(tree, &sequences, &start, &GtrModel::JC69(JC69Params::default()))?;
        let output = run(
          &optimize_params(representation, indel_policy),
          input,
          &NoopCancel,
          &NoopProgress,
          &NoopProgress,
        )?;
        require_same_edges(&tree.graph, &output.graph)?;
        Ok(output)
      };
      let exclude = optimize(IndelPolicy::Exclude)?;
      let include = optimize(IndelPolicy::Include)?;

      let indel_rate = get_exactly_one(&include.alignments)?
        .indel_rate
        .ok_or_else(|| make_report!("The Include run reports no indel rate"))?;
      let events_per_edge = events.iter().counts_by(|event| event.edge);
      let root = tree.graph.get_exactly_one_root()?.key();
      let edge_ratio = |edge: &Edge| -> Result<Option<EdgeRatio>, Report> {
        let gapped_sites = sequences[&edge.target()].iter().filter(|&&state| state == GAP).count();
        let effective_length = (SITES - gapped_sites).to_f64_exact();
        let exclude_substitutions = effective_length * length_of(&exclude, edge.key())?;
        let include_substitutions = (effective_length + indel_rate) * length_of(&include, edge.key())?;
        Ok(match events_per_edge.get(&edge.key()) {
          Some(&count) => Some(EdgeRatio::WithGap(
            (include_substitutions / (exclude_substitutions + count.to_f64_exact())).ln(),
          )),
          None if exclude_substitutions >= MIN_SUBSTITUTIONS => Some(EdgeRatio::WithoutGap(
            (include_substitutions / exclude_substitutions).ln(),
          )),
          None => None,
        })
      };
      let ratios = tree
        .graph
        .get_edges()
        .filter(|edge| edge.source() != root)
        .map(edge_ratio)
        .filter_map_ok(|ratio| ratio)
        .try_collect()?;
      Ok(PairedRuns {
        expected_indel_rate: events.len().to_f64_exact() / start.values().sum::<f64>(),
        indel_rate,
        ratios,
      })
    }

    fn length_of(output: &OptimizeOutput, edge: GraphEdgeKey) -> Result<f64, Report> {
      output
        .branch_lengths
        .get(edge)
        .copied()
        .flatten()
        .ok_or_else(|| make_report!("Edge {edge} has no optimized branch length"))
    }

    fn require_same_edges(input: &Graph, output: &Graph) -> Result<(), Report> {
      let endpoints = |graph: &Graph| -> BTreeMap<_, _> {
        graph
          .get_edges()
          .map(|edge| (edge.key(), (edge.source(), edge.target())))
          .collect()
      };
      if endpoints(input) != endpoints(output) {
        return make_error!("The optimize run changed the topology, but the comparison needs the input edges");
      }
      Ok(())
    }

    fn forty_tip_sampling() -> Vec<SampleBatch> {
      (0..BATCH_COUNT)
        .map(|batch| SampleBatch {
          tau: SAMPLING_WINDOW * batch.to_f64_exact() / (BATCH_COUNT - 1).to_f64_exact(),
          count: TIPS_PER_BATCH,
        })
        .collect()
    }

    fn optimize_params(representation: Representation, indel_policy: IndelPolicy) -> OptimizeParams {
      OptimizeParams {
        dense: Some(matches!(representation, Representation::Dense)),
        initial_guess: InitialGuessMode::Never,
        reroot_spec: None,
        loop_settings: OptimizeLoopSettings {
          max_iter: 30,
          dp: 1e-6,
          damping: 0.0,
          opt_method: BranchOptMethod::Brent,
          indel_policy,
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

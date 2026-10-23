#[cfg(test)]
pub(crate) mod tests {
  use eyre::Report;
  use helpers::{Calibration, require_mean_not_rejected};
  use statrs::stats_tests::Alternative;

  #[test]
  fn test_sim_ancestral_calibration_dense_posterior_predicts_the_error_of_the_most_likely_state() -> Result<(), Report>
  {
    let calibration = Calibration::simulate()?;

    let errors = calibration.site_errors(&calibration.dense_posteriors(calibration.true_gtr()?)?)?;

    calibration.require_error_rate_in_range(&errors)?;
    require_mean_not_rejected(
      "dense: wrong most likely states minus their posterior error probability, per site",
      &errors.deviations(),
      Alternative::TwoSided,
    )
  }

  #[test]
  fn test_sim_ancestral_calibration_sparse_posterior_predicts_the_error_of_the_most_likely_state() -> Result<(), Report>
  {
    let calibration = Calibration::simulate()?;

    let errors = calibration.site_errors(&calibration.sparse_posteriors(calibration.true_gtr()?)?)?;

    require_mean_not_rejected(
      "sparse: wrong most likely states minus their posterior error probability, per site",
      &errors.deviations(),
      Alternative::TwoSided,
    )
  }

  #[test]
  fn test_sim_ancestral_calibration_dense_posterior_with_coarsened_tips_predicts_the_error_of_the_most_likely_state()
  -> Result<(), Report> {
    let calibration = Calibration::simulate_coarsened()?;

    let errors = calibration.site_errors(&calibration.dense_posteriors(calibration.true_gtr()?)?)?;

    calibration.require_error_rate_in_range(&errors)?;
    require_mean_not_rejected(
      "dense, coarsened tips: wrong most likely states minus their posterior error probability, per site",
      &errors.deviations(),
      Alternative::TwoSided,
    )
  }

  #[test]
  fn test_sim_ancestral_calibration_sparse_posterior_with_coarsened_tips_predicts_the_error_of_the_most_likely_state()
  -> Result<(), Report> {
    let calibration = Calibration::simulate_coarsened()?;

    let errors = calibration.site_errors(&calibration.sparse_posteriors(calibration.true_gtr()?)?)?;

    require_mean_not_rejected(
      "sparse, coarsened tips: wrong most likely states minus their posterior error probability, per site",
      &errors.deviations(),
      Alternative::TwoSided,
    )
  }

  #[test]
  fn test_sim_ancestral_calibration_marginal_makes_no_more_errors_than_parsimony() -> Result<(), Report> {
    let calibration = Calibration::simulate()?;

    let marginal = calibration.site_errors(&calibration.dense_posteriors(calibration.true_gtr()?)?)?;
    let parsimony = calibration.parsimony_wrong_states()?;

    require_mean_not_rejected(
      "parsimony minus marginal wrong internal states, per site",
      &(parsimony - &marginal.wrong),
      Alternative::Less,
    )
  }

  pub(crate) mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::ancestral::params::{AncestralParams, MethodAncestral};
    use crate::ancestral::pipeline;
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::GtrModel;
    use crate::gtr::gtr::GTR;
    use crate::partition::marginal::posterior_states::most_likely_state;
    use crate::partition::marginal::sample::SampleMode;
    use crate::progress::NoopProgress;
    use crate::seq::sink::SeqStream;
    use crate::test_utils::{
      CaseSpec, ModelFamily, RecordingSeqSink, SimulatedCase, ancestral_input, gtr_from_model, run_dense_marginal,
      run_sparse_marginal,
    };
    use eyre::Report;
    use itertools::Itertools;
    use ndarray::{Array1, Array2, aview1};
    use statrs::stats_tests::ttest_onesample::ttest_onesample;
    use statrs::stats_tests::{Alternative, NaNPolicy};
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::node::GraphNodeKey;
    use treetime_primitives::{AlignmentRecord, Seq};
    use treetime_simulation::FALSE_FAIL_RATE;
    use treetime_simulation::branch_lengths::exact;
    use treetime_simulation::coalescent::SampleBatch;
    use treetime_simulation::evolve::leaf_records;
    use treetime_simulation::observe::{TipObservation, observe_tips};
    use treetime_simulation::tc_schedule::TcSchedule;
    use treetime_utils::cast::ToF64Exact;
    use treetime_utils::sync::random::get_random_number_generator;
    use treetime_utils::{make_error, make_report};

    const ANCESTRAL_MU: f64 = 0.4;
    const ANCESTRAL_TC: f64 = 5.0;
    const BATCH_COUNT: usize = 8;
    const TIPS_PER_BATCH: usize = 5;
    const SAMPLING_WINDOW: f64 = 2.0;
    const CALIBRATION_SEED: u64 = 1;
    const CALIBRATION_LENGTH: usize = 5000;
    const MIN_ERROR_RATE: f64 = 0.05;
    const MAX_ERROR_RATE: f64 = 0.2;
    const COARSENING_SEED: u64 = 2;

    pub(crate) fn ancestral_case(seed: u64, length: usize) -> Result<CaseSpec, Report> {
      let gaps = (BATCH_COUNT - 1).to_f64_exact();
      Ok(CaseSpec {
        seed,
        sampling: (0..BATCH_COUNT)
          .map(|batch| SampleBatch {
            tau: SAMPLING_WINDOW * batch.to_f64_exact() / gaps,
            count: TIPS_PER_BATCH,
          })
          .collect(),
        schedule: TcSchedule::constant(ANCESTRAL_TC)?,
        model: ModelFamily::RandomGtr,
        mu: ANCESTRAL_MU,
        length,
      })
    }

    pub(super) fn require_mean_not_rejected(
      name: &str,
      values: &Array1<f64>,
      alternative: Alternative,
    ) -> Result<(), Report> {
      let (_, p_value) = ttest_onesample(values.to_vec(), 0.0, alternative, NaNPolicy::Error)
        .map_err(|error| make_report!("{name}: one-sample t-test against 0: {error}"))?;
      if p_value > FALSE_FAIL_RATE {
        return Ok(());
      }
      let mean = values.mean().ok_or_else(|| make_report!("{name}: no values"))?;
      make_error!(
        "{name}: the {alternative:?} t-test rejects a mean of 0 with p = {p_value:.3e}, at or below {FALSE_FAIL_RATE:.0e}: n = {}, mean = {mean:.6}, standard deviation = {:.6}",
        values.len(),
        values.std(1.0)
      )
    }

    pub(super) struct Calibration {
      case: SimulatedCase,
      lengths: BTreeMap<GraphEdgeKey, f64>,
      records: Vec<AlignmentRecord>,
      internal_nodes: Vec<GraphNodeKey>,
    }

    impl Calibration {
      pub(super) fn simulate() -> Result<Self, Report> {
        let case = ancestral_case(CALIBRATION_SEED, CALIBRATION_LENGTH)?.simulate()?;
        let records = case.sequences.leaf_records(&case.tree.graph, &case.tree.names);
        Self::new(case, records)
      }

      pub(super) fn simulate_coarsened() -> Result<Self, Report> {
        let case = ancestral_case(CALIBRATION_SEED, CALIBRATION_LENGTH)?.simulate()?;
        let observed = observe_tips(
          &case.tree.graph,
          &case.sequences.sequences,
          &coarsening(),
          &mut get_random_number_generator(COARSENING_SEED),
        )?;
        let records = leaf_records(&case.tree.graph, &case.tree.names, &observed);
        Self::new(case, records)
      }

      fn new(case: SimulatedCase, records: Vec<AlignmentRecord>) -> Result<Self, Report> {
        let lengths = exact(&case.tree.edge_durations(), &case.model)?;
        let internal_nodes = case.tree.graph.get_internal_nodes().map(|node| node.key()).collect();
        Ok(Self {
          case,
          lengths,
          records,
          internal_nodes,
        })
      }

      pub(super) fn true_gtr(&self) -> Result<GTR, Report> {
        gtr_from_model(&self.case.model, 1.0)
      }

      pub(super) fn dense_posteriors(&self, gtr: GTR) -> Result<Posteriors, Report> {
        let (_, dense) = run_dense_marginal(
          &self.case.tree.graph,
          &self.optional_lengths(),
          &self.case.tree.names,
          &self.records,
          gtr,
        )?;
        let nodes = self
          .internal_nodes
          .iter()
          .map(|&node| (node, dense.node_states[&node].profile.dis.clone()))
          .collect();
        Ok(Posteriors {
          alphabet: dense.partition.alphabet,
          nodes,
        })
      }

      pub(super) fn sparse_posteriors(&self, gtr: GTR) -> Result<Posteriors, Report> {
        let (_, sparse) = run_sparse_marginal(
          &self.case.tree.graph,
          &self.optional_lengths(),
          &self.case.tree.names,
          &self.records,
          gtr,
        )?;
        let alphabet = sparse.partition.alphabet.clone();
        let nodes = self
          .internal_nodes
          .iter()
          .map(|&node| {
            let state = &sparse.node_states[&node];
            let mut posterior = Array2::zeros((CALIBRATION_LENGTH, alphabet.n_canonical()));
            for (pos, mut row) in posterior.rows_mut().into_iter().enumerate() {
              let dis = if let Some(var) = state.profile.variable.get(pos) {
                var.dis
              } else {
                let fixed_state = state.sequence[pos];
                state
                  .profile
                  .fixed
                  .get(fixed_state)
                  .ok_or_else(|| {
                    make_report!(
                      "Node {} has no fixed profile of state {fixed_state} at position {pos}",
                      &self.case.tree.names[node]
                    )
                  })?
                  .dis
                  .as_slice()
              };
              row.assign(&aview1(dis));
            }
            Ok((node, posterior))
          })
          .collect::<Result<_, Report>>()?;
        Ok(Posteriors { alphabet, nodes })
      }

      pub(super) fn parsimony_wrong_states(&self) -> Result<Array1<f64>, Report> {
        let params = AncestralParams {
          method: MethodAncestral::Parsimony,
          dense: None,
          include_leaves: false,
          report_ambiguous: false,
          impute_missing_data: false,
          seed: 0,
          sample_from_profile: SampleMode::Argmax,
        };
        let input = ancestral_input(
          &self.case.tree,
          &self.case.sequences.sequences,
          &self.lengths,
          &GtrModel::Infer,
        )?;
        let mut sink = RecordingSeqSink::default();
        pipeline::run(
          &params,
          input,
          Some(SeqStream::new(&mut sink, &self.case.tree.names)),
          &NoopCancel,
          &NoopProgress,
          &NoopProgress,
        )?;
        let sequences: BTreeMap<GraphNodeKey, Seq> = sink
          .items
          .into_iter()
          .filter(|(_, emitted, _)| *emitted)
          .map(|(node, _, seq)| (node, seq))
          .collect();
        self
          .internal_nodes
          .iter()
          .try_fold(Array1::zeros(CALIBRATION_LENGTH), |mut wrong: Array1<f64>, &node| {
            let reconstructed = sequences.get(&node).ok_or_else(|| {
              make_report!(
                "Parsimony reconstruction has no sequence for node {}",
                &self.case.tree.names[node]
              )
            })?;
            let truth = &self.case.sequences.sequences[&node];
            for (count, (&state, &true_state)) in wrong.iter_mut().zip(reconstructed.iter().zip_eq(truth.iter())) {
              *count += f64::from(u8::from(state != true_state));
            }
            Ok(wrong)
          })
      }

      pub(super) fn site_errors(&self, posteriors: &Posteriors) -> Result<SiteErrors, Report> {
        let mut errors = SiteErrors {
          wrong: Array1::zeros(CALIBRATION_LENGTH),
          expected: Array1::zeros(CALIBRATION_LENGTH),
        };
        for &node in &self.internal_nodes {
          let truth = &self.case.sequences.sequences[&node];
          for (pos, (row, &true_state)) in posteriors.nodes[&node]
            .rows()
            .into_iter()
            .zip_eq(truth.iter())
            .enumerate()
          {
            let most_likely = most_likely_state(&row).ok_or_else(|| {
              make_report!(
                "Posterior of node {} at position {pos} is empty",
                &self.case.tree.names[node]
              )
            })?;
            errors.wrong[pos] += f64::from(u8::from(posteriors.alphabet.char(most_likely) != true_state));
            errors.expected[pos] += 1.0 - row[most_likely];
          }
        }
        Ok(errors)
      }

      fn error_rate(&self, errors: &SiteErrors) -> f64 {
        errors.wrong.sum() / (CALIBRATION_LENGTH * self.internal_nodes.len()).to_f64_exact()
      }

      pub(super) fn require_error_rate_in_range(&self, errors: &SiteErrors) -> Result<(), Report> {
        let rate = self.error_rate(errors);
        if (MIN_ERROR_RATE..=MAX_ERROR_RATE).contains(&rate) {
          return Ok(());
        }
        make_error!(
          "The mean error of the most likely internal states is {rate:.4}, outside the range [{MIN_ERROR_RATE}, {MAX_ERROR_RATE}] that the calibration tests need"
        )
      }

      fn optional_lengths(&self) -> BTreeMap<GraphEdgeKey, Option<f64>> {
        self
          .lengths
          .iter()
          .map(|(&edge, &length)| (edge, Some(length)))
          .collect()
      }
    }

    pub(super) struct Posteriors {
      alphabet: Alphabet,
      nodes: BTreeMap<GraphNodeKey, Array2<f64>>,
    }

    pub(super) struct SiteErrors {
      pub wrong: Array1<f64>,
      pub expected: Array1<f64>,
    }

    impl SiteErrors {
      pub(super) fn deviations(&self) -> Array1<f64> {
        &self.wrong - &self.expected
      }
    }

    fn coarsening() -> TipObservation {
      TipObservation {
        ambiguity_rate: 0.1,
        unknown_run_rate: 0.25,
        unknown_run_lengths: 10..=500,
      }
    }
  }
}

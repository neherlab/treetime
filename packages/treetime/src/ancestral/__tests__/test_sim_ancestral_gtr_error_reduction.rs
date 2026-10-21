#[cfg(test)]
mod tests {
  use crate::test_utils::{decrease_check, require_passed, run_replicates};
  use eyre::Report;
  use helpers::{BASE_SEED, REPLICATE_COUNT, replicate_log_ratios};
  use treetime_simulation::FALSE_FAIL_RATE;

  #[test]
  #[ignore = "kb/issues/M-gtr-fitch-fit-frequencies-do-not-improve-with-more-sites.md"]
  fn test_sim_ancestral_gtr_error_falls_with_more_sites() -> Result<(), Report> {
    let log_ratios = run_replicates(BASE_SEED, REPLICATE_COUNT, replicate_log_ratios)?;

    let (w_ratios, pi_ratios): (Vec<f64>, Vec<f64>) =
      log_ratios.iter().map(|ratios| (ratios.w_rms, ratios.pi_max)).unzip();
    let checks = [
      decrease_check(
        "log ratio of the RMS error of log w at all and at the first sites",
        &w_ratios,
        FALSE_FAIL_RATE / 2.0,
      )?,
      decrease_check(
        "log ratio of the largest error of pi at all and at the first sites",
        &pi_ratios,
        FALSE_FAIL_RATE / 2.0,
      )?,
    ];
    require_passed(&checks)
  }

  mod helpers {
    use crate::ancestral::__tests__::test_sim_ancestral_calibration::tests::helpers::ancestral_case;
    use crate::ancestral::params::{AncestralParams, MethodAncestral};
    use crate::ancestral::pipeline;
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::GtrModelName;
    use crate::gtr::gtr::avg_transition;
    use crate::partition::marginal::sample::SampleMode;
    use crate::progress::NoopProgress;
    use crate::test_utils::{SimulatedCase, ancestral_input, first_columns};
    use eyre::Report;
    use itertools::Itertools;
    use ndarray::{Array1, Array2, aview1};
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_simulation::branch_lengths::exact;
    use treetime_simulation::substitution_model::NUCLEOTIDES;
    use treetime_utils::array::ndarray::max_or;
    use treetime_utils::{make_error, make_report};

    pub(super) const REPLICATE_COUNT: usize = 200;
    pub(super) const BASE_SEED: u64 = 1;
    const FEW_SITES: usize = 1000;
    const ALL_SITES: usize = 10_000;

    pub(super) fn replicate_log_ratios(seed: u64) -> Result<LogRatios, Report> {
      log_ratios(&ancestral_case(seed, ALL_SITES)?.simulate()?)
    }

    fn log_ratios(case: &SimulatedCase) -> Result<LogRatios, Report> {
      let lengths = exact(&case.tree.edge_durations(), &case.model)?;
      let few = gtr_errors(case, &lengths, FEW_SITES)?;
      let all = gtr_errors(case, &lengths, ALL_SITES)?;
      Ok(LogRatios {
        w_rms: (all.w_rms / few.w_rms).ln(),
        pi_max: (all.pi_max / few.pi_max).ln(),
      })
    }

    #[derive(Clone, Copy, Debug)]
    pub(super) struct LogRatios {
      pub w_rms: f64,
      pub pi_max: f64,
    }

    fn gtr_errors(
      case: &SimulatedCase,
      lengths: &BTreeMap<GraphEdgeKey, f64>,
      sites: usize,
    ) -> Result<GtrErrors, Report> {
      let params = AncestralParams {
        method: MethodAncestral::Marginal,
        model: GtrModelName::Infer,
        dense: Some(true),
        include_leaves: false,
        report_ambiguous: false,
        impute_missing_data: false,
        site_specific_gtr: false,
        seed: 0,
        sample_from_profile: SampleMode::Argmax,
      };
      let input = ancestral_input(&case.tree, &first_columns(&case.sequences.sequences, sites)?, lengths)?;
      let output = pipeline::run(&params, input, None, &NoopCancel, &NoopProgress, &NoopProgress)?;
      let gtr = output
        .results
        .iter()
        .exactly_one()
        .map_err(|error| make_report!("Expected one partition: {error}"))?
        .gtr
        .as_ref()
        .ok_or_else(|| make_report!("The partition has no inferred model"))?;
      if gtr.W().shape() != [NUCLEOTIDES.len(), NUCLEOTIDES.len()] {
        return make_error!("The inferred model has shape {:?}, expected 4x4", gtr.W().shape());
      }
      let w_inferred = normalized_w(gtr.W(), gtr.pi());
      let w_true = normalized_w(case.model.w(), case.model.pi());
      let log_w_errors = (0..NUCLEOTIDES.len())
        .tuple_combinations()
        .map(|(i, j)| (w_inferred[[i, j]] / w_true[[i, j]]).ln())
        .collect_vec();
      let w_rms = aview1(&log_w_errors)
        .mapv(|error| error * error)
        .mean()
        .ok_or_else(|| make_report!("No exchangeabilities"))?
        .sqrt();
      let pi_max = max_or(&(gtr.pi() - case.model.pi()).mapv(f64::abs), f64::NAN);
      Ok(GtrErrors { w_rms, pi_max })
    }

    fn normalized_w(w: &Array2<f64>, pi: &Array1<f64>) -> Array2<f64> {
      w / avg_transition(w, pi)
    }

    struct GtrErrors {
      w_rms: f64,
      pi_max: f64,
    }
  }
}

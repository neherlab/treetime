#[cfg(test)]
mod tests {
  use crate::test_utils::{recovery_check, require_passed, run_replicates};
  use eyre::Report;
  use helpers::{BASE_SEED, MARGIN, REPLICATE_COUNT, SIMULATED_SHAPES, replicate_shape};
  use itertools::Itertools;
  use treetime_simulation::FALSE_FAIL_RATE;
  use treetime_utils::cast::ToF64Exact;

  #[test]
  fn test_sim_ancestral_gamma_shape_recovers_the_simulated_shape() -> Result<(), Report> {
    let checks = SIMULATED_SHAPES
      .iter()
      .map(|&shape| {
        let estimates = run_replicates(BASE_SEED, REPLICATE_COUNT, |seed| replicate_shape(seed, shape))?;
        let log_ratios = estimates.iter().map(|estimate| (estimate / shape).ln()).collect_vec();
        recovery_check(
          &format!("log ratio of the estimated and the simulated gamma shape {shape}"),
          &log_ratios,
          MARGIN,
          FALSE_FAIL_RATE / SIMULATED_SHAPES.len().to_f64_exact(),
        )
      })
      .try_collect::<_, Vec<_>, _>()?;
    require_passed(&checks)
  }

  mod helpers {
    use crate::ancestral::params::{AncestralParams, MethodAncestral};
    use crate::ancestral::pipeline::run;
    use crate::cancel::NoopCancel;
    use crate::gtr::gamma::GammaRates;
    use crate::gtr::get_gtr::{GtrModel, JC69Params};
    use crate::partition::marginal::sample::SampleMode;
    use crate::progress::NoopProgress;
    use crate::test_utils::{TestDataset, ancestral_input_for_tree, simulated_sequence_source};
    use eyre::Report;
    use itertools::Itertools;
    use treetime_graph::key_map::EdgeMap;
    use treetime_io::nwk::nwk_read;
    use treetime_simulation::evolve::{draw_site_rates, evolve_sequences_at_site_rates};
    use treetime_simulation::substitution_model::SubstitutionModel;
    use treetime_utils::sync::random::get_random_number_generator;
    use treetime_utils::{make_internal_error, make_internal_report, make_report};

    pub(super) const SIMULATED_SHAPES: [f64; 2] = [0.5, 2.0];
    pub(super) const REPLICATE_COUNT: usize = 371;
    pub(super) const BASE_SEED: u64 = 1;
    pub(super) const MARGIN: f64 = 0.05;
    const SITES: usize = 3000;
    const CATEGORIES: usize = 4;
    const TREE_SCALE: f64 = 20.0;

    pub(super) fn replicate_shape(seed: u64, shape: f64) -> Result<f64, Report> {
      let mut parse = nwk_read(TestDataset::Zika86.tree().as_bytes())?;
      let durations: EdgeMap<f64> = parse
        .branch_lengths
        .iter()
        .map(|(edge, &length)| {
          length
            .map(|length| (edge, length * TREE_SCALE))
            .ok_or_else(|| make_report!("Edge {edge} of the Zika tree has no branch length"))
        })
        .try_collect()?;
      parse.branch_lengths = durations
        .iter()
        .map(|(edge, &duration)| (edge, Some(duration)))
        .collect();
      let mut rng = get_random_number_generator(seed);
      let site_rates = draw_site_rates(GammaRates::new(CATEGORIES, shape)?.rates(), SITES, &mut rng)?;
      let sequences = evolve_sequences_at_site_rates(
        &parse.graph,
        &durations,
        &SubstitutionModel::jc69(1.0)?,
        &site_rates,
        &mut rng,
      )?;
      let source = simulated_sequence_source(&parse.graph, &parse.names, &sequences.sequences)?;
      let mut input = ancestral_input_for_tree(parse, vec![source], &GtrModel::JC69(JC69Params::default()))?;
      input.sequences.models.gamma_categories = CATEGORIES;
      let params = AncestralParams {
        method: MethodAncestral::Marginal,
        dense: Some(false),
        include_leaves: false,
        report_ambiguous: false,
        impute_missing_data: false,
        seed,
        sample_from_profile: SampleMode::Argmax,
      };
      let output = run(&params, input, None, &NoopCancel, &NoopProgress, &NoopProgress)?;
      let [partition] = output.partitions.as_slice() else {
        return make_internal_error!("A run on one alignment gave {} partitions", output.partitions.len());
      };
      partition
        .model
        .as_ref()
        .and_then(|model| model.gamma().shape())
        .ok_or_else(|| make_internal_report!("The partition has no estimated gamma shape"))
    }
  }
}

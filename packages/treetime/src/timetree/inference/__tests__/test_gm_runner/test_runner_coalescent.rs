#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{SequenceMode, marginal_setup};
  use crate::coalescent::mean_field::{CoalescentPrior, LineageExpectation};
  use crate::optimize::indel::IndelPolicy;
  use crate::timetree::inference::result::PosteriorSummaries;
  use eyre::Report;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_distribution::Distribution;

  #[rustfmt::skip]
  #[rstest]
  #[case::tc_0_1( 0.1)]
  #[case::tc_1_0( 1.0)]
  #[case::tc_10_0(10.0)]
  #[trace]
  fn test_runner_coalescent_completes(#[case] tc: f64) -> Result<(), Report> {
    let setup = marginal_setup("flu_h3n2_20", SequenceMode::Dense)?;
    let survival = PosteriorSummaries {
      lineage_survival: true,
      ..PosteriorSummaries::default()
    };
    let prior_free = setup.run(IndelPolicy::Include, None, survival)?;
    let expectation = LineageExpectation::of_inference(&setup.graph, &prior_free)?;
    let coalescent = CoalescentPrior::new(&expectation, &Distribution::constant(tc))?;

    let inference = setup.run(IndelPolicy::Include, Some(&coalescent), PosteriorSummaries::default())?;

    let times = setup.node_times(&inference);
    assert_eq!(setup.graph.num_nodes(), times.len(), "every node must be dated");
    let non_finite: BTreeMap<&String, &f64> = times.iter().filter(|(_, time)| !time.is_finite()).collect();
    assert_eq!(BTreeMap::new(), non_finite, "every node time must be finite");
    Ok(())
  }
}

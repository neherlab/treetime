#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::clock::clock_model::{ClockModel, ClockModelStats, RegressionStats};
  use crate::clock::date_constraints::DateConstraints;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::partition::create::{Representation, build_marginal_partition};
  use crate::progress::NoopProgress;
  use crate::test_utils::leaf_seq_inputs;
  use crate::test_utils::{find_node_key_by_name, point_date_constraints};
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::confidence::{
    RateSusceptibility, compute_rate_susceptibility, date_uncertainty_due_to_rate, determine_rate_std,
    quantile_to_zscore,
  };
  use crate::timetree::inference::bad_branches::bad_leaves;
  use crate::timetree::inference::runner::{TimeInferenceInputs, run_timetree};
  use crate::timetree::optimization::relaxed_clock::unit_gammas;
  use approx::assert_relative_eq;
  use eyre::Report;
  use indoc::indoc;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_grid::MaxGridPoints;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::assert_error;

  const FIXTURE_CLOCK_RATE: f64 = 0.002;

  const RATE_STD_FRACTION: f64 = 0.5;

  const FIXTURE_DATES: [(&str, f64); 4] = [("A", 2010.0), ("B", 2012.0), ("C", 2011.0), ("D", 2013.0)];

  #[rustfmt::skip]
  #[rstest]
  #[case::lower_2_5pct(0.025,  -1.959964)]
  #[case::lower_5pct(  0.05,   -1.644854)]
  #[case::median(       0.5,    0.0      )]
  #[case::upper_95pct(  0.95,   1.644854 )]
  #[case::upper_97_5pct(0.975,  1.959964 )]
  #[case::boundary_zero(0.0,    0.0      )]
  #[case::boundary_one( 1.0,    0.0      )]
  #[trace]
  fn test_quantile_to_zscore(#[case] p: f64, #[case] expected: f64) {
    let z = quantile_to_zscore(p);
    assert_relative_eq!(z, expected, epsilon = 1e-4);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::symmetric_1sigma(   [9.0, 10.0, 11.0],  (0.025, 0.975), (8.040036, 11.959964))]
  #[case::symmetric_narrow(   [9.5, 10.0, 10.5],  (0.025, 0.975), (9.020018, 10.979982))]
  #[case::asymmetric(         [8.0, 10.0, 11.0],  (0.025, 0.975), (6.080072, 11.959964))]
  #[case::boundary_quantiles( [9.0, 10.0, 11.0],  (0.0,   1.0),   (10.0,     10.0)     )]
  #[case::equal_dates(        [10.0, 10.0, 10.0], (0.025, 0.975), (10.0,     10.0)     )]
  #[trace]
  fn test_date_uncertainty_due_to_rate(
    #[case] dates: [f64; 3],
    #[case] interval: (f64, f64),
    #[case] (expected_lower, expected_upper): (f64, f64),
  ) {
    let (lower, upper) = date_uncertainty_due_to_rate(dates, interval);
    assert_relative_eq!(lower, expected_lower, epsilon = 1e-4);
    assert_relative_eq!(upper, expected_upper, epsilon = 1e-4);
  }

  #[test]
  fn test_determine_rate_std_explicit_clock_std_dev() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(Some(0.001), false, &clock_model, &NoopProgress).unwrap();
    assert_relative_eq!(result.unwrap(), 0.001);
  }

  #[test]
  fn test_determine_rate_std_rejects_negative() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(Some(-0.001), false, &clock_model, &NoopProgress);
    assert_error!(result, "--clock-std-dev must be positive, got -0.001");
  }

  #[test]
  fn test_determine_rate_std_rejects_zero() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(Some(0.0), false, &clock_model, &NoopProgress);
    assert_error!(result, "--clock-std-dev must be positive, got 0");
  }

  #[test]
  fn test_determine_rate_std_none_without_covariation() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(None, false, &clock_model, &NoopProgress).unwrap();
    assert!(result.is_none());
  }

  #[test]
  fn test_determine_rate_std_from_covariance_matrix() {
    let clock_model = ClockModel::for_testing_with_stats(
      0.003,
      0.0,
      ClockModelStats::Estimated(RegressionStats {
        chisq: 0.0,
        r_val: 0.9,
        hessian: array![[1.0, 0.0], [0.0, 1.0]],
        cov: array![[1e-6, 0.0], [0.0, 1.0]],
      }),
    );
    let result = determine_rate_std(None, true, &clock_model, &NoopProgress).unwrap();
    assert_relative_eq!(result.unwrap(), 1e-3, epsilon = 1e-10);
  }

  #[test]
  fn test_determine_rate_std_none_for_fixed_clock_with_covariation() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(None, true, &clock_model, &NoopProgress).unwrap();
    assert!(result.is_none());
  }

  #[test]
  fn test_compute_rate_susceptibility_brackets_each_date_and_keeps_the_central_inference() -> Result<(), Report> {
    let fixture = helpers::SusceptibilityFixture::new()?;
    let inputs = fixture.inputs();
    let rate_std = RATE_STD_FRACTION * FIXTURE_CLOCK_RATE;

    let RateSusceptibility { dates, central } = compute_rate_susceptibility(&inputs, None, rate_std, &NoopProgress)?;

    let direct = run_timetree(&inputs, None, &NoopProgress)?;
    assert_eq!(direct, central);

    let dated_nodes: BTreeSet<GraphNodeKey> = central
      .posterior
      .iter()
      .filter_map(|(key, posterior)| posterior.time.map(|_| *key))
      .collect();
    assert_eq!(dated_nodes, dates.keys().copied().collect::<BTreeSet<_>>());
    for (key, triple) in &dates {
      assert!(
        triple[0] <= triple[1] && triple[1] <= triple[2],
        "dates of node {key} must be sorted: {triple:?}"
      );
      let central_time = central.posterior[key].time.expect("dated node");
      assert!(
        triple.contains(&central_time),
        "dates of node {key} must contain the central time {central_time}: {triple:?}"
      );
    }
    for (name, date) in FIXTURE_DATES {
      assert_eq!(
        [date, date, date].map(f64::to_bits),
        dates[&fixture.key(name)].map(f64::to_bits),
        "exact date of leaf {name}"
      );
    }
    let root = dates[&fixture.key("root")];
    assert!(
      root[0] < root[2],
      "a slower and a faster clock must move the root apart: {root:?}"
    );
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) struct SusceptibilityFixture {
      graph: Graph,
      names: BTreeMap<GraphNodeKey, Option<String>>,
      constraints: DateConstraints,
      leaf_bad_branches: BTreeMap<GraphNodeKey, bool>,
      gammas: BTreeMap<GraphEdgeKey, f64>,
      branch_model: BranchModel,
      branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
      clock_model: ClockModel,
    }

    impl SusceptibilityFixture {
      pub(super) fn new() -> Result<Self, Report> {
        let nwk_parsed = nwk_read(b"((A:0.01,B:0.02)AB:0.01,(C:0.015,D:0.01)CD:0.02)root;".as_slice())?;
        let names = nwk_parsed.names();
        let graph = nwk_parsed.graph;
        let branch_lengths = nwk_parsed.branch_lengths;
        let alphabet = Alphabet::new(AlphabetName::Nuc)?;
        let aln: Vec<AlignmentRecord> = fasta_read(
          indoc! {"
            >A
            ACGTACGTACGTACGTACGTACGT
            >B
            ACGTACGTACGAACGTACGTACCT
            >C
            ACGTTCGTACGTACGTACGTACGT
            >D
            ACGTTCGTACGTACGTAGGTACGT
          "}
          .as_bytes(),
          &alphabet,
        )?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
        let reconstruction = build_marginal_partition(
          Representation::Dense,
          GtrModelName::JC69,
          &graph,
          alphabet,
          leaf_seq_inputs(&graph, &names, aln),
          &branch_lengths_or_zero(&branch_lengths),
          &NoopProgress,
        )?
        .marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?
        .0;
        let constraints = point_date_constraints(&graph, &names, &FIXTURE_DATES);
        let leaf_bad_branches = bad_leaves(&graph, &constraints, &BTreeSet::new());
        let gammas = unit_gammas(&graph);
        Ok(Self {
          graph,
          names,
          constraints,
          leaf_bad_branches,
          gammas,
          branch_model: BranchModel::Marginal(reconstruction),
          branch_lengths,
          clock_model: ClockModel::for_testing(FIXTURE_CLOCK_RATE, 0.0),
        })
      }

      pub(super) fn key(&self, name: &str) -> GraphNodeKey {
        find_node_key_by_name(&self.graph, &self.names, name).expect("fixture node must exist")
      }

      pub(super) fn inputs(&self) -> TimeInferenceInputs<'_> {
        TimeInferenceInputs {
          graph: &self.graph,
          date_constraints: &self.constraints,
          leaf_bad_branches: &self.leaf_bad_branches,
          gammas: &self.gammas,
          branch_model: &self.branch_model,
          branch_lengths: &self.branch_lengths,
          names: &self.names,
          clock_model: &self.clock_model,
          clock_rate_fixed: false,
          no_indels: false,
          max_grid_points: MaxGridPoints::default(),
        }
      }
    }
  }
}

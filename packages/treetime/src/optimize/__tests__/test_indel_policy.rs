#[cfg(test)]
mod tests {
  use crate::optimize::__tests__::test_optimize_indel::tests::inject_indels_on_first_edge;
  use crate::optimize::dispatch::guess_initial_branch_lengths;
  use crate::optimize::gather::InitialGuessEvidence;
  use crate::optimize::indel::{EdgeIndels, IndelPolicy, IndelTerm};
  use crate::optimize::params::ExistingBranchLengths;
  use crate::test_utils::deletion;
  use eyre::Report;
  use helpers::{BRANCH_LENGTH, Fixture};
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use statrs::function::factorial::ln_factorial;
  use treetime_primitives::Seq;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::include(IndelPolicy::Include, 4.0)]
  #[case::exclude(IndelPolicy::Exclude, 0.0)]
  #[trace]
  fn test_indel_policy_rate_is_poisson_estimate_or_zero(
    #[case] policy: IndelPolicy,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let fixture = Fixture::new()?;
    let indels = fixture.edge_indels(policy);
    pretty_assert_ulps_eq!(expected, indels.rate());
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::include_indel_edge(IndelPolicy::Include, 0, (3, 4.0))]
  #[case::include_other_edge(IndelPolicy::Include, 1, (0, 4.0))]
  #[case::exclude_indel_edge(IndelPolicy::Exclude, 0, (0, 0.0))]
  #[trace]
  fn test_indel_policy_term_carries_edge_count_and_rate(
    #[case] policy: IndelPolicy,
    #[case] edge_index: usize,
    #[case] (count, rate): (usize, f64),
  ) -> Result<(), Report> {
    let fixture = Fixture::new()?;
    let term = fixture.edge_indels(policy).term(fixture.edge(edge_index));
    assert_eq!(count, term.count);
    pretty_assert_ulps_eq!(rate, term.rate);
    pretty_assert_ulps_eq!(0.0, term.held_length);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::include(IndelPolicy::Include, (3, 4.0))]
  #[case::exclude(IndelPolicy::Exclude, (0, 0.0))]
  #[trace]
  fn test_indel_policy_chain_term_sums_counts_and_keeps_held_length(
    #[case] policy: IndelPolicy,
    #[case] (count, rate): (usize, f64),
  ) -> Result<(), Report> {
    let fixture = Fixture::new()?;
    let chain = [fixture.edge(0), fixture.edge(1)];
    let term = fixture.edge_indels(policy).chain_term(&chain, 0.25);
    assert_eq!(count, term.count);
    pretty_assert_ulps_eq!(rate, term.rate);
    pretty_assert_ulps_eq!(0.25, term.held_length);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::include(IndelPolicy::Include, (3.0 * 2.0_f64.ln() - 2.0 - ln_factorial(3), 2.0, -12.0))]
  #[case::exclude(IndelPolicy::Exclude, (0.0,                                         0.0,  0.0))]
  #[trace]
  fn test_indel_policy_objective_term_is_poisson_or_zero(
    #[case] policy: IndelPolicy,
    #[case] (log_lh, derivative, second_derivative): (f64, f64, f64),
  ) -> Result<(), Report> {
    let fixture = Fixture::new()?;
    let metrics = fixture.edge_indels(policy).term(fixture.edge(0)).metrics(0.5)?;
    pretty_assert_abs_diff_eq!(log_lh, metrics.log_lh.value(), epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(derivative, metrics.derivative, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(second_derivative, metrics.second_derivative, epsilon = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::include(IndelPolicy::Include, 3.0 * 0.5_f64.ln() - 3.0 - ln_factorial(3))]
  #[case::exclude(IndelPolicy::Exclude, 0.0)]
  #[trace]
  fn test_indel_policy_total_log_lh_sums_poisson_terms_or_is_zero(
    #[case] policy: IndelPolicy,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let fixture = Fixture::new()?;
    let total = fixture.edge_indels(policy).total_log_lh(&fixture.graph, &fixture.branch_lengths)?;
    pretty_assert_abs_diff_eq!(expected, total.value(), epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_indel_policy_exclude_total_log_lh_accepts_zero_length_indel_edge() -> Result<(), Report> {
    let mut fixture = Fixture::new()?;
    fixture.branch_lengths.insert(fixture.edge(0), Some(0.0));
    let total = fixture
      .edge_indels(IndelPolicy::Exclude)
      .total_log_lh(&fixture.graph, &fixture.branch_lengths)?;
    pretty_assert_ulps_eq!(0.0, total.value());
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::include_zero_length_indel_edge(IndelPolicy::Include, 0.0,           true )]
  #[case::include_positive_indel_edge(   IndelPolicy::Include, BRANCH_LENGTH, false)]
  #[case::exclude_zero_length_indel_edge(IndelPolicy::Exclude, 0.0,           false)]
  #[trace]
  fn test_indel_policy_zero_length_indel_edge_detection(
    #[case] policy: IndelPolicy,
    #[case] indel_edge_length: f64,
    #[case] expected: bool,
  ) -> Result<(), Report> {
    let mut fixture = Fixture::new()?;
    fixture.branch_lengths.insert(fixture.edge(0), Some(indel_edge_length));
    let indels = fixture.edge_indels(policy);
    assert_eq!(expected, indels.has_indel_edge_with_zero_length(&fixture.graph, &fixture.branch_lengths));
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::include(IndelPolicy::Include, array![0.75, 0.0, 0.0, 0.0, 0.0, 0.0])]
  #[case::exclude(IndelPolicy::Exclude, array![0.0,  0.0, 0.0, 0.0, 0.0, 0.0])]
  #[trace]
  fn test_indel_policy_initial_guess(#[case] policy: IndelPolicy, #[case] expected: Array1<f64>) -> Result<(), Report> {
    let mut fixture = Fixture::new()?;
    let indels = fixture.edge_indels(policy);
    guess_initial_branch_lengths(
      &fixture.graph,
      &indels,
      &InitialGuessEvidence::gather(&fixture.graph, &fixture.reconstruction)?,
      ExistingBranchLengths::Overwrite,
      &mut fixture.branch_lengths,
    )?;
    pretty_assert_abs_diff_eq!(expected, fixture.edge_lengths(), epsilon = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_indels(IndelTerm { count: 0, rate: 4.0, held_length: 0.0 })]
  #[case::zero_rate(IndelTerm { count: 3, rate: 0.0, held_length: 0.0 })]
  #[trace]
  fn test_indel_policy_term_without_poisson_optimum_has_no_length_estimate(#[case] term: IndelTerm) {
    assert_eq!(None, term.length_estimate());
  }

  #[test]
  fn test_indel_policy_term_length_estimate_is_poisson_optimum() {
    let term = IndelTerm {
      count: 3,
      rate: 4.0,
      held_length: 0.0,
    };
    let expected = 0.75;
    pretty_assert_ulps_eq!(expected, term.length_estimate().expect("positive count and rate"));
  }

  #[test]
  fn test_indel_policy_include_recount_keeps_rate_and_reads_new_counts() -> Result<(), Report> {
    let mut fixture = Fixture::new()?;
    let indels = fixture.edge_indels(IndelPolicy::Include);
    inject_indels_on_first_edge(
      &fixture.graph,
      &mut fixture.reconstruction,
      &[
        deletion((0, 1), Seq::try_from_str("A")?),
        deletion((4, 5), Seq::try_from_str("A")?),
      ],
    );
    let recounted = indels.recount(&fixture.graph, &fixture.reconstruction);
    let term = recounted.term(fixture.edge(0));
    assert_eq!(2, term.count);
    pretty_assert_ulps_eq!(4.0, term.rate);
    Ok(())
  }

  #[test]
  fn test_indel_policy_exclude_recount_stays_excluded() -> Result<(), Report> {
    let fixture = Fixture::new()?;
    let recounted = fixture
      .edge_indels(IndelPolicy::Exclude)
      .recount(&fixture.graph, &fixture.reconstruction);
    assert!(matches!(recounted, EdgeIndels::Excluded), "got {recounted:?}");
    Ok(())
  }

  mod helpers {
    use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::TREE_NEWICK;
    use crate::optimize::__tests__::test_optimize_indel::tests::{
      inject_indels_on_first_edge, setup_identical_reconstruction,
    };
    use crate::optimize::indel::{EdgeIndels, IndelPolicy};
    use crate::partition::marginal::reconstruction::MarginalReconstruction;
    use crate::test_utils::deletion;
    use eyre::Report;
    use ndarray::Array1;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::Seq;

    pub(super) const BRANCH_LENGTH: f64 = 0.125;

    pub(super) struct Fixture {
      pub graph: Graph,
      pub reconstruction: MarginalReconstruction,
      pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
    }

    impl Fixture {
      pub(super) fn new() -> Result<Self, Report> {
        let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
        let names = nwk_parsed.names();
        let graph = nwk_parsed.graph;
        let mut branch_lengths = nwk_parsed.branch_lengths;
        let mut reconstruction = setup_identical_reconstruction(&graph, &names, &mut branch_lengths)?;
        inject_indels_on_first_edge(
          &graph,
          &mut reconstruction,
          &[
            deletion((0, 1), Seq::try_from_str("A")?),
            deletion((4, 5), Seq::try_from_str("A")?),
            deletion((8, 9), Seq::try_from_str("A")?),
          ],
        );
        let branch_lengths = branch_lengths
          .into_keys()
          .map(|edge_key| (edge_key, Some(BRANCH_LENGTH)))
          .collect();
        Ok(Self {
          graph,
          reconstruction,
          branch_lengths,
        })
      }

      pub(super) fn edge_indels(&self, policy: IndelPolicy) -> EdgeIndels {
        policy.edge_indels(&self.graph, &self.reconstruction, &self.branch_lengths)
      }

      pub(super) fn edge(&self, index: usize) -> GraphEdgeKey {
        self
          .graph
          .get_edges()
          .nth(index)
          .expect("edge index within the tree")
          .key()
      }

      pub(super) fn edge_lengths(&self) -> Array1<f64> {
        self
          .graph
          .get_edges()
          .map(|edge| self.branch_lengths[&edge.key()].expect("every edge has a length"))
          .collect()
      }
    }
  }
}

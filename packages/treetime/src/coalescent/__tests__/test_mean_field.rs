#[cfg(test)]
mod tests {
  use crate::coalescent::mean_field::{CoalescentPrior, LineageExpectation, NodeSurvival, SegmentStatistics};
  use crate::test_utils::{empty_time_inference, find_node_key_by_name, point_coalescent_prior};
  use crate::timetree::inference::result::{NodePosterior, PosteriorSummaries, TimeInference};
  use eyre::Report;
  use ndarray::Array1;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::root_far_back(   2.0, 1990.0)]
  #[case::root_near_tips(  2.0, 2009.0)]
  #[case::long_time_scale(25.0, 2000.0)]
  #[trace]
  fn test_mean_field_cherry_root_density_is_kingman(#[case] tc: f64, #[case] root_time: f64) -> Result<(), Report> {
    let tip = 2010.0;
    let parsed = nwk_read(b"(A:1,B:1)root;".as_slice())?;
    let names = parsed.names();
    let graph = parsed.graph;
    let prior = point_coalescent_prior(&graph, &names, &[("A", tip), ("B", tip), ("root", root_time)], tc)?;
    let root = find_node_key_by_name(&graph, &names, "root").expect("root must exist");

    let reference = 2005.0;
    for t in [1960.0, 1995.0, 2004.0, 2009.9] {
      let actual = prior.node_contribution(root, 2, t)? - prior.node_contribution(root, 2, reference)?;
      pretty_assert_abs_diff_eq!((reference - t) / tc, actual, epsilon = 1e-9);
    }
    Ok(())
  }

  #[test]
  fn test_mean_field_internal_slope_counts_the_lineages_without_the_own_merger() -> Result<(), Report> {
    let parsed = nwk_read(b"((A:1,B:1)I:1,C:1)root;".as_slice())?;
    let names = parsed.names();
    let graph = parsed.graph;
    let times = [
      ("A", 2015.0),
      ("B", 2014.0),
      ("C", 2013.0),
      ("I", 2011.0),
      ("root", 2005.0),
    ];
    let prior = point_coalescent_prior(&graph, &names, &times, 1.0)?;
    let internal = find_node_key_by_name(&graph, &names, "I").expect("I must exist");
    let slope = |a: f64, b: f64| -> Result<f64, Report> {
      Ok((prior.node_contribution(internal, 2, b)? - prior.node_contribution(internal, 2, a)?) / (b - a))
    };

    pretty_assert_abs_diff_eq!(1.0 - 3.0, slope(2006.0, 2012.0)?, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(1.0 - 2.0, slope(2013.2, 2013.8)?, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(1.0 - 1.0, slope(2014.2, 2014.8)?, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(1.0 - 2.0, slope(1990.0, 2000.0)?, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_mean_field_segment_statistics_of_point_times_are_those_of_the_lineage_counts() -> Result<(), Report> {
    let parsed = nwk_read(b"((A:1,B:1)I:1,C:1)root;".as_slice())?;
    let names = parsed.names();
    let graph = parsed.graph;
    let times = [
      ("A", 2015.0),
      ("B", 2014.0),
      ("C", 2013.0),
      ("I", 2011.0),
      ("root", 2005.0),
    ];
    let expectation = helpers::expectation(&graph, &names, &times, &[])?;

    let SegmentStatistics { pair_integral, mergers } = expectation.segment_statistics(&[2005.0, 2015.0]);

    pretty_assert_abs_diff_eq!(13.0, pair_integral[0], epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(2.0, mergers[0], epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_mean_field_segment_statistics_integrate_over_the_root_posterior() -> Result<(), Report> {
    let parsed = nwk_read(b"(A:1,B:1)root;".as_slice())?;
    let names = parsed.names();
    let graph = parsed.graph;
    let flat_root = Distribution::function(Array1::linspace(2000.0, 2004.0, 401), Array1::zeros(401))?;
    let expectation = helpers::expectation(&graph, &names, &[("A", 2010.0), ("B", 2010.0)], &[("root", flat_root)])?;

    let SegmentStatistics { pair_integral, mergers } = expectation.segment_statistics(&[2000.0, 2010.0]);

    pretty_assert_abs_diff_eq!(8.0, pair_integral[0], epsilon = 1e-11);
    pretty_assert_abs_diff_eq!(1.0, mergers[0], epsilon = 1e-9);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::far_truncation(   2000.0, 2.0, 2000.0)]
  #[case::strong_truncation(2015.0, 1.0, 2015.0 + statrs::distribution::ContinuousCDF::inverse_cdf(&statrs::distribution::Normal::standard(), statrs::distribution::ContinuousCDF::cdf(&statrs::distribution::Normal::standard(), -5.0) / 2.0))]
  #[case::exponential_tail( 2050.0, 1.0, 2010.0 - std::f64::consts::LN_2 / 40.0)]
  #[trace]
  fn test_mean_field_gaussian_seed_median_is_that_of_the_truncated_normal(
    #[case] mean: f64,
    #[case] std: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let (_, _, expectation) = helpers::gaussian_cherry(mean, std, &[])?;

    let SegmentStatistics { mergers, .. } = expectation.segment_statistics(&[expected - 100.0, expected, expected + 100.0]);

    pretty_assert_abs_diff_eq!(0.5, mergers[0], epsilon = 1e-6);
    Ok(())
  }

  #[test]
  fn test_mean_field_gaussian_seed_splits_the_merger_at_the_mean() -> Result<(), Report> {
    let (_, _, expectation) = helpers::gaussian_cherry(2000.0, 2.0, &[])?;

    let SegmentStatistics { mergers, .. } = expectation.segment_statistics(&[1990.0, 2000.0, 2020.0]);

    pretty_assert_abs_diff_eq!(0.5, mergers[0], epsilon = 1e-6);
    pretty_assert_abs_diff_eq!(0.5, mergers[1], epsilon = 1e-6);
    Ok(())
  }

  #[test]
  fn test_mean_field_gaussian_seed_truncates_at_the_earliest_included_tip() -> Result<(), Report> {
    let (_, _, expectation) = helpers::gaussian_cherry_dated(2005.0, 1.0, [2010.0, 2000.0], &["B"])?;

    let SegmentStatistics { mergers, .. } = expectation.segment_statistics(&[1990.0, 2000.0, 2020.0]);

    pretty_assert_abs_diff_eq!(1.0, mergers[1], epsilon = 1e-6);
    Ok(())
  }

  #[test]
  fn test_mean_field_own_share_is_removed_at_its_weight_in_the_expectation() -> Result<(), Report> {
    let parsed = nwk_read(b"(A:1,B:1,C:1)root;".as_slice())?;
    let names = parsed.names();
    let graph = parsed.graph;
    let tc = 2.0;
    let prior = point_coalescent_prior(
      &graph,
      &names,
      &[("A", 2010.0), ("B", 2010.0), ("C", 2010.0), ("root", 2000.0)],
      tc,
    )?;
    let root = find_node_key_by_name(&graph, &names, "root").expect("root must exist");
    let slope = |x: f64| -> Result<f64, Report> {
      let h = 1e-3;
      Ok((prior.node_contribution(root, 2, x + h)? - prior.node_contribution(root, 2, x - h)?) / (2.0 * h))
    };

    pretty_assert_abs_diff_eq!(-2.0 / tc, slope(1995.0)?, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(-2.0 / tc, slope(2005.0)?, epsilon = 1e-9);
    Ok(())
  }

  #[test]
  fn test_mean_field_with_topology_equals_a_fresh_expectation() -> Result<(), Report> {
    let parsed = nwk_read(b"((A:1,B:1)N:1,C:1)root;".as_slice())?;
    let names = parsed.names();
    let graph = parsed.graph;
    let key = |name| find_node_key_by_name(&graph, &names, name).expect("node must exist");
    let times = [
      ("A", 2010.0),
      ("B", 2009.0),
      ("C", 2012.0),
      ("root", 2000.0),
      ("N", 2005.0),
    ];
    let point = |time: f64| NodePosterior {
      likely_time: Some(time),
      time: Some(time),
      ..NodePosterior::default()
    };
    let full: BTreeMap<GraphNodeKey, NodePosterior> = times.iter().map(|(n, t)| (key(n), point(*t))).collect();
    let mut without_n = full.clone();
    without_n.remove(&key("N"));
    let bad_branches = BTreeMap::new();
    let tc = Distribution::constant(3.0);

    let extended = LineageExpectation::from_posterior(&graph, &without_n, &bad_branches)?
      .with_topology(&graph, &BTreeMap::from([(key("N"), 2005.0)]))?;
    let fresh = LineageExpectation::from_posterior(&graph, &full, &bad_branches)?;
    let (extended, fresh) = (
      CoalescentPrior::new(&extended, &tc)?,
      CoalescentPrior::new(&fresh, &tc)?,
    );

    for (name, n_children) in [("root", 2), ("N", 2), ("A", 0)] {
      for x in [1995.0, 2003.0, 2007.0] {
        pretty_assert_abs_diff_eq!(
          fresh.node_contribution(key(name), n_children, x)?,
          extended.node_contribution(key(name), n_children, x)?,
          epsilon = 1e-12
        );
      }
    }
    Ok(())
  }

  #[test]
  fn test_mean_field_of_inference_rejects_an_inference_without_survival_functions() -> Result<(), Report> {
    let graph = nwk_read(b"(A:1,B:1)root;".as_slice())?.graph;
    let dated = NodePosterior {
      likely_time: Some(2010.0),
      time: Some(2010.0),
      ..NodePosterior::default()
    };
    let inference = TimeInference {
      posterior: graph.get_nodes().map(|node| (node.key(), dated.clone())).collect(),
      posterior_summaries: PosteriorSummaries {
        regions: true,
        lineage_survival: false,
      },
      ..empty_time_inference(&graph)
    };

    let result = LineageExpectation::of_inference(&graph, &inference);

    assert_error!(
      result,
      "The time inference holds no survival functions, so it cannot give expected lineage counts. This is an \
       internal error. Please report it to developers."
    );
    Ok(())
  }

  mod helpers {
    use crate::clock::gaussian::NodeTimeGaussian;

    pub(super) fn gaussian_cherry(
      mean: f64,
      std: f64,
      excluded: &[&str],
    ) -> Result<(Graph, BTreeMap<GraphNodeKey, Option<String>>, Arc<LineageExpectation>), Report> {
      gaussian_cherry_dated(mean, std, [2010.0, 2010.0], excluded)
    }

    pub(super) fn gaussian_cherry_dated(
      mean: f64,
      std: f64,
      dates: [f64; 2],
      excluded: &[&str],
    ) -> Result<(Graph, BTreeMap<GraphNodeKey, Option<String>>, Arc<LineageExpectation>), Report> {
      let parsed = nwk_read(b"(A:1,B:1)root;".as_slice())?;
      let names = parsed.names();
      let graph = parsed.graph;
      let key = |name| find_node_key_by_name(&graph, &names, name).expect("node must exist");
      let marginals = BTreeMap::from([(key("root"), NodeTimeGaussian { mean, std })]);
      let tip_times = BTreeMap::from([
        (key("A"), Some(dates[0])),
        (key("B"), Some(dates[1])),
        (key("root"), None),
      ]);
      let excluded = excluded.iter().map(|name| (key(name), true)).collect();
      let expectation = LineageExpectation::from_gaussian(&graph, &marginals, &tip_times, &excluded)?;
      Ok((graph, names, expectation))
    }

    use super::*;

    pub(super) fn expectation(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      points: &[(&str, f64)],
      distributions: &[(&str, Distribution<NegLog>)],
    ) -> Result<Arc<LineageExpectation>, Report> {
      let key = |name: &str| find_node_key_by_name(graph, names, name).expect("node must exist");
      let posterior = points
        .iter()
        .map(|&(name, time)| {
          let entry = NodePosterior {
            time: Some(time),
            ..NodePosterior::default()
          };
          (key(name), entry)
        })
        .chain(distributions.iter().map(|(name, distribution)| {
          let entry = NodePosterior {
            survival: NodeSurvival::of_posterior(distribution)
              .expect("fixture distribution must have a survival function")
              .map(Arc::new),
            ..NodePosterior::default()
          };
          (key(name), entry)
        }))
        .collect();
      let bad_branches = graph.get_nodes().map(|node| (node.key(), false)).collect();
      LineageExpectation::from_posterior(graph, &posterior, &bad_branches)
    }
  }
}

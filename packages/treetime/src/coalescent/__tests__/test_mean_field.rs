#[cfg(test)]
mod tests {
  use crate::clock::gaussian::NodeTimeGaussian;
  use crate::coalescent::mean_field::{CoalescentPrior, LineageExpectation, NodeSurvival, SegmentStatistics};
  use crate::test_utils::{empty_time_inference, find_node_key_by_name, point_coalescent_prior};
  use crate::timetree::inference::result::{NodePosterior, PosteriorSummaries, TimeInference};
  use eyre::Report;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::{BoundaryBehavior, Distribution, DistributionFunction, NegLog};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_graph::node_names::NodeNames;
  use treetime_grid::{Side, SoftTailLaw};
  use treetime_io::nwk::nwk_read;
  use treetime_utils::cast::ToF64Exact;
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
    let names = parsed.names.clone();
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

  #[rustfmt::skip]
  #[rstest]
  #[case::binary(
    "((A:1,B:1)N:1,C:1)root;",
    &[("A", 2010.0), ("B", 2009.0), ("C", 2012.0), ("N", 2005.0), ("root", 2000.0)],
    (1.0 + 4.0 * 3.0 + 5.0) / 3.0 + 2.0 * 3.0_f64.ln(),
  )]
  #[case::polytomy(
    "(A:1,B:1,C:1)root;",
    &[("A", 2010.0), ("B", 2009.0), ("C", 2012.0), ("root", 2000.0)],
    (1.0 + 9.0 * 3.0) / 3.0 + 2.0 * 3.0_f64.ln(),
  )]
  #[trace]
  fn test_mean_field_normalized_prior_at_its_point_times_is_the_kingman_density(
    #[case] tree: &str,
    #[case] times: &[(&str, f64)],
    #[case] kingman_neg_log_density: f64,
  ) -> Result<(), Report> {
    let parsed = nwk_read(tree.as_bytes())?;
    let names = parsed.names.clone();
    let graph = parsed.graph;
    let prior = point_coalescent_prior(&graph, &names, times, 3.0)?;

    let mut factors = 0.0;
    for &(name, time) in times {
      let node = graph
        .get_node(find_node_key_by_name(&graph, &names, name).expect("node must exist"))
        .expect("node must exist");
      factors += prior.node_contribution(node.key(), node.outbound().len(), time)?;
    }

    pretty_assert_abs_diff_eq!(kingman_neg_log_density, factors - prior.log_normalization(), epsilon = 1e-9);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::root(   "root", 2, (-0.5, 0.5))]
  #[case::leaf(   "A",    0, ( 0.0, 0.0))]
  #[trace]
  fn test_mean_field_outward_slopes_continue_the_contribution_past_the_lineage_grid(
    #[case] name: &str,
    #[case] n_children: usize,
    #[case] (left, right): (f64, f64),
  ) -> Result<(), Report> {
    let tc = 2.0;
    let parsed = nwk_read(b"(A:1,B:1)root;".as_slice())?;
    let names = parsed.names.clone();
    let graph = parsed.graph;
    let prior = point_coalescent_prior(&graph, &names, &[("A", 2010.0), ("B", 2010.0), ("root", 2000.0)], tc)?;
    let key = find_node_key_by_name(&graph, &names, name).expect("node must exist");
    let secant = |a: f64, b: f64| -> Result<f64, Report> {
      Ok((prior.node_contribution(key, n_children, b)? - prior.node_contribution(key, n_children, a)?) / (b - a))
    };

    let actual = (
      prior.node_contribution_outward_slope(key, n_children, Side::Left, 1950.0),
      prior.node_contribution_outward_slope(key, n_children, Side::Right, 2050.0),
    );

    assert_eq!((left, right), actual);
    pretty_assert_abs_diff_eq!(secant(1900.0, 1950.0)?, actual.0, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(secant(2050.0, 2100.0)?, actual.1, epsilon = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::left_inside_the_grid_before_the_merger( Side::Left,  2010.0, 1.0 - 2.0 - 1.0)]
  #[case::left_at_the_first_grid_point(           Side::Left,  2005.0, 1.0 - 1.0 - 1.0)]
  #[case::left_beyond_the_grid(                   Side::Left,  2000.0, 1.0 - 1.0 - 1.0)]
  #[case::right_after_the_merger(                 Side::Right, 2012.0, 1.0 - 3.0      )]
  #[case::right_after_the_first_tip(              Side::Right, 2013.5, 1.0 - 2.0      )]
  #[case::right_beyond_the_grid(                  Side::Right, 2020.0, 1.0 - 0.0      )]
  #[trace]
  fn test_mean_field_outward_slope_is_taken_just_outside_the_edge(
    #[case] side: Side,
    #[case] edge: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let parsed = nwk_read(b"((A:1,B:1)I:1,C:1)root;".as_slice())?;
    let names = parsed.names.clone();
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

    let actual = prior.node_contribution_outward_slope(internal, 2, side, edge);

    pretty_assert_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_mean_field_internal_slope_counts_the_lineages_without_the_own_merger() -> Result<(), Report> {
    let parsed = nwk_read(b"((A:1,B:1)I:1,C:1)root;".as_slice())?;
    let names = parsed.names.clone();
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
    let names = parsed.names.clone();
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
    let names = parsed.names.clone();
    let graph = parsed.graph;
    let flat_root = Distribution::function(Array1::linspace(2000.0, 2004.0, 401), Array1::zeros(401))?;
    let expectation = helpers::expectation(&graph, &names, &[("A", 2010.0), ("B", 2010.0)], &[("root", flat_root)])?;

    let SegmentStatistics { pair_integral, mergers } = expectation.segment_statistics(&[2000.0, 2010.0]);

    pretty_assert_abs_diff_eq!(8.0, pair_integral[0], epsilon = 1e-11);
    pretty_assert_abs_diff_eq!(1.0, mergers[0], epsilon = 1e-9);
    Ok(())
  }

  #[test]
  fn test_mean_field_segment_statistics_read_the_survival_inside_a_long_earlier_tail() -> Result<(), Report> {
    let parsed = nwk_read(b"(A:1,B:1)root;".as_slice())?;
    let names = parsed.names.clone();
    let graph = parsed.graph;
    let unit_box_with_tail_mass_one = Distribution::Function(
      DistributionFunction::<f64, NegLog>::from_range_values((-1.0, 0.0), array![0.0, 0.0])?
        .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -1.0 }))?
        .with_right_extrap(BoundaryBehavior::Hard)?,
    );
    let expectation = helpers::expectation(
      &graph,
      &names,
      &[("A", 1.0), ("B", 1.0)],
      &[("root", unit_box_with_tail_mass_one)],
    )?;

    let SegmentStatistics { mergers, .. } = expectation.segment_statistics(&[-30.0, -2.0, 1.0]);

    let mass_before_minus_2 = (-1.0_f64).exp() / 2.0;
    pretty_assert_abs_diff_eq!(mass_before_minus_2, mergers[0], epsilon = 1e-3);
    pretty_assert_abs_diff_eq!(1.0 - mass_before_minus_2, mergers[1], epsilon = 1e-3);
    Ok(())
  }

  #[test]
  fn test_mean_field_segment_statistics_read_the_survival_of_a_narrow_posterior_with_a_long_soft_tail()
  -> Result<(), Report> {
    let parsed = nwk_read(b"(A:1,B:1)root;".as_slice())?;
    let names = parsed.names.clone();
    let graph = parsed.graph;
    let tent = Array1::linspace(0.0, 1.0, 11).mapv(|t: f64| 1.0 - 1.8 * (t - 0.5).abs());
    let rate = 0.2;
    let root = Distribution::Function(
      DistributionFunction::<f64, NegLog>::from_range_values((0.0, 1.0), tent.mapv(|p| -p.ln()))?
        .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -rate }))?
        .with_right_extrap(BoundaryBehavior::Hard)?,
    );
    let expectation = helpers::expectation(&graph, &names, &[("A", 2.0), ("B", 2.0)], &[("root", root)])?;
    let boundaries = Array1::linspace(-110.0, 1.0, 2221).to_vec();

    let SegmentStatistics { mergers, .. } = expectation.segment_statistics(&boundaries);

    let tail_mass = 0.1 / rate;
    let total = tail_mass + 0.55;
    let mass_below = |t: f64| {
      if t <= 0.0 {
        tail_mass * (rate * t).exp()
      } else if t <= 0.5 {
        tail_mass + 0.1 * t + 0.9 * t * t
      } else {
        tail_mass + 0.275 + (t - 0.5) - 0.9 * (t - 0.5).powi(2)
      }
    };
    let read_survival = mergers.iter().scan(1.0, |survival, merger| {
      *survival -= merger;
      Some(*survival)
    });
    for (&t, survival) in boundaries[1..].iter().zip(read_survival) {
      pretty_assert_abs_diff_eq!(1.0 - mass_below(t) / total, survival, epsilon = 1e-2);
    }
    Ok(())
  }

  #[test]
  fn test_mean_field_own_integral_follows_a_tc_step_inside_a_knot_interval() -> Result<(), Report> {
    let parsed = nwk_read(b"(A:1,B:1)root;".as_slice())?;
    let names = parsed.names.clone();
    let graph = parsed.graph;
    let density = [1.0, 3.0, 2.0, 4.0, 1.0];
    let root = Distribution::Function(
      DistributionFunction::<f64, NegLog>::from_range_values(
        (0.0, 4.0),
        Array1::from_iter(density.map(|p: f64| -p.ln())),
      )?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Hard)?,
    );
    let step = 0.5;
    let expectation = helpers::expectation(&graph, &names, &[("A", step), ("B", 3.0)], &[("root", root)])?;
    let (tc_before, tc_after) = (2.0, 0.5);
    let tc_grid = Array1::linspace(-1.0, 5.0, 6 * 2_usize.pow(14) + 1);
    let tc = Distribution::function(
      tc_grid.clone(),
      tc_grid.mapv(|t| if t < step { tc_before } else { tc_after }),
    )?;
    let prior = CoalescentPrior::new(&expectation, &tc)?;
    let key = find_node_key_by_name(&graph, &names, "root").expect("root must exist");

    let survival_integral = |a: f64, b: f64| helpers::trapezoid_survival_integral(&density, a, b) / 10.0;
    let start = f64::midpoint(-1.0, (1.0 + 4e-8_f64).sqrt());
    let own =
      |x: f64| survival_integral(start, x.min(step)) / tc_before + survival_integral(step, x.max(step)) / tc_after;
    for x in [0.25, 0.75, 1.5, 2.5, 3.99] {
      pretty_assert_abs_diff_eq!(own(x), prior.own_integral_at(key, x), epsilon = 1e-3);
    }
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

    pretty_assert_abs_diff_eq!(1.0, mergers[1], epsilon = 1e-5);
    Ok(())
  }

  #[test]
  fn test_mean_field_gaussian_seed_counts_an_undated_tip_by_its_marginal() -> Result<(), Report> {
    let parsed = nwk_read(b"(A:1,B:1)root;".as_slice())?;
    let names = parsed.names.clone();
    let graph = parsed.graph;
    let key = |name| find_node_key_by_name(&graph, &names, name).expect("node must exist");
    let marginals = BTreeMap::from([
      (key("root"), NodeTimeGaussian { mean: 2000.0, std: 0.5 }),
      (key("B"), NodeTimeGaussian { mean: 2006.0, std: 0.2 }),
    ]);
    let tip_times = BTreeMap::from([(key("A"), Some(2010.0)), (key("B"), None), (key("root"), None)]);

    let expectation = LineageExpectation::from_gaussian(&graph, &marginals, &tip_times, &BTreeMap::new())?;

    let SegmentStatistics { pair_integral, .. } =
      expectation.segment_statistics(&[1990.0, 2002.0, 2005.0, 2008.0, 2010.0]);
    pretty_assert_abs_diff_eq!(2.0, pair_integral[0], epsilon = 1e-5);
    pretty_assert_abs_diff_eq!(3.0, pair_integral[1], epsilon = 1e-5);
    pretty_assert_abs_diff_eq!(0.0, pair_integral[3], epsilon = 1e-5);
    Ok(())
  }

  #[test]
  fn test_mean_field_own_share_is_removed_at_its_weight_in_the_expectation() -> Result<(), Report> {
    let parsed = nwk_read(b"(A:1,B:1,C:1)root;".as_slice())?;
    let names = parsed.names.clone();
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
    let names = parsed.names.clone();
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
    let tc = Distribution::constant(3.0);

    let extended = LineageExpectation::from_posterior(&graph, &without_n)?
      .with_topology(&graph, &BTreeMap::from([(key("N"), 2005.0)]))?;
    let fresh = LineageExpectation::from_posterior(&graph, &full)?;
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
  fn test_mean_field_of_inference_counts_undated_tips_by_their_posterior() -> Result<(), Report> {
    let parsed = nwk_read(b"(A:1,B:1,C:1,D:1)root;".as_slice())?;
    let names = parsed.names.clone();
    let graph = parsed.graph;
    let key = |name| find_node_key_by_name(&graph, &names, name).expect("node must exist");
    let dated = NodePosterior {
      time: Some(2015.0),
      ..NodePosterior::default()
    };
    let flat = |from: f64, to: f64| -> Result<NodePosterior, Report> {
      let distribution = Distribution::function(Array1::linspace(from, to, 201), Array1::zeros(201))?;
      Ok(NodePosterior {
        survival: NodeSurvival::of_posterior(&distribution)?.map(Arc::new),
        ..NodePosterior::default()
      })
    };
    let inference = TimeInference {
      bad_branches: BTreeMap::from([
        (key("A"), false),
        (key("B"), false),
        (key("C"), true),
        (key("D"), true),
        (key("root"), false),
      ]),
      posterior: BTreeMap::from([
        (key("A"), dated.clone()),
        (key("B"), dated),
        (key("C"), flat(2006.0, 2008.0)?),
        (key("D"), flat(2006.0, 2008.0)?),
        (key("root"), flat(2000.0, 2004.0)?),
      ]),
      posterior_summaries: PosteriorSummaries {
        regions: false,
        lineage_survival: true,
      },
      ..empty_time_inference(&graph)
    };

    let expectation = LineageExpectation::of_inference(&graph, &inference)?;

    let SegmentStatistics { pair_integral, .. } =
      expectation.segment_statistics(&[2000.0, 2004.0, 2006.0, 2008.0, 2015.0]);
    pretty_assert_abs_diff_eq!(12.0, pair_integral[0], epsilon = 1e-5);
    pretty_assert_abs_diff_eq!(12.0, pair_integral[1], epsilon = 1e-5);
    pretty_assert_abs_diff_eq!(7.0, pair_integral[3], epsilon = 1e-5);
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

  #[test]
  fn test_mean_field_prior_reports_a_tc_that_cannot_be_evaluated_at_the_lineage_grid_edge() -> Result<(), Report> {
    let parsed = nwk_read(b"(A:1,B:1)root;".as_slice())?;
    let names = parsed.names.clone();
    let graph = parsed.graph;
    let flat_root = Distribution::function(Array1::linspace(2000.0, 2004.0, 401), Array1::zeros(401))?;
    let expectation = helpers::expectation(&graph, &names, &[("A", 2010.0), ("B", 2010.0)], &[("root", flat_root)])?;
    let tc = Distribution::function(Array1::linspace(2000.001, 2011.0, 11), Array1::from_elem(11, 2.0))?;

    let result = CoalescentPrior::new(&expectation, &tc);

    assert_error!(
      result,
      "When evaluating coalescent Tc at calendar time 2.000000e3: GridFn evaluated at 2000.000000004, below the support boundary 2000.001, but no extrapolation policy is set for that side"
    );
    Ok(())
  }

  mod helpers {

    pub(super) fn gaussian_cherry(
      mean: f64,
      std: f64,
      excluded: &[&str],
    ) -> Result<(Graph, NodeNames, Arc<LineageExpectation>), Report> {
      gaussian_cherry_dated(mean, std, [2010.0, 2010.0], excluded)
    }

    pub(super) fn gaussian_cherry_dated(
      mean: f64,
      std: f64,
      dates: [f64; 2],
      excluded: &[&str],
    ) -> Result<(Graph, NodeNames, Arc<LineageExpectation>), Report> {
      let parsed = nwk_read(b"(A:1,B:1)root;".as_slice())?;
      let names = parsed.names.clone();
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

    pub(super) fn trapezoid_survival_integral(density: &[f64], a: f64, b: f64) -> f64 {
      let cells = density.len() - 1;
      let later = |cell: usize| {
        (cell + 1..cells)
          .map(|k| f64::midpoint(density[k], density[k + 1]))
          .sum::<f64>()
      };
      let primitive = |cell: usize, u: f64| {
        let (p0, p1) = (density[cell], density[cell + 1]);
        (f64::midpoint(p0, p1) + later(cell)) * u - p0 * u * u / 2.0 - (p1 - p0) * u * u * u / 6.0
      };
      (0..cells)
        .map(|cell| {
          let x0 = cell.to_f64_exact();
          let (lo, hi) = (a.max(x0), b.min(x0 + 1.0));
          if lo < hi {
            primitive(cell, hi - x0) - primitive(cell, lo - x0)
          } else {
            0.0
          }
        })
        .sum()
    }

    pub(super) fn expectation(
      graph: &Graph,
      names: &NodeNames,
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
      LineageExpectation::from_posterior(graph, &posterior)
    }
  }
}

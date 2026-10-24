#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::coalescent::mean_field::NodeSurvival;
  use crate::progress::{LogSink, NoopProgress};
  use crate::test_utils::{
    FullPosterior, RecordingLog, find_node_key_by_name, full_posterior, parent_edge_key, unknown_branches,
  };
  use crate::timetree::confidence::extract_confidence_intervals;
  use crate::timetree::inference::forward_pass::{ForwardPassInputs, committed_time, propagate_distributions_forward};
  use crate::timetree::inference::result::{
    BranchLikelihood, NodeBackward, NodePosterior, PosteriorSummaries, PosteriorView, TimeInference,
  };
  use crate::timetree::inference::runner::{EPS, GRID_POINTS, TRIM_EPS};
  use crate::timetree::max_posterior_region::MaxPosteriorRegion;
  use crate::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use ndarray::{Array1, array};
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use rayon::ThreadPoolBuilder;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use std::sync::Weak;
  use treetime_distribution::{
    Distribution, DistributionFunction, MassProfile, NegLog, TailMass, convolve_across_edge,
    distribution_multiplication,
  };
  use treetime_graph::key_map::NodeMap;
  use treetime_graph::keymap;

  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_graph::node_names::NodeNames;
  use treetime_grid::{BoundaryBehavior, MaxGridPoints, Side, SoftTailLaw};
  use treetime_io::nwk::nwk_read;
  use treetime_utils::assert_error;
  use treetime_utils::cast::ToF64Exact;

  const CAVITY_GRID_POINTS: usize = 1001;

  #[rustfmt::skip]
  #[rstest]
  #[case::no_parent(      None,      5.0)]
  #[case::parent_earlier( Some(3.0), 5.0)]
  #[case::parent_clamps(  Some(8.0), 8.0)]
  #[trace]
  fn test_forward_pass_committed_time_clamps_the_likely_time_to_the_parent(
    #[case] parent_time: Option<f64>,
    #[case] expected: f64,
  ) {
    let committed = committed_time(5.0, parent_time);
    pretty_assert_ulps_eq!(committed, expected, max_ulps = 4);
  }

  #[test]
  fn test_forward_pass_leaves_internal_node_with_empty_distribution_undated() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:2.5)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;

    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(&mut inputs, root_key, Distribution::empty());
    set_date(&mut inputs, leaf_key, Distribution::point(2013.0, 0.0));

    let log = RecordingLog::default();
    let posterior = run_forward_pass_with_log(&graph, &inputs, &names, &log)?;

    let root_time = node_time(&posterior, root_key);

    assert_eq!(None, root_time);
    let leaf_time = node_time(&posterior, leaf_key).expect("leaf A should keep its observed date");
    pretty_assert_ulps_eq!(leaf_time, 2013.0, max_ulps = 4);
    let expected_warnings = vec![
      "Timetree forward pass: node 'root' has an empty time distribution; no date was assigned. The messages \
       meeting at this node leave no time with any probability: the dates below it and the times the rest of the \
       tree implies have disjoint support."
        .to_owned(),
    ];
    assert_eq!(expected_warnings, log.warnings());

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_thread(  1)]
  #[case::four_threads(4)]
  #[trace]
  fn test_forward_pass_warns_about_undated_nodes_in_preorder(#[case] threads: usize) -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:1.0,B:1.0,C:1.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");

    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(&mut inputs, key("root"), Distribution::empty());
    for name in ["A", "B", "C"] {
      set_date(&mut inputs, key(name), Distribution::range((2010.0, 2011.0), 0.0));
      set_time_distribution(&mut inputs, key(name), Distribution::empty());
    }

    let log = RecordingLog::default();
    ThreadPoolBuilder::new()
      .num_threads(threads)
      .build()?
      .install(|| run_forward_pass_with_log(&graph, &inputs, &names, &log))?;

    let expected = ["root", "A", "B", "C"].map(|name| {
      format!(
        "Timetree forward pass: node '{name}' has an empty time distribution; no date was assigned. The messages \
         meeting at this node leave no time with any probability: the dates below it and the times the rest of the \
         tree implies have disjoint support."
      )
    });
    assert_eq!(expected.to_vec(), log.warnings());
    Ok(())
  }

  #[test]
  fn test_forward_pass_refines_uncertain_leaf_date_from_parent() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:1.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(&mut inputs, root_key, Distribution::point(2009.0, 0.0));
    set_date(&mut inputs, leaf_key, Distribution::range((2009.5, 2013.5), 0.0));
    set_branch_length_distribution(&graph, &mut inputs, leaf_key, 1.0);

    let posterior = run_forward_pass(&graph, &inputs, &names)?;

    let leaf_time = node_time(&posterior, leaf_key).expect("leaf A should be dated");
    pretty_assert_ulps_eq!(leaf_time, 2010.0, max_ulps = 4);

    let leaf_dist = leaf_time_distribution(&posterior, leaf_key).expect("leaf A should have a distribution");
    assert_eq!(Distribution::point(2010.0, 0.0), leaf_dist);

    Ok(())
  }

  #[test]
  fn test_forward_pass_keeps_exact_leaf_date_earlier_than_its_parent() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:1.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(&mut inputs, root_key, Distribution::point(2009.0, 0.0));
    set_date(&mut inputs, leaf_key, Distribution::point(2008.0, 0.0));
    set_branch_length_distribution(&graph, &mut inputs, leaf_key, 1.0);

    let posterior = run_forward_pass(&graph, &inputs, &names)?;

    let leaf_time = node_time(&posterior, leaf_key).expect("leaf A should keep its observed date");
    pretty_assert_ulps_eq!(leaf_time, 2008.0, max_ulps = 4);

    let leaf_dist = leaf_time_distribution(&posterior, leaf_key).expect("leaf A should have a distribution");
    assert_eq!(Distribution::point(2008.0, 0.0), leaf_dist);

    Ok(())
  }

  #[test]
  fn test_forward_pass_clamps_uncertain_leaf_date_to_parent_time() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:1.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(&mut inputs, root_key, Distribution::point(2009.0, 0.0));
    set_date(&mut inputs, leaf_key, Distribution::range((2005.0, 2007.0), 0.0));

    let posterior = run_forward_pass(&graph, &inputs, &names)?;

    let leaf_time = node_time(&posterior, leaf_key).expect("leaf A should be dated");
    pretty_assert_ulps_eq!(leaf_time, 2009.0, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_forward_pass_keeps_uncertain_leaf_date_the_tree_contradicts() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:1.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(&mut inputs, root_key, Distribution::point(2009.0, 0.0));
    let given = Distribution::range((2005.0, 2007.0), 0.0);
    set_date(&mut inputs, leaf_key, given.clone());
    set_branch_length_distribution(&graph, &mut inputs, leaf_key, 1.0);

    let log = RecordingLog::default();
    let posterior = run_forward_pass_with_log(&graph, &inputs, &names, &log)?;

    let expected = FullPosterior {
      distribution: Some(Arc::new(given)),
      likely_time: Some(2006.0),
      time: Some(2009.0),
      contradicted: true,
      subtree_support: Some((2005.0, 2007.0)),
    };
    assert_eq!(expected, posterior[&leaf_key]);
    let expected_warnings = vec![
      "Timetree forward pass: 1 node(s) carry a date that the rest of the tree gives no probability at all, so \
       their posterior came out empty and each kept the date it was given, unrefined. The usual cause is a \
       sequence whose divergence implies a date far from the one it is stamped with, which the clock filter \
       reports separately. Run with `--verbosity=debug` to see which nodes and where the tree puts each of them."
        .to_owned(),
    ];
    assert_eq!(expected_warnings, log.warnings());

    Ok(())
  }

  #[test]
  fn test_forward_pass_refinement_is_unchanged_by_the_reachable_window() -> Result<(), Report> {
    let wide_parent = {
      let t = Array1::linspace(1950.0, 2050.0, 4001);
      let y = t.mapv(|t: f64| 0.5 * ((t - 2009.0) / 3.0).powi(2));
      Distribution::function(t, y)?
    };

    let refine = |parent: Distribution<NegLog>| -> Result<f64, Report> {
      let nwk_parsed = nwk_read(b"(A:1.0)root;".as_slice())?;
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
      let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
      let mut inputs = ForwardInputs::new(&graph);
      set_time_distribution(&mut inputs, root_key, parent);
      set_date(&mut inputs, leaf_key, Distribution::range((2010.5, 2010.6), 0.0));
      set_branch_length_distribution(&graph, &mut inputs, leaf_key, 1.0);
      let posterior = run_forward_pass(&graph, &inputs, &names)?;
      node_time(&posterior, leaf_key).ok_or_else(|| eyre::eyre!("leaf A should be dated"))
    };

    let Distribution::Function(wide_function) = &wide_parent else {
      return Err(eyre::eyre!("the wide parent must be a distribution function"));
    };
    let windowed = Distribution::Function(wide_function.resample_range_dx(
      (2009.4, 2009.7),
      wide_function.dx(),
      MaxGridPoints::default(),
    )?);

    pretty_assert_abs_diff_eq!(refine(wide_parent)?, refine(windowed)?, epsilon = 1e-9);

    Ok(())
  }

  #[test]
  fn test_forward_pass_refines_a_date_range_narrower_than_the_parent_grid() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:1.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let coarse_parent = {
      let t = Array1::linspace(1900.0, 2020.0, 121);
      let y = t.mapv(|t: f64| 0.5 * ((t - 2009.0) / 5.0).powi(2));
      Distribution::function(t, y)?
    };
    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(&mut inputs, root_key, coarse_parent);
    set_date(&mut inputs, leaf_key, Distribution::range((2010.500, 2010.503), 0.0));
    set_branch_length_distribution(&graph, &mut inputs, leaf_key, 1.0);

    let posterior = run_forward_pass(&graph, &inputs, &names)?;

    let leaf_time = node_time(&posterior, leaf_key).expect("leaf A should be dated");
    assert!(
      (2010.500..=2010.503).contains(&leaf_time),
      "leaf A should be dated within the day it was given, got {leaf_time}"
    );

    Ok(())
  }

  #[test]
  fn test_forward_pass_divides_the_childs_own_message_out_of_the_parent_posterior() -> Result<(), Report> {
    let grid = Array1::linspace(2000.0, 2010.0, CAVITY_GRID_POINTS);
    let parent_posterior = Distribution::function(grid.clone(), quadratic_neglog(&grid, 2005.0, 4.0))?;
    let message_from_child = Distribution::function(grid.clone(), quadratic_neglog(&grid, 2004.0, 2.0))?;
    let analytic_cavity = Distribution::function(grid.clone(), quadratic_neglog(&grid, 2006.0, 2.0))?;

    let with_message = refine_child_of(parent_posterior.clone(), Some(message_from_child))?;
    let from_analytic_cavity = refine_child_of(analytic_cavity, None)?;
    let without_division = refine_child_of(parent_posterior, None)?;

    assert_eq!(from_analytic_cavity, with_message);
    let time = with_message.expect("the child must be dated");
    assert!(
      (time - 2007.0).abs() < (time - 2006.0).abs(),
      "the cavity peaks at 2006 and the branch adds 1, so the child must sit near 2007, got {time}"
    );
    assert_ne!(without_division, with_message);
    Ok(())
  }

  #[test]
  fn test_forward_pass_trims_each_refined_posterior_to_its_mass() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:1.0,B:1.0)N:2.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("node not found");
    let (root_key, node_key, a_key, b_key) = (key("root"), key("N"), key("A"), key("B"));

    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(
      &mut inputs,
      root_key,
      hard_gaussian((1000.0, 2010.0), 20_001, 2000.0, 2.0)?,
    );
    set_time_distribution(
      &mut inputs,
      node_key,
      hard_gaussian((1000.0, 2010.0), 20_001, 2003.0, 1.5)?,
    );
    set_date(&mut inputs, a_key, Distribution::range((2004.5, 2005.0), 0.0));
    set_date(&mut inputs, b_key, Distribution::range((2005.0, 2005.5), 0.0));
    set_branch_distribution(
      &graph,
      &mut inputs,
      node_key,
      hard_gaussian((0.0, 10.0), 1001, 2.0, 0.5)?,
    );
    set_branch_distribution(&graph, &mut inputs, a_key, hard_gaussian((0.0, 10.0), 1001, 1.5, 0.5)?);
    set_branch_distribution(&graph, &mut inputs, b_key, hard_gaussian((0.0, 10.0), 1001, 1.5, 0.5)?);

    let posterior = run_forward_pass(&graph, &inputs, &names)?;

    for (parent, child) in [(root_key, node_key), (node_key, a_key), (node_key, b_key)] {
      let parent_time = posterior[&parent].likely_time.expect("every node is dated");
      let child_time = posterior[&child].likely_time.expect("every node is dated");
      assert!(
        child_time >= parent_time,
        "child {child} at {child_time} before parent {parent} at {parent_time}"
      );

      let parent_dist = posterior[&parent]
        .distribution
        .as_deref()
        .expect("the parent has a posterior");
      let branch = inputs.branches[&parent_edge_key(&graph, child)].distribution.as_deref();
      let branch = branch.expect("the branch has a distribution");
      let subtree = inputs.backward[&child]
        .subtree
        .as_deref()
        .expect("the child has a subtree distribution");
      let from_parent = convolve_across_edge(
        parent_dist,
        branch,
        Side::Right,
        EPS,
        GRID_POINTS,
        MaxGridPoints::default(),
      )?;
      let untrimmed = distribution_multiplication(&from_parent, subtree, MaxGridPoints::default())?;
      let Distribution::Function(untrimmed) = untrimmed else {
        return Err(eyre::eyre!("the untrimmed posterior of {child} must be a function"));
      };
      let profile = MassProfile::new(&untrimmed)?;
      let fractions = profile
        .fractions(TailMass::Required)
        .expect("the untrimmed posterior has a finite mass");
      let (lo, hi) = posterior[&child]
        .distribution
        .as_deref()
        .and_then(Distribution::time_bounds)
        .expect("dated");
      let inside = fractions.mass_between(lo, hi)? / fractions.total();
      assert!(
        inside >= 1.0 - 2.0 * TRIM_EPS,
        "the trimmed posterior of {child} holds {inside} of the untrimmed mass"
      );
      if child == node_key {
        assert!(
          inside < 1.0,
          "the internal node N reaches back to the root's hard edge and must be trimmed"
        );
      }
    }
    Ok(())
  }

  #[test]
  fn test_forward_pass_child_interval_is_not_clipped_at_the_cut_of_its_parent() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(((A:1.0)C:1.0)P:2.0)R;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("node not found");
    let (parent_key, child_key) = (key("P"), key("C"));

    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(
      &mut inputs,
      key("R"),
      hard_gaussian((1000.0, 2010.0), 20_001, 2000.0, 2.0)?,
    );
    set_time_distribution(
      &mut inputs,
      parent_key,
      hard_gaussian((1000.0, 2010.0), 20_001, 2003.0, 1.5)?,
    );
    set_time_distribution(
      &mut inputs,
      child_key,
      soft_left_gaussian((1990.0, 2010.0), 4001, 1995.0, 1.0)?,
    );
    set_branch_distribution(
      &graph,
      &mut inputs,
      parent_key,
      hard_gaussian((0.0, 10.0), 1001, 2.0, 0.5)?,
    );
    set_branch_length_distribution(&graph, &mut inputs, child_key, 1.0);

    let posterior = run_forward_pass(&graph, &inputs, &names)?;

    let (parent_lower, _) = posterior[&parent_key]
      .distribution
      .as_deref()
      .and_then(Distribution::time_bounds)
      .expect("P has a posterior");
    let (child_lower, _) = posterior[&child_key]
      .distribution
      .as_deref()
      .and_then(Distribution::time_bounds)
      .expect("C has a posterior");
    pretty_assert_abs_diff_eq!(parent_lower + 1.0, child_lower, epsilon = 1e-9);
    let subtree = inputs.backward[&child_key].subtree.as_deref();
    let (_, subtree_upper) = subtree
      .and_then(Distribution::time_bounds)
      .expect("C has a subtree distribution");
    assert_eq!(
      Some((f64::NEG_INFINITY, subtree_upper)),
      posterior[&child_key].subtree_support
    );

    let summaries = PosteriorSummaries {
      regions: true,
      ..PosteriorSummaries::default()
    };
    let inference = TimeInference {
      posterior: run_forward_pass_summaries(&graph, &inputs, &names, summaries)?,
      posterior_summaries: summaries,
      ..TimeInference::default()
    };
    let date = inference.posterior[&child_key].time.expect("C is dated");
    let rate_dates = keymap! { child_key => [date - 50.0, date, date + 50.0] };
    let intervals = extract_confidence_intervals(&graph, &inference, &rate_dates, &names)?;
    let unclipped = {
      let mut inference = inference;
      inference
        .posterior
        .get_mut(child_key)
        .expect("C has a posterior")
        .subtree_support = None;
      extract_confidence_intervals(&graph, &inference, &rate_dates, &names)?
    };
    let child = intervals
      .iter()
      .find(|interval| interval.key == child_key)
      .expect("C has an interval");
    let child_unclipped = unclipped
      .iter()
      .find(|interval| interval.key == child_key)
      .expect("C has an interval");

    assert!(
      child.lower < child_lower - 10.0,
      "the lower bound {} reaches past the parent's cut",
      child.lower
    );
    pretty_assert_ulps_eq!(child_unclipped.lower, child.lower, max_ulps = 0);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::posterior_as_wide_as_the_message( 4001, 3.0, false)]
  #[case::posterior_narrower_than_the_cap(  4001, 0.1, true )]
  #[case::uncapped_spacing_bounds_the_count( 401, 0.1, false)]
  #[trace]
  fn test_forward_pass_message_resolves_the_posterior_despite_the_spacing_cap(
    #[case] root_points: usize,
    #[case] subtree_sigma: f64,
    #[case] refined: bool,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(N:2.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("node not found");
    let (root_key, node_key) = (key("root"), key("N"));

    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(&mut inputs, root_key, hard_gaussian((1998.0, 2002.0), root_points, 2000.0, 0.5)?);
    set_time_distribution(&mut inputs, node_key, hard_gaussian((1990.0, 2020.0), 601, 2002.0, subtree_sigma)?);
    set_branch_distribution(&graph, &mut inputs, node_key, hard_gaussian((0.0, 4.0), 401, 2.0, 0.5)?);

    let posterior = run_forward_pass(&graph, &inputs, &names)?;

    let root = posterior[&root_key].distribution.as_deref().expect("the root has a posterior");
    let branch = inputs.branches[&parent_edge_key(&graph, node_key)].distribution.as_deref();
    let branch = branch.expect("the branch has a distribution");
    let capped = convolve_across_edge(root, branch, Side::Right, EPS, GRID_POINTS, MaxGridPoints::default())?;
    let (Distribution::Function(capped), Distribution::Function(root), Distribution::Function(branch)) =
      (capped, root, branch)
    else {
      return Err(eyre::eyre!("the capped message and its operands must be functions"));
    };
    let fine_dx = root.dx().min(branch.dx());
    let intervals = (GRID_POINTS - 1).to_f64_exact();
    let uncapped_dx = ((capped.x_max() - capped.x_min()) / intervals).min(fine_dx);
    let Some(Distribution::Function(result)) = posterior[&node_key].distribution.as_deref() else {
      return Err(eyre::eyre!("N must have a function posterior"));
    };
    let cells = ((result.x_max() - result.x_min()) / result.dx()).round();
    let required = ((result.x_max() - result.x_min()) / uncapped_dx).floor().min(intervals);

    pretty_assert_abs_diff_eq!(0.01, capped.dx(), epsilon = 1e-12);
    assert!(cells >= required, "the posterior spans {cells} cells, the uncapped path {required}");
    if refined {
      assert!(result.dx() < capped.dx(), "spacing {} not below the cap", result.dx());
      assert!(result.dx() >= fine_dx, "spacing {} below the convolution spacing", result.dx());
    } else {
      pretty_assert_abs_diff_eq!(capped.dx(), result.dx(), epsilon = 1e-12);
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::trimmed_product_with_the_subtree(Some(Distribution::range((-0.01, 2.01), 0.0)))]
  #[case::message_without_a_subtree(       None                                         )]
  #[trace]
  fn test_forward_pass_message_resolves_the_posterior_it_forms_after_a_resample(
    #[case] subtree: Option<Distribution<NegLog>>,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(N:1.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("node not found");
    let (root_key, node_key) = (key("root"), key("N"));

    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(&mut inputs, root_key, hard_gaussian((-0.01, 0.01), 2001, 0.0, 0.001)?);
    if let Some(subtree) = subtree {
      set_time_distribution(&mut inputs, node_key, subtree);
    }
    let branch = DistributionFunction::from_start_dx_values(0.0, 1.0, array![5e7, 0.0, 5e7])?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Hard)?;
    set_branch_distribution(&graph, &mut inputs, node_key, Distribution::Function(branch));

    let posterior = run_forward_pass(&graph, &inputs, &names)?;

    let Some(Distribution::Function(result)) = posterior[&node_key].distribution.as_deref() else {
      return Err(eyre::eyre!("N must have a function posterior"));
    };
    let profile = MassProfile::new(result)?;
    let fractions = profile.fractions(TailMass::Required).expect("the posterior has a finite mass");
    let mass_width = fractions.survival_inverse(TRIM_EPS) - fractions.quantile(TRIM_EPS);
    let cells = mass_width / result.dx();
    assert!(
      cells >= (GRID_POINTS - 1).to_f64_exact(),
      "the message resolves the posterior mass window {mass_width} with {cells} cells of {}",
      result.dx()
    );
    Ok(())
  }

  #[test]
  fn test_forward_pass_missing_backward_output_is_an_internal_error() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:1.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut inputs = ForwardInputs::new(&graph);
    inputs.backward.remove(leaf_key);

    assert_error!(
      run_forward_pass(&graph, &inputs, &names),
      format!("Scan input is missing node {leaf_key}. This is an internal error. Please report it to developers.")
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::not_requested(false, None)]
  #[case::requested(    true,  Some(MaxPosteriorRegion::Point))]
  #[trace]
  fn test_forward_pass_stores_the_posterior_region_only_when_requested(
    #[case] regions: bool,
    #[case] expected: Option<MaxPosteriorRegion>,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:1.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(&mut inputs, root_key, Distribution::point(2009.0, 0.0));
    set_date(&mut inputs, leaf_key, Distribution::point(2008.0, 0.0));
    set_branch_length_distribution(&graph, &mut inputs, leaf_key, 1.0);
    let summaries = PosteriorSummaries {
      regions,
      ..PosteriorSummaries::default()
    };

    let posterior = run_forward_pass_summaries(&graph, &inputs, &names, summaries)?;

    assert_eq!(expected, posterior[&leaf_key].region);
    Ok(())
  }

  #[test]
  fn test_forward_pass_stores_no_region_for_an_undated_node() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:1.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(&mut inputs, root_key, Distribution::point(2009.0, 0.0));
    let summaries = PosteriorSummaries {
      regions: true,
      ..PosteriorSummaries::default()
    };

    let posterior = run_forward_pass_summaries(&graph, &inputs, &names, summaries)?;

    assert_eq!(NodePosterior::default(), posterior[&leaf_key]);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::not_requested(false, false)]
  #[case::requested(    true,  true )]
  #[trace]
  fn test_forward_pass_summary_stores_the_survival_only_when_requested(
    #[case] lineage_survival: bool,
    #[case] stored: bool,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:1.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let leaf_key = find_node_key_by_name(&nwk_parsed.graph, &names, "A").expect("leaf A not found");
    let distribution = Arc::new(Distribution::range((2004.0, 2006.0), 0.0));
    let view = PosteriorView {
      distribution: Some(&distribution),
      likely_time: Some(2005.0),
      time: Some(2005.0),
      contradicted: false,
      subtree_support: None,
    };
    let summaries = PosteriorSummaries {
      lineage_survival,
      ..PosteriorSummaries::default()
    };

    let posterior = NodePosterior::summarize(leaf_key, view, summaries)?;

    let expected = if stored {
      NodeSurvival::of_posterior(&distribution)?.map(Arc::new)
    } else {
      None
    };
    assert_eq!(expected, posterior.survival);
    Ok(())
  }

  #[test]
  fn test_forward_pass_frees_every_posterior_distribution_when_it_returns() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:1.0,B:1.0)I:1.0,C:1.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");

    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(&mut inputs, key("root"), Distribution::range((2000.0, 2002.0), 0.0));
    for (name, date) in [("A", 2010.0), ("B", 2009.0), ("C", 2011.0)] {
      set_date(&mut inputs, key(name), Distribution::point(date, 0.0));
    }
    for name in ["A", "B", "C", "I"] {
      set_branch_length_distribution(&graph, &mut inputs, key(name), 1.0);
    }

    let posterior = propagate_distributions_forward(
      &graph,
      &forward_inputs(&inputs, &names),
      inputs.backward.clone(),
      |_, view| Ok(view.distribution.map(Arc::downgrade)),
      &NoopProgress,
    )?;

    let live = posterior
      .into_values()
      .map(|weak| weak.as_ref().map(Weak::strong_count))
      .collect_vec();
    assert_eq!(vec![Some(0); 5], live);
    Ok(())
  }

  #[test]
  fn test_forward_pass_frees_the_grids_of_a_finished_sibling_subtree_before_the_next_sibling() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:1.0,B:1.0)I:1.0,C:1.0)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");

    let mut inputs = ForwardInputs::new(&graph);
    set_time_distribution(&mut inputs, key("root"), Distribution::range((2000.0, 2002.0), 0.0));
    for (name, date) in [("A", 2010.0), ("B", 2009.0), ("C", 2011.0)] {
      set_date(&mut inputs, key(name), Distribution::point(date, 0.0));
    }
    for name in ["A", "B", "C", "I"] {
      set_branch_length_distribution(&graph, &mut inputs, key(name), 1.0);
    }
    let grids: Mutex<NodeMap<Weak<Distribution<NegLog>>>> = Mutex::new(NodeMap::new());

    let live_before_each_node = ThreadPoolBuilder::new().num_threads(1).build()?.install(|| {
      propagate_distributions_forward(
        &graph,
        &forward_inputs(&inputs, &names),
        inputs.backward.clone(),
        |key, view| {
          let mut grids = grids.lock();
          let live = grids
            .iter()
            .map(|(key, grid)| (key, grid.strong_count()))
            .collect::<BTreeMap<_, _>>();
          if let Some(distribution) = view.distribution {
            grids.insert(key, Arc::downgrade(distribution));
          }
          Ok(live)
        },
        &NoopProgress,
      )
    })?;

    let expected = btreemap! { key("root") => 1, key("I") => 0, key("A") => 0, key("B") => 0 };
    assert_eq!(expected, live_before_each_node[&key("C")]);
    Ok(())
  }

  mod helpers {
    use super::*;
    use treetime_graph::key_map::{EdgeMap, NodeMap};

    pub(super) struct ForwardInputs {
      pub constraints: DateConstraints,
      pub bad_branches: NodeMap<bool>,
      pub branches: EdgeMap<BranchLikelihood>,
      pub backward: NodeMap<NodeBackward>,
    }

    impl ForwardInputs {
      pub(super) fn new(graph: &Graph) -> Self {
        Self {
          constraints: DateConstraints::default(),
          bad_branches: graph.get_nodes().map(|node| (node.key(), false)).collect(),
          branches: unknown_branches(graph),
          backward: graph
            .get_nodes()
            .map(|node| (node.key(), NodeBackward::default()))
            .collect(),
        }
      }
    }

    pub(super) fn run_forward_pass(
      graph: &Graph,
      inputs: &ForwardInputs,
      names: &NodeNames,
    ) -> Result<NodeMap<FullPosterior>, Report> {
      run_forward_pass_with_log(graph, inputs, names, &NoopProgress)
    }

    pub(super) fn run_forward_pass_with_log(
      graph: &Graph,
      inputs: &ForwardInputs,
      names: &NodeNames,
      log: &dyn LogSink,
    ) -> Result<NodeMap<FullPosterior>, Report> {
      propagate_distributions_forward(
        graph,
        &forward_inputs(inputs, names),
        inputs.backward.clone(),
        full_posterior,
        log,
      )
    }

    pub(super) fn run_forward_pass_summaries(
      graph: &Graph,
      inputs: &ForwardInputs,
      names: &NodeNames,
      summaries: PosteriorSummaries,
    ) -> Result<NodeMap<NodePosterior>, Report> {
      propagate_distributions_forward(
        graph,
        &forward_inputs(inputs, names),
        inputs.backward.clone(),
        |key, view| NodePosterior::summarize(key, view, summaries),
        &NoopProgress,
      )
    }

    pub(super) fn forward_inputs<'a>(inputs: &'a ForwardInputs, names: &'a NodeNames) -> ForwardPassInputs<'a> {
      ForwardPassInputs {
        constraints: &inputs.constraints,
        coalescent: None,
        bad_branches: &inputs.bad_branches,
        names,
        branches: &inputs.branches,
        max_points: MaxGridPoints::default(),
      }
    }

    pub(super) fn set_date(inputs: &mut ForwardInputs, key: GraphNodeKey, dist: Distribution<NegLog>) {
      let dist = Arc::new(dist);
      inputs.constraints.by_node.insert(key, Some(Arc::clone(&dist)));
      inputs.backward.entry(key).or_default().subtree = Some(dist);
    }

    pub(super) fn set_time_distribution(inputs: &mut ForwardInputs, key: GraphNodeKey, dist: Distribution<NegLog>) {
      inputs.backward.entry(key).or_default().subtree = Some(Arc::new(dist));
    }

    pub(super) fn set_branch_length_distribution(
      graph: &Graph,
      inputs: &mut ForwardInputs,
      target_key: GraphNodeKey,
      branch_length: f64,
    ) {
      let branch = BranchLikelihood {
        distribution: Some(Arc::new(Distribution::point(branch_length, 0.0))),
        time_length: Some(branch_length),
      };
      inputs.branches.insert(parent_edge_key(graph, target_key), branch);
    }

    pub(super) fn set_branch_distribution(
      graph: &Graph,
      inputs: &mut ForwardInputs,
      target_key: GraphNodeKey,
      dist: Distribution<NegLog>,
    ) {
      let time_length = dist.likely_time().expect("the branch distribution has a peak");
      let branch = BranchLikelihood {
        distribution: Some(Arc::new(dist)),
        time_length,
      };
      inputs.branches.insert(parent_edge_key(graph, target_key), branch);
    }

    pub(super) fn hard_gaussian(
      (x_min, x_max): (f64, f64),
      n_points: usize,
      mean: f64,
      sigma: f64,
    ) -> Result<Distribution<NegLog>, Report> {
      let grid = Array1::linspace(x_min, x_max, n_points);
      let y = quadratic_neglog(&grid, mean, 1.0 / (sigma * sigma));
      let f = DistributionFunction::from_start_dx_values(x_min, grid[1] - grid[0], y)?
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(BoundaryBehavior::Hard)?;
      Ok(Distribution::Function(f))
    }

    pub(super) fn soft_left_gaussian(
      (x_min, x_max): (f64, f64),
      n_points: usize,
      mean: f64,
      sigma: f64,
    ) -> Result<Distribution<NegLog>, Report> {
      let grid = Array1::linspace(x_min, x_max, n_points);
      let y = quadratic_neglog(&grid, mean, 1.0 / (sigma * sigma));
      let law = SoftTailLaw {
        slope: (x_min - mean) / sigma.powi(2),
      };
      let f = DistributionFunction::from_start_dx_values(x_min, grid[1] - grid[0], y)?
        .with_left_extrap(BoundaryBehavior::Linear(law))?
        .with_right_extrap(BoundaryBehavior::Hard)?;
      Ok(Distribution::Function(f))
    }

    pub(super) fn node_time(posterior: &NodeMap<FullPosterior>, key: GraphNodeKey) -> Option<f64> {
      posterior[&key].time
    }

    pub(super) fn quadratic_neglog(grid: &Array1<f64>, mean: f64, precision: f64) -> Array1<f64> {
      grid.mapv(|t| 0.5 * precision * (t - mean).powi(2))
    }

    pub(super) fn refine_child_of(
      parent_posterior: Distribution<NegLog>,
      message_from_child: Option<Distribution<NegLog>>,
    ) -> Result<Option<f64>, Report> {
      let nwk_parsed = nwk_read(b"(N:1.0)root;".as_slice())?;
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
      let child_key = find_node_key_by_name(&graph, &names, "N").expect("node N not found");
      let mut inputs = ForwardInputs::new(&graph);
      set_time_distribution(&mut inputs, root_key, parent_posterior);
      set_date(&mut inputs, child_key, Distribution::range((2003.0, 2010.0), 0.0));
      set_branch_length_distribution(&graph, &mut inputs, child_key, 1.0);
      inputs.backward.entry(child_key).or_default().to_parent = message_from_child.map(Arc::new);
      let posterior = run_forward_pass(&graph, &inputs, &names)?;
      Ok(node_time(&posterior, child_key))
    }

    pub(super) fn leaf_time_distribution(
      posterior: &NodeMap<FullPosterior>,
      key: GraphNodeKey,
    ) -> Option<Distribution<NegLog>> {
      posterior[&key].distribution.as_deref().cloned()
    }
  }

  use helpers::*;
}

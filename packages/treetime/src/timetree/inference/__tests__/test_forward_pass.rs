#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::pretty_assert_ulps_eq;
  use crate::progress::{LogSink, NoopProgress};
  use crate::test_utils::{RecordingLog, find_node_key_by_name, parent_edge_key, unknown_branches};
  use crate::timetree::inference::forward_pass::{committed_time, propagate_distributions_forward};
  use crate::timetree::inference::result::{BranchLikelihood, NodePosterior, TimeBackward};
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::assert_error;

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
    let nwk_parsed = nwk_read_str("(A:2.5)root;")?;
    let names = nwk_parsed.names();
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

  #[test]
  fn test_forward_pass_refines_uncertain_leaf_date_from_parent() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:1.0)root;")?;
    let names = nwk_parsed.names();
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
    let nwk_parsed = nwk_read_str("(A:1.0)root;")?;
    let names = nwk_parsed.names();
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
    let nwk_parsed = nwk_read_str("(A:1.0)root;")?;
    let names = nwk_parsed.names();
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
    let nwk_parsed = nwk_read_str("(A:1.0)root;")?;
    let names = nwk_parsed.names();
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

    let expected = NodePosterior {
      distribution: Some(Arc::new(given)),
      likely_time: Some(2006.0),
      time: Some(2009.0),
      contradicted: true,
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
      let nwk_parsed = nwk_read_str("(A:1.0)root;")?;
      let names = nwk_parsed.names();
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
    let windowed = Distribution::Function(wide_function.resample_range_dx((2009.4, 2009.7), wide_function.dx())?);

    pretty_assert_ulps_eq!(refine(wide_parent)?, refine(windowed)?, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_forward_pass_refines_a_date_range_narrower_than_the_parent_grid() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:1.0)root;")?;
    let names = nwk_parsed.names();
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
  fn test_forward_pass_missing_backward_output_is_an_internal_error() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut inputs = ForwardInputs::new(&graph);
    inputs.backward.subtree.remove(&leaf_key);

    assert_error!(
      run_forward_pass(&graph, &inputs, &names),
      format!(
        "Backward pass output is missing node {leaf_key}. This is an internal error. Please report it to developers."
      )
    );
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) struct ForwardInputs {
      pub constraints: DateConstraints,
      pub branches: BTreeMap<GraphEdgeKey, BranchLikelihood>,
      pub backward: TimeBackward,
    }

    impl ForwardInputs {
      pub(super) fn new(graph: &Graph) -> Self {
        let backward = TimeBackward {
          subtree: graph.get_nodes().map(|node| (node.key(), None)).collect(),
          messages: graph.get_edges().map(|edge| (edge.key(), None)).collect(),
        };
        Self {
          constraints: DateConstraints::default(),
          branches: unknown_branches(graph),
          backward,
        }
      }
    }

    pub(super) fn run_forward_pass(
      graph: &Graph,
      inputs: &ForwardInputs,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
    ) -> Result<BTreeMap<GraphNodeKey, NodePosterior>, Report> {
      run_forward_pass_with_log(graph, inputs, names, &NoopProgress)
    }

    pub(super) fn run_forward_pass_with_log(
      graph: &Graph,
      inputs: &ForwardInputs,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      log: &dyn LogSink,
    ) -> Result<BTreeMap<GraphNodeKey, NodePosterior>, Report> {
      propagate_distributions_forward(
        graph,
        &inputs.constraints,
        names,
        &inputs.branches,
        &inputs.backward,
        log,
      )
    }

    pub(super) fn set_date(inputs: &mut ForwardInputs, key: GraphNodeKey, dist: Distribution<NegLog>) {
      let dist = Arc::new(dist);
      inputs.constraints.by_node.insert(key, Some(Arc::clone(&dist)));
      inputs.backward.subtree.insert(key, Some(dist));
    }

    pub(super) fn set_time_distribution(inputs: &mut ForwardInputs, key: GraphNodeKey, dist: Distribution<NegLog>) {
      inputs.backward.subtree.insert(key, Some(Arc::new(dist)));
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

    pub(super) fn node_time(posterior: &BTreeMap<GraphNodeKey, NodePosterior>, key: GraphNodeKey) -> Option<f64> {
      posterior[&key].time
    }

    pub(super) fn quadratic_neglog(grid: &Array1<f64>, mean: f64, precision: f64) -> Array1<f64> {
      grid.mapv(|t| 0.5 * precision * (t - mean).powi(2))
    }

    pub(super) fn refine_child_of(
      parent_posterior: Distribution<NegLog>,
      message_from_child: Option<Distribution<NegLog>>,
    ) -> Result<Option<f64>, Report> {
      let nwk_parsed = nwk_read_str("(N:1.0)root;")?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let root_key = find_node_key_by_name(&graph, &names, "root").expect("root not found");
      let child_key = find_node_key_by_name(&graph, &names, "N").expect("node N not found");
      let mut inputs = ForwardInputs::new(&graph);
      set_time_distribution(&mut inputs, root_key, parent_posterior);
      set_date(&mut inputs, child_key, Distribution::range((2003.0, 2010.0), 0.0));
      set_branch_length_distribution(&graph, &mut inputs, child_key, 1.0);
      inputs
        .backward
        .messages
        .insert(parent_edge_key(&graph, child_key), message_from_child.map(Arc::new));
      let posterior = run_forward_pass(&graph, &inputs, &names)?;
      Ok(node_time(&posterior, child_key))
    }

    pub(super) fn leaf_time_distribution(
      posterior: &BTreeMap<GraphNodeKey, NodePosterior>,
      key: GraphNodeKey,
    ) -> Option<Distribution<NegLog>> {
      posterior[&key].distribution.as_deref().cloned()
    }
  }

  use helpers::*;
}

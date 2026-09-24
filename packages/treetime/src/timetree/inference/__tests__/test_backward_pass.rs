#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::coalescent::coalescent::CoalescentModel;
  use crate::pretty_assert_ulps_eq;
  use treetime_utils::assert_error;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::timetree_state::TimetreeState;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use ndarray::array;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::node::GraphNodeKey;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_backward_pass_computes_internal_node_time() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:2.5)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let mut state = TimetreeState::new(&graph);
    set_leaf_time(&mut state, leaf_key, 2013.0);
    set_edge_branch_dist(&graph, &mut state, leaf_key, 2.5);

    let state = run_backward_pass(&graph, &DateConstraints::default(), state, None)?;

    let internal_key = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
    let time_dist = node_time_distribution(&state, internal_key)
      .expect("internal node should have time distribution after backward pass");
    let likely_time = time_dist
      .likely_time()
      .expect("time distribution should have likely_time");

    pretty_assert_ulps_eq!(likely_time, 2010.5, max_ulps = 4);
    assert!(likely_time < 2013.0, "Parent should be older than child");

    Ok(())
  }

  #[test]
  fn test_backward_pass_nan_leaf_likelihood_error_reaches_caller() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:2.5)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let internal_key = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
    let mut state = TimetreeState::new(&graph);
    state.node_mut(leaf_key).time_distribution = Some(Arc::new(Distribution::point(2013.0, f64::NAN)));
    set_edge_branch_dist(&graph, &mut state, leaf_key, 2.5);

    let edge_key = graph
      .get_edges()
      .into_iter()
      .find(|edge| edge.target() == leaf_key)
      .expect("edge above leaf A not found")
      .key();
    let result = propagate_distributions_backward(&graph, &DateConstraints::default(), None, &mut state);

    assert_error!(
      result,
      format!("When sending the time message backward along edge {edge_key}: Cannot normalize a distribution point: its peak negative log-likelihood is NaN")
    );
    assert_eq!(None, node_time_distribution(&state, internal_key));
    Ok(())
  }

  #[test]
  fn test_backward_pass_multiplies_child_messages() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0,B:2.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");

    let mut state = TimetreeState::new(&graph);
    set_leaf_time(&mut state, leaf_a_key, 2015.0);
    set_leaf_time(&mut state, leaf_b_key, 2014.0);
    set_edge_branch_dist(&graph, &mut state, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, &mut state, leaf_b_key, 2.0);

    let state = run_backward_pass(&graph, &DateConstraints::default(), state, None)?;

    let internal_key = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
    let time_dist = node_time_distribution(&state, internal_key).expect("internal node should have time distribution");
    let likely_time = time_dist
      .likely_time()
      .expect("time distribution should have likely_time");

    pretty_assert_ulps_eq!(likely_time, 2012.0, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_backward_pass_preserves_leaf_time_distribution_with_coalescent() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0,B:2.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");

    let date_a = 2015.0;
    let date_b = 2014.0;

    let mut state = TimetreeState::new(&graph);
    set_leaf_time(&mut state, leaf_a_key, date_a);
    set_leaf_time(&mut state, leaf_b_key, date_b);
    set_edge_branch_dist(&graph, &mut state, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, &mut state, leaf_b_key, 2.0);

    let coalescent_model = coalescent_model(0.01)?;

    propagate_distributions_backward(&graph, &DateConstraints::default(), Some(&coalescent_model), &mut state)?;
    propagate_distributions_backward(&graph, &DateConstraints::default(), Some(&coalescent_model), &mut state)?;

    {
      let time_dist = node_time_distribution(&state, leaf_a_key).expect("leaf A should have time distribution");
      let expected = Distribution::point(date_a, 0.0);
      assert_eq!(&expected, time_dist.as_ref());
    }

    {
      let time_dist = node_time_distribution(&state, leaf_b_key).expect("leaf B should have time distribution");
      let expected = Distribution::point(date_b, 0.0);
      assert_eq!(&expected, time_dist.as_ref());
    }

    Ok(())
  }

  #[test]
  fn test_backward_pass_preserves_internal_time_with_strong_coalescent() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0,B:2.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");
    let internal_key = find_node_key_by_name(&graph, &names, "I").expect("internal I not found");

    let mut state = TimetreeState::new(&graph);
    set_leaf_time(&mut state, leaf_a_key, 2015.0);
    set_leaf_time(&mut state, leaf_b_key, 2014.0);
    set_edge_branch_dist(&graph, &mut state, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, &mut state, leaf_b_key, 2.0);

    let coalescent_model = coalescent_model(1e-6)?;

    let state = run_backward_pass(&graph, &DateConstraints::default(), state, Some(&coalescent_model))?;

    let actual = node_time_distribution(&state, internal_key).and_then(|distribution| distribution.likely_time());
    let expected = Some(2012.0);
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_backward_pass_restores_leaf_time_distribution_from_the_date_constraint() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let constraint = Distribution::range((2014.0, 2015.0), 0.0);
    let mut constraints = DateConstraints::default();
    constraints
      .date_constraints
      .insert(leaf_key, Some(Arc::new(constraint.clone())));
    let mut state = TimetreeState::new(&graph);
    {
      let node = state.node_mut(leaf_key);
      node.time_distribution = Some(Arc::new(Distribution::point(2014.2, 0.0)));
    }
    set_edge_branch_dist(&graph, &mut state, leaf_key, 3.0);

    let state = run_backward_pass(&graph, &constraints, state, None)?;

    let actual = node_time_distribution(&state, leaf_key).expect("leaf A should have a time distribution");
    assert_eq!(&constraint, actual.as_ref());

    Ok(())
  }

  #[test]
  fn test_backward_pass_applies_internal_node_date_constraint() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0,B:2.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");
    let internal_key = find_node_key_by_name(&graph, &names, "I").expect("internal I not found");

    let mut constraints = DateConstraints::default();
    set_date_constraint(&mut constraints, leaf_a_key, Distribution::range((2014.0, 2016.0), 0.0));
    set_date_constraint(&mut constraints, leaf_b_key, Distribution::range((2013.0, 2015.0), 0.0));
    set_date_constraint(
      &mut constraints,
      internal_key,
      Distribution::range((2012.0, 2014.0), 0.0),
    );
    let mut state = TimetreeState::new(&graph);
    set_edge_branch_dist(&graph, &mut state, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, &mut state, leaf_b_key, 2.0);

    let state = run_backward_pass(&graph, &constraints, state, None)?;

    let actual = node_time_distribution(&state, internal_key).expect("internal node should have a time distribution");
    let expected = Distribution::range((2012.0, 2013.0), 0.0);
    assert_eq!(&expected, actual.as_ref());

    Ok(())
  }

  #[test]
  fn test_backward_pass_sets_edge_messages() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:2.5)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut state = TimetreeState::new(&graph);
    set_leaf_time(&mut state, leaf_key, 2013.0);
    set_edge_branch_dist(&graph, &mut state, leaf_key, 2.5);

    propagate_distributions_backward(&graph, &DateConstraints::default(), None, &mut state)?;

    for edge in graph.get_edges() {
      let edge_read = edge;
      if edge_read.target() == leaf_key {
        let msg = state
          .edge(edge_read.key())
          .msg_to_parent
          .as_ref()
          .expect("edge should have msg_to_parent after backward pass");
        let msg_time = msg.likely_time().expect("message should have likely_time");
        pretty_assert_ulps_eq!(msg_time, 2010.5, max_ulps = 4);
      }
    }

    Ok(())
  }

  #[test]
  fn test_backward_pass_skips_bad_branch_children() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0,B:2.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");

    let mut state = TimetreeState::new(&graph);
    set_leaf_time(&mut state, leaf_a_key, 2015.0);
    set_leaf_time(&mut state, leaf_b_key, 2014.0);

    state.node_mut(leaf_b_key).bad_branch = true;

    set_edge_branch_dist(&graph, &mut state, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, &mut state, leaf_b_key, 2.0);

    let state = run_backward_pass(&graph, &DateConstraints::default(), state, None)?;

    let internal_key = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
    let time_dist = node_time_distribution(&state, internal_key).expect("internal node should have time distribution");
    let likely_time = time_dist
      .likely_time()
      .expect("time distribution should have likely_time");

    pretty_assert_ulps_eq!(likely_time, 2012.0, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_backward_pass_bad_branch_equivalent_to_removal() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0)I:1.0)root;")?;
    let ref_names = nwk_parsed.names();
    let ref_graph = nwk_parsed.graph;
    let ref_a_key = find_node_key_by_name(&ref_graph, &ref_names, "A").expect("leaf A not found");
    let mut ref_state = TimetreeState::new(&ref_graph);
    set_leaf_time(&mut ref_state, ref_a_key, 2015.0);
    set_edge_branch_dist(&ref_graph, &mut ref_state, ref_a_key, 3.0);
    let ref_state = run_backward_pass(&ref_graph, &DateConstraints::default(), ref_state, None)?;

    let ref_internal_key = find_node_key_by_name(&ref_graph, &ref_names, "I").expect("internal I not found");
    let ref_time = node_time_distribution(&ref_state, ref_internal_key)
      .expect("should have time dist")
      .likely_time()
      .expect("should have likely_time");

    let nwk_parsed = nwk_read_str("((A:3.0,B:2.0)I:1.0)root;")?;
    let test_names = nwk_parsed.names();
    let test_graph = nwk_parsed.graph;
    let test_a_key = find_node_key_by_name(&test_graph, &test_names, "A").expect("leaf A not found");
    let test_b_key = find_node_key_by_name(&test_graph, &test_names, "B").expect("leaf B not found");
    let mut test_state = TimetreeState::new(&test_graph);
    set_leaf_time(&mut test_state, test_a_key, 2015.0);
    set_leaf_time(&mut test_state, test_b_key, 2014.0);
    set_edge_branch_dist(&test_graph, &mut test_state, test_a_key, 3.0);
    set_edge_branch_dist(&test_graph, &mut test_state, test_b_key, 2.0);

    test_state.node_mut(test_b_key).bad_branch = true;

    let test_state = run_backward_pass(&test_graph, &DateConstraints::default(), test_state, None)?;

    let test_internal_key = find_node_key_by_name(&test_graph, &test_names, "I").expect("internal I not found");
    let test_time = node_time_distribution(&test_state, test_internal_key)
      .expect("should have time dist")
      .likely_time()
      .expect("should have likely_time");

    pretty_assert_ulps_eq!(ref_time, test_time, max_ulps = 4);

    Ok(())
  }

  #[test]
  #[ignore = "fold now receives mass-windowed messages; Gaussian-product oracle no longer holds"]
  fn test_backward_pass_sums_function_children_to_gaussian_product() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.0,B:0.0,C:0.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let a = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let b = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");
    let c = find_node_key_by_name(&graph, &names, "C").expect("leaf C not found");

    let x = Array1::linspace(2000.0, 2010.0, 2001);
    let mut state = TimetreeState::new(&graph);
    set_leaf_function(&mut state, a, &x, gaussian_neglog(&x, 2002.0, 1.0))?;
    set_leaf_function(&mut state, b, &x, gaussian_neglog(&x, 2008.0, 1.0))?;
    set_leaf_function(&mut state, c, &x, gaussian_neglog(&x, 2005.0, 2.0))?;
    set_edge_branch_dist(&graph, &mut state, a, 0.0);
    set_edge_branch_dist(&graph, &mut state, b, 0.0);
    set_edge_branch_dist(&graph, &mut state, c, 0.0);

    let state = run_backward_pass(&graph, &DateConstraints::default(), state, None)?;

    let internal = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
    let dist = node_time_distribution(&state, internal).expect("internal node should have a time distribution");

    let grid = dist.t();
    let spacing = grid[1] - grid[0];
    let peak = dist.likely_time().expect("distribution should have a likely_time");
    assert_abs_diff_eq!(peak, 2005.0, epsilon = spacing);

    for t in [2004.0_f64, 2005.0, 2006.0] {
      let expected = 2.0 * (t - 2005.0).powi(2);
      assert_abs_diff_eq!(dist.eval(t)?, expected, epsilon = 1e-4);
    }

    Ok(())
  }

  #[test]
  fn test_backward_pass_fan_out_result_independent_of_child_order() -> Result<(), Report> {
    let x = Array1::linspace(2000.0, 2010.0, 11);
    let ya = gaussian_neglog(&x, 2002.0, 1.0);
    let yb = gaussian_neglog(&x, 2008.0, 1.0);
    let yc = gaussian_neglog(&x, 2005.0, 2.0);

    let fold_in_order = |newick: &str| -> Result<(Array1<f64>, f64), Report> {
      let nwk_parsed = nwk_read_str(newick)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let mut state = TimetreeState::new(&graph);
      for (name, y) in [("A", &ya), ("B", &yb), ("C", &yc)] {
        let key = find_node_key_by_name(&graph, &names, name).expect("leaf not found");
        set_leaf_function(&mut state, key, &x, y.clone())?;
        set_edge_branch_dist(&graph, &mut state, key, 0.0);
      }
      let state = run_backward_pass(&graph, &DateConstraints::default(), state, None)?;
      let internal = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
      let dist = node_time_distribution(&state, internal).expect("internal node should have a time distribution");
      Ok((
        dist.y()?,
        dist.likely_time().expect("distribution should have a likely_time"),
      ))
    };

    let (y_abc, peak_abc) = fold_in_order("((A:0.0,B:0.0,C:0.0)I:1.0)root;")?;
    let (y_cab, peak_cab) = fold_in_order("((C:0.0,A:0.0,B:0.0)I:1.0)root;")?;

    pretty_assert_ulps_eq!(y_abc, y_cab, max_ulps = 8);
    pretty_assert_ulps_eq!(peak_abc, peak_cab, max_ulps = 4);

    Ok(())
  }

  #[test]
  #[ignore = "fold now receives mass-windowed messages; precision-weighted-mean oracle no longer holds"]
  fn test_backward_pass_function_children_peak_at_precision_weighted_mean() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.0,B:0.0,C:0.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let a = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let b = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");
    let c = find_node_key_by_name(&graph, &names, "C").expect("leaf C not found");

    let x = Array1::linspace(2000.0, 2010.0, 11);
    let mut state = TimetreeState::new(&graph);
    set_leaf_function(&mut state, a, &x, gaussian_neglog(&x, 2002.0, 1.0))?;
    set_leaf_function(&mut state, b, &x, gaussian_neglog(&x, 2008.0, 1.0))?;
    set_leaf_function(&mut state, c, &x, gaussian_neglog(&x, 2005.0, 2.0))?;
    set_edge_branch_dist(&graph, &mut state, a, 0.0);
    set_edge_branch_dist(&graph, &mut state, b, 0.0);
    set_edge_branch_dist(&graph, &mut state, c, 0.0);

    let state = run_backward_pass(&graph, &DateConstraints::default(), state, None)?;

    let internal = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
    let dist = node_time_distribution(&state, internal).expect("internal node should have a time distribution");
    let likely_time = dist.likely_time().expect("distribution should have a likely_time");

    let grid = dist.t();
    let spacing = grid[1] - grid[0];
    assert_abs_diff_eq!(likely_time, 2005.0, epsilon = spacing);

    Ok(())
  }

  mod helpers {
    use super::*;
    use treetime_graph::graph::Graph;

    pub(super) fn run_backward_pass(
      graph: &Graph,
      constraints: &DateConstraints,
      mut state: TimetreeState,
      coalescent_model: Option<&CoalescentModel>,
    ) -> Result<TimetreeState, Report> {
      propagate_distributions_backward(graph, constraints, coalescent_model, &mut state)?;
      Ok(state)
    }

    pub(super) fn node_time_distribution(
      state: &TimetreeState,
      key: GraphNodeKey,
    ) -> Option<Arc<Distribution<NegLog>>> {
      state.nodes.get(&key).and_then(|node| node.time_distribution.clone())
    }

    pub(super) fn set_date_constraint(
      constraints: &mut DateConstraints,
      key: GraphNodeKey,
      dist: Distribution<NegLog>,
    ) {
      constraints.date_constraints.insert(key, Some(Arc::new(dist)));
    }

    pub(super) fn set_leaf_time(state: &mut TimetreeState, key: GraphNodeKey, time: f64) {
      state.node_mut(key).time_distribution = Some(Arc::new(Distribution::point(time, 0.0)));
    }

    pub(super) fn set_edge_branch_dist(graph: &Graph, state: &mut TimetreeState, target_key: GraphNodeKey, bl: f64) {
      for edge in graph.get_edges() {
        let edge_read = edge;
        if edge_read.target() == target_key {
          state.edge_mut(edge_read.key()).branch_length_distribution = Some(Arc::new(Distribution::point(bl, 0.0)));
        }
      }
    }

    pub(super) fn gaussian_neglog(x: &Array1<f64>, mean: f64, precision: f64) -> Array1<f64> {
      x.mapv(|t| 0.5 * precision * (t - mean).powi(2))
    }

    pub(super) fn set_leaf_function(
      state: &mut TimetreeState,
      key: GraphNodeKey,
      x: &Array1<f64>,
      y: Array1<f64>,
    ) -> Result<(), Report> {
      let dist = Distribution::function(x.clone(), y)?;
      state.node_mut(key).time_distribution = Some(Arc::new(dist));
      Ok(())
    }

    pub(super) fn coalescent_model(tc: f64) -> Result<CoalescentModel, Report> {
      let lineage_counts = PiecewiseConstantFn::new(array![1900.0, 2100.0], array![1.0, 2.0, 0.0]);
      CoalescentModel::new(&lineage_counts, &Distribution::constant(tc))
    }
  }

  use helpers::*;
}
